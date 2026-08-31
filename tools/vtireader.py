"""
Read Palabos .vti (VTK XML ImageData) files with numpy and the standard library.

WHY THIS EXISTS

The obvious way to read a .vti is to import vtk. That works, and this module
uses it when it is available. But `pip install vtk` is a 100 MB dependency that
often will not build on a login node, and a post-processing script that cannot
run on the machine that produced the data is not much use. So this module falls
back to parsing the file itself, which needs nothing beyond numpy.

WHAT IT HANDLES

VTK XML has more encodings than anyone would like:

    format="ascii"                 numbers as text
    format="binary"                base64, with a length header
    format="appended"              one raw blob at the end of the file,
                                   offsets into it, base64 or raw
    compressor="vtkZLibDataCompressor"    any of the above, zlib-deflated
                                          in blocks
    header_type="UInt32" or "UInt64"      the width of the length prefix

All of those are supported. Anything else raises with a message naming what it
found, rather than returning silently wrong numbers -- which is the failure
mode that matters here, because a post-processor that quietly misreads an array
produces a plot that looks fine and is not.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import base64
import os
import re
import zlib

import numpy as np

VTK_TO_NUMPY = {
    "Int8": np.int8, "UInt8": np.uint8,
    "Int16": np.int16, "UInt16": np.uint16,
    "Int32": np.int32, "UInt32": np.uint32,
    "Int64": np.int64, "UInt64": np.uint64,
    "Float32": np.float32, "Float64": np.float64,
}


class VtiError(Exception):
    pass


def _attr(tag, name, default=None):
    m = re.search(r'\b%s\s*=\s*"([^"]*)"' % name, tag)
    return m.group(1) if m else default


def _decode_blocks(raw, dtype, header_dtype, compressed):
    """Undo VTK's block structure.

    Uncompressed: [nbytes][data].
    Compressed:   [nblocks][blocksize][last_blocksize][c1][c2]...[cn] then the
                  n zlib streams back to back. The header is itself sometimes
                  compressed in newer writers, but never in the ones that
                  produce Palabos output, so it is read plain here.
    """
    hsize = np.dtype(header_dtype).itemsize
    if not compressed:
        nbytes = int(np.frombuffer(raw[:hsize], dtype=header_dtype, count=1)[0])
        return np.frombuffer(raw[hsize:hsize + nbytes], dtype=dtype)

    nblocks = int(np.frombuffer(raw[:hsize], dtype=header_dtype, count=1)[0])
    head = np.frombuffer(raw[:hsize * (3 + nblocks)], dtype=header_dtype,
                         count=3 + nblocks)
    csizes = head[3:]
    off = hsize * (3 + nblocks)
    out = bytearray()
    for c in csizes:
        c = int(c)
        out += zlib.decompress(bytes(raw[off:off + c]))
        off += c
    return np.frombuffer(bytes(out), dtype=dtype)


def _b64_len(nbytes):
    """How many base64 characters encode nbytes bytes, padding included."""
    return ((nbytes + 2) // 3) * 4


def _decode_b64_pair(text, dtype, header_dtype, compressed):
    """Decode a VTK base64 payload.

    VTK has TWO conventions here and getting them confused is the reason a
    naive reader returns garbage on compressed files:

      uncompressed  the 1-word header and the data are encoded TOGETHER as one
                    base64 stream
      compressed    the header is its own base64 stream, followed by a second
                    stream holding the data. It has to be, because the reader
                    cannot know how long the data is until it has read the
                    header.
    """
    text = "".join(text.split())
    hsize = np.dtype(header_dtype).itemsize

    if not compressed:
        return _decode_blocks(base64.b64decode(text), dtype, header_dtype, False)

    # first word of the header tells us how many blocks, hence how long the
    # header is; decode just enough base64 to see it
    first = base64.b64decode(text[:_b64_len(hsize)] + "==")
    nblocks = int(np.frombuffer(first[:hsize], dtype=header_dtype, count=1)[0])
    hbytes = hsize * (3 + nblocks)
    hchars = _b64_len(hbytes)
    header = base64.b64decode(text[:hchars])[:hbytes]
    body = base64.b64decode(text[hchars:])
    return _decode_blocks(header + body, dtype, header_dtype, True)


def _read_inline(text, dtype, fmt, header_dtype, compressed):
    if fmt == "ascii":
        return np.array(text.split(), dtype=dtype)
    if fmt == "binary":
        return _decode_b64_pair(text, dtype, header_dtype, compressed)
    raise VtiError("inline DataArray has format='%s', which is not ascii or "
                   "binary" % fmt)


def read_vti(path):
    """Return (dims, spacing, origin, {name: ndarray}).

    Arrays come back shaped (nz, ny, nx) -- z slowest -- which is VTK's own
    ordering and matches how CompLB3D writes its geometry.
    """
    try:
        return _read_with_vtk(path)
    except ImportError:
        pass
    return _read_manual(path)


def _read_with_vtk(path):
    import vtk                                    # noqa: F401
    from vtk.util.numpy_support import vtk_to_numpy

    r = vtk.vtkXMLImageDataReader()
    r.SetFileName(path)
    r.Update()
    img = r.GetOutput()
    nx, ny, nz = img.GetDimensions()
    arrays = {}
    pd = img.GetPointData()
    for i in range(pd.GetNumberOfArrays()):
        a = pd.GetArray(i)
        if a is None:
            continue
        arrays[a.GetName() or ("array%d" % i)] = \
            vtk_to_numpy(a).reshape((nz, ny, nx))
    return (nx, ny, nz), img.GetSpacing(), img.GetOrigin(), arrays


def _read_manual(path):
    with open(path, "rb") as f:
        blob = f.read()

    # The appended section is raw bytes and must not go through a text decode.
    ap = blob.find(b"<AppendedData")
    if ap >= 0:
        start = blob.find(b"_", ap) + 1
        end = blob.find(b"</AppendedData>", start)
        appended = blob[start:end]
        header_text = blob[:ap].decode("ascii", "replace")
        ap_tag = blob[ap:blob.find(b">", ap)].decode("ascii", "replace")
        ap_encoding = _attr(ap_tag, "encoding", "raw")
    else:
        appended = b""
        header_text = blob.decode("ascii", "replace")
        ap_encoding = "raw"

    vt = re.search(r"<VTKFile[^>]*>", header_text)
    if not vt:
        raise VtiError("%s does not start with a <VTKFile> element; is it a "
                       ".vti?" % os.path.basename(path))
    vtag = vt.group(0)
    if _attr(vtag, "type") != "ImageData":
        raise VtiError("%s has VTKFile type='%s'; this reader handles "
                       "ImageData only" % (os.path.basename(path),
                                           _attr(vtag, "type")))
    compressed = "ZLib" in (_attr(vtag, "compressor") or "")
    header_dtype = np.uint64 if _attr(vtag, "header_type") == "UInt64" else np.uint32
    if _attr(vtag, "byte_order", "LittleEndian") != "LittleEndian":
        raise VtiError("byte_order is BigEndian; not supported")

    im = re.search(r"<ImageData[^>]*>", header_text)
    ext = [int(v) for v in _attr(im.group(0), "WholeExtent", "0 0 0 0 0 0").split()]
    dims = (ext[1] - ext[0] + 1, ext[3] - ext[2] + 1, ext[5] - ext[4] + 1)
    spacing = tuple(float(v) for v in _attr(im.group(0), "Spacing", "1 1 1").split())
    origin = tuple(float(v) for v in _attr(im.group(0), "Origin", "0 0 0").split())

    arrays = {}
    for m in re.finditer(r"<DataArray([^>]*?)(/>|>(.*?)</DataArray>)",
                         header_text, re.S):
        attrs, closing, body = m.group(1), m.group(2), m.group(3)
        name = _attr(attrs, "Name")
        vtype = _attr(attrs, "type")
        fmt = _attr(attrs, "format", "ascii")
        if name is None or vtype not in VTK_TO_NUMPY:
            continue
        dtype = VTK_TO_NUMPY[vtype]

        if fmt == "appended":
            off = int(_attr(attrs, "offset", "0"))
            raw = appended[off:]
            if ap_encoding == "base64":
                data = _decode_b64_pair(raw.decode("ascii", "replace"),
                                        dtype, header_dtype, compressed)
            else:
                data = _decode_blocks(raw, dtype, header_dtype, compressed)
        else:
            if closing == "/>":
                continue
            data = _read_inline(body, dtype, fmt, header_dtype, compressed)

        ncomp = int(_attr(attrs, "NumberOfComponents", "1"))
        want = dims[0] * dims[1] * dims[2] * ncomp
        if data.size < want:
            raise VtiError("array '%s' has %d values but the extent needs %d"
                           % (name, data.size, want))
        data = data[:want]
        shape = (dims[2], dims[1], dims[0]) if ncomp == 1 else \
                (dims[2], dims[1], dims[0], ncomp)
        arrays[name] = data.reshape(shape)

    if not arrays:
        raise VtiError("no readable DataArray found in %s"
                       % os.path.basename(path))
    return dims, spacing, origin, arrays


# ---------------------------------------------------------------------------
def write_vti(path, array, spacing=(1.0, 1.0, 1.0), name="value"):
    """Write a (nz, ny, nx) array as an ASCII .vti.

    Used by geometry.py to make a pore structure viewable in ParaView before a
    single simulation step has been run. ASCII is slower and larger than
    binary, but a geometry is written once and read by eye, and a format
    anybody can open in a text editor is worth more here than compactness.
    """
    a = np.asarray(array)
    nz, ny, nx = a.shape
    dt = "Float64" if a.dtype.kind == "f" else "Int32"
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="ImageData" version="0.1" '
                'byte_order="LittleEndian">\n')
        f.write('  <ImageData WholeExtent="0 %d 0 %d 0 %d" Origin="0 0 0" '
                'Spacing="%g %g %g">\n' % (nx - 1, ny - 1, nz - 1, *spacing))
        f.write('    <Piece Extent="0 %d 0 %d 0 %d">\n' % (nx - 1, ny - 1, nz - 1))
        f.write('      <PointData Scalars="%s">\n' % name)
        f.write('        <DataArray type="%s" Name="%s" format="ascii">\n' % (dt, name))
        flat = a.ravel()
        for i in range(0, flat.size, 20):
            f.write("          " + " ".join(str(v) for v in flat[i:i + 20]) + "\n")
        f.write('        </DataArray>\n')
        f.write('      </PointData>\n    </Piece>\n  </ImageData>\n</VTKFile>\n')
