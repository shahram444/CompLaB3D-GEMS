/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_geometry.hh  --  BUILD OR IMPORT THE PORE SPACE FROM THE XML
 * ================================================================================================
 *
 *  WHAT THIS REPLACES
 *
 *  makeGeometry.py and `geometry.py create`. Instead of generating a .dat file beforehand and
 *  pointing at it, the domain block can describe the pore space directly:
 *
 *      <domain>
 *          <generate>random</generate>       <!-- or spheres, cylinders, fracture, layered, channel -->
 *          <porosity>0.35</porosity>
 *          <seed>1</seed>
 *          <walls>y</walls>
 *          <!-- or, as before:  <filename>geometry.dat</filename> -->
 *      </domain>
 *
 *  A raw binary volume can be read directly too, which covers segmented CT data exported from
 *  ImageJ or similar:
 *
 *          <import_raw>scan.raw</import_raw>
 *          <raw_dtype>uint8</raw_dtype>
 *          <threshold>128</threshold>
 *
 *  PNG and TIFF stacks are NOT read here and will not be: that needs an image codec, and a
 *  half-working one is worse than none. tools/geometry.py handles those, and writes a .dat or a
 *  raw volume this can read.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE INSPECTION IS THE POINT
 *
 *  Whether the geometry was generated, imported or read from a file, inspect() runs on it before
 *  the first time step and reports porosity, whether the pore space PERCOLATES from inlet to
 *  outlet, and how much of it is ISOLATED. A sealed domain makes the pressure boundary
 *  ill-posed, and a disconnected pocket holds solute that can never leave -- both quietly wreck
 *  domain-averaged results, and both are cheap to detect and expensive to discover later.
 *
 *  ------------------------------------------------------------------------------------------------
 *  ORDERING. The .dat format is x fastest, then y, then z, with no header, so nothing in the file
 *  records its own shape. That contract is implemented once, here, and the same ordering is used
 *  by tools/geometry.py -- verified against it rather than assumed.
 * ================================================================================================
 */
#ifndef COMPLAB3D_GEOMETRY_HH
#define COMPLAB3D_GEOMETRY_HH

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace complab_geom {

const int SOLID = 0, WALL = 1, PORE = 2;

struct Grid {
    int nx, ny, nz;
    std::vector<int> v;                       // index(x,y,z) = x + nx*(y + ny*z)
    Grid() : nx(0), ny(0), nz(0) {}
    Grid(int X, int Y, int Z, int fill = PORE) : nx(X), ny(Y), nz(Z),
        v((size_t) X * (size_t) Y * (size_t) Z, fill) {}
    int &at(int x, int y, int z) { return v[(size_t) x + (size_t) nx * ((size_t) y + (size_t) ny * (size_t) z)]; }
    int at(int x, int y, int z) const { return v[(size_t) x + (size_t) nx * ((size_t) y + (size_t) ny * (size_t) z)]; }
    size_t size() const { return v.size(); }
};

/* Deterministic RNG, for the same reason as in the surrogate header: <random>'s distributions are
 * not portable across standard libraries, and a geometry that changes with the compiler is not a
 * reproducible input. */
struct Rng {
    unsigned long long s;
    explicit Rng(unsigned long long seed = 1) : s(seed ? seed : 1) {}
    unsigned long long next() { s ^= s >> 12; s ^= s << 25; s ^= s >> 27; return s * 2685821657736338717ULL; }
    double uniform() { return (double) (next() >> 11) / 9007199254740992.0; }
    int below(int n) { return n > 0 ? (int) (uniform() * (double) n) : 0; }
};

/* ------------------------------------------------------------------------------------------------
 *  Walls. Never on x by default: x carries the pressure boundary, and walling it off leaves a
 *  domain with no inlet.
 * ------------------------------------------------------------------------------------------------ */
inline void addWalls(Grid &g, const std::string &axes)
{
    const bool wx = axes.find('x') != std::string::npos || axes == "all";
    const bool wy = axes.find('y') != std::string::npos || axes == "all";
    const bool wz = axes.find('z') != std::string::npos || axes == "all";
    for (int z = 0; z < g.nz; ++z)
        for (int y = 0; y < g.ny; ++y)
            for (int x = 0; x < g.nx; ++x) {
                const bool onX = (x == 0 || x == g.nx - 1);
                const bool onY = (y == 0 || y == g.ny - 1);
                const bool onZ = (z == 0 || z == g.nz - 1);
                if ((wx && onX) || (wy && onY) || (wz && onZ)) g.at(x, y, z) = WALL;
            }
}

/* ================================================================================================
 *  Generators
 * ================================================================================================ */
struct GenOptions {
    std::string kind;
    double porosity, radius, aperture, roughness;
    int count, layers;
    unsigned long long seed;
    std::string walls;
    GenOptions() : kind("channel"), porosity(0.4), radius(6), aperture(6), roughness(2),
                   count(20), layers(4), seed(1), walls("y") {}
};

inline void smooth(std::vector<double> &f, int nx, int ny, int nz, int passes)
{
    std::vector<double> t(f.size());
    for (int p = 0; p < passes; ++p) {
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
                for (int x = 0; x < nx; ++x) {
                    double s = 0.0; int n = 0;
                    for (int dz = -1; dz <= 1; ++dz)
                        for (int dy = -1; dy <= 1; ++dy)
                            for (int dx = -1; dx <= 1; ++dx) {
                                const int X = (x + dx + nx) % nx, Y = (y + dy + ny) % ny, Z = (z + dz + nz) % nz;
                                s += f[(size_t) X + (size_t) nx * ((size_t) Y + (size_t) ny * (size_t) Z)];
                                ++n;
                            }
                    t[(size_t) x + (size_t) nx * ((size_t) y + (size_t) ny * (size_t) z)] = s / (double) n;
                }
        f.swap(t);
    }
}

inline Grid generate(int nx, int ny, int nz, const GenOptions &o, std::string *err = 0)
{
    Grid g(nx, ny, nz, PORE);
    Rng rng(o.seed);

    if (o.kind == "channel") {
        /* nothing to carve */
    }
    else if (o.kind == "spheres" || o.kind == "cylinders") {
        const bool cyl = (o.kind == "cylinders");
        for (int k = 0; k < o.count; ++k) {
            const int cx = rng.below(nx), cy = rng.below(ny);
            const int cz = cyl ? 0 : rng.below(nz);
            for (int z = 0; z < nz; ++z)
                for (int y = 0; y < ny; ++y)
                    for (int x = 0; x < nx; ++x) {
                        const double dx = x - cx, dy = y - cy, dz = cyl ? 0.0 : (double) (z - cz);
                        if (dx * dx + dy * dy + dz * dz <= o.radius * o.radius) g.at(x, y, z) = SOLID;
                    }
        }
    }
    else if (o.kind == "random") {
        /* Threshold a smoothed random field at the quantile that gives the requested porosity, so
         * the porosity comes out exact rather than approximately right. */
        std::vector<double> f(g.size());
        for (size_t i = 0; i < f.size(); ++i) f[i] = rng.uniform();
        smooth(f, nx, ny, nz, 2);
        std::vector<double> sorted(f);
        std::sort(sorted.begin(), sorted.end());
        const size_t cut = (size_t) ((1.0 - o.porosity) * (double) sorted.size());
        const double thr = sorted[cut < sorted.size() ? cut : sorted.size() - 1];
        for (size_t i = 0; i < f.size(); ++i) g.v[i] = (f[i] >= thr) ? PORE : SOLID;
    }
    else if (o.kind == "fracture") {
        for (size_t i = 0; i < g.size(); ++i) g.v[i] = SOLID;
        const double mid = 0.5 * ny;
        for (int z = 0; z < nz; ++z) {
            std::vector<double> wob((size_t) nx, 0.0);
            for (int k = 1; k <= 3; ++k) {
                const double ph = rng.uniform() * 6.283185307179586;
                for (int x = 0; x < nx; ++x)
                    wob[(size_t) x] += (o.roughness / k) * std::sin(6.283185307179586 * k * x / nx + ph);
            }
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                    if (std::fabs(y - (mid + wob[(size_t) x])) <= 0.5 * o.aperture) g.at(x, y, z) = PORE;
        }
    }
    else if (o.kind == "layered") {
        for (int i = 0; i < o.layers; ++i) {
            const int y0 = (int) ((double) i * ny / o.layers);
            const int y1 = (int) ((double) (i + 1) * ny / o.layers);
            const double p = o.porosity * ((i % 2 == 0) ? 1.4 : 0.5);
            for (int z = 0; z < nz; ++z)
                for (int y = y0; y < y1; ++y)
                    for (int x = 0; x < nx; ++x)
                        if (rng.uniform() > (p < 1.0 ? p : 1.0)) g.at(x, y, z) = SOLID;
        }
    }
    else {
        if (err) *err = "<generate>" + o.kind + "</generate> is not one of: channel, spheres, "
                        "cylinders, random, fracture, layered";
        return Grid();
    }

    if (o.walls != "none" && !o.walls.empty()) addWalls(g, o.walls);
    return g;
}

/* ================================================================================================
 *  Import and export
 * ================================================================================================ */
inline bool importRaw(Grid &g, const std::string &path, int nx, int ny, int nz,
                      const std::string &dtype, double threshold, bool invert,
                      const std::string &walls, std::string *err = 0)
{
    std::FILE *f = std::fopen(path.c_str(), "rb");
    if (!f) { if (err) *err = "cannot open raw volume '" + path + "'"; return false; }

    size_t esize = 1;
    if (dtype == "uint8" || dtype == "int8") esize = 1;
    else if (dtype == "uint16" || dtype == "int16") esize = 2;
    else if (dtype == "float32" || dtype == "uint32" || dtype == "int32") esize = 4;
    else if (dtype == "float64") esize = 8;
    else { std::fclose(f); if (err) *err = "raw_dtype '" + dtype + "' is not one of uint8, int8, "
                                           "uint16, int16, uint32, int32, float32, float64"; return false; }

    const size_t n = (size_t) nx * (size_t) ny * (size_t) nz;
    std::vector<unsigned char> buf(n * esize);
    const size_t got = std::fread(&buf[0], 1, buf.size(), f);
    std::fclose(f);
    if (got != buf.size()) {
        if (err) {
            char b[256];
            std::snprintf(b, sizeof(b), "'%s' holds %lu bytes but nx*ny*nz*%lu = %lu",
                          path.c_str(), (unsigned long) got, (unsigned long) esize,
                          (unsigned long) buf.size());
            *err = b;
        }
        return false;
    }

    g = Grid(nx, ny, nz, PORE);
    for (size_t i = 0; i < n; ++i) {
        double val = 0.0;
        const unsigned char *p = &buf[i * esize];
        if (dtype == "uint8")        val = (double) *p;
        else if (dtype == "int8")    val = (double) *(const signed char *) p;
        else if (dtype == "uint16")  { unsigned short t; std::memcpy(&t, p, 2); val = (double) t; }
        else if (dtype == "int16")   { short t;          std::memcpy(&t, p, 2); val = (double) t; }
        else if (dtype == "uint32")  { unsigned int t;   std::memcpy(&t, p, 4); val = (double) t; }
        else if (dtype == "int32")   { int t;            std::memcpy(&t, p, 4); val = (double) t; }
        else if (dtype == "float32") { float t;          std::memcpy(&t, p, 4); val = (double) t; }
        else                         { double t;         std::memcpy(&t, p, 8); val = t; }
        bool pore = (val >= threshold);
        if (invert) pore = !pore;
        g.v[i] = pore ? PORE : SOLID;
    }
    if (walls != "none" && !walls.empty()) addWalls(g, walls);
    return true;
}

inline bool writeDat(const Grid &g, const std::string &path, std::string *err = 0)
{
    std::FILE *f = std::fopen(path.c_str(), "w");
    if (!f) { if (err) *err = "cannot write '" + path + "'"; return false; }
    for (int z = 0; z < g.nz; ++z)
        for (int y = 0; y < g.ny; ++y)
            for (int x = 0; x < g.nx; ++x)
                std::fprintf(f, "%d\n", g.at(x, y, z));
    std::fclose(f);
    return true;
}

/* ================================================================================================
 *  Inspection
 * ================================================================================================ */
struct Inspection {
    double porosity;
    long openVoxels, isolatedVoxels;
    int clusters, spanningClusters;
    bool percolates, quasi2D;
    std::vector<std::pair<int, long> > counts;      // material number, voxel count
    Inspection() : porosity(0), openVoxels(0), isolatedVoxels(0), clusters(0),
                   spanningClusters(0), percolates(false), quasi2D(false) {}
};

/* Flood fill with an explicit stack. Recursion would overflow on a real CT volume: a 500^3 scan
 * has 125 million voxels and a single connected pore network through most of them. */
inline Inspection inspect(const Grid &g)
{
    Inspection R;
    const size_t n = g.size();

    std::vector<int> mat;
    std::vector<long> cnt;
    for (size_t i = 0; i < n; ++i) {
        const int m = g.v[i];
        size_t k = 0;
        while (k < mat.size() && mat[k] != m) ++k;
        if (k == mat.size()) { mat.push_back(m); cnt.push_back(0); }
        ++cnt[k];
    }
    for (size_t k = 0; k < mat.size(); ++k) R.counts.push_back(std::make_pair(mat[k], cnt[k]));

    std::vector<char> open(n, 0);
    for (size_t i = 0; i < n; ++i) if (g.v[i] >= PORE) { open[i] = 1; ++R.openVoxels; }
    R.porosity = (double) R.openVoxels / (double) n;

    std::vector<int> lab(n, 0);
    int cur = 0;
    std::vector<size_t> stack;
    for (size_t start = 0; start < n; ++start) {
        if (!open[start] || lab[start]) continue;
        ++cur;
        stack.clear();
        stack.push_back(start);
        lab[start] = cur;
        while (!stack.empty()) {
            const size_t i = stack.back(); stack.pop_back();
            const int x = (int) (i % (size_t) g.nx);
            const int y = (int) ((i / (size_t) g.nx) % (size_t) g.ny);
            const int z = (int) (i / ((size_t) g.nx * (size_t) g.ny));
            const int dx[6] = {1, -1, 0, 0, 0, 0}, dy[6] = {0, 0, 1, -1, 0, 0}, dz[6] = {0, 0, 0, 0, 1, -1};
            for (int d = 0; d < 6; ++d) {
                const int X = x + dx[d], Y = y + dy[d], Z = z + dz[d];
                if (X < 0 || Y < 0 || Z < 0 || X >= g.nx || Y >= g.ny || Z >= g.nz) continue;
                const size_t j = (size_t) X + (size_t) g.nx * ((size_t) Y + (size_t) g.ny * (size_t) Z);
                if (open[j] && !lab[j]) { lab[j] = cur; stack.push_back(j); }
            }
        }
    }
    R.clusters = cur;

    std::vector<char> atInlet((size_t) cur + 1, 0), atOutlet((size_t) cur + 1, 0);
    for (int z = 0; z < g.nz; ++z)
        for (int y = 0; y < g.ny; ++y) {
            const size_t a = (size_t) 0 + (size_t) g.nx * ((size_t) y + (size_t) g.ny * (size_t) z);
            const size_t b = (size_t) (g.nx - 1) + (size_t) g.nx * ((size_t) y + (size_t) g.ny * (size_t) z);
            if (lab[a]) atInlet[(size_t) lab[a]] = 1;
            if (lab[b]) atOutlet[(size_t) lab[b]] = 1;
        }
    std::vector<char> spanning((size_t) cur + 1, 0);
    for (int c = 1; c <= cur; ++c)
        if (atInlet[(size_t) c] && atOutlet[(size_t) c]) { spanning[(size_t) c] = 1; ++R.spanningClusters; }
    R.percolates = R.spanningClusters > 0;

    for (size_t i = 0; i < n; ++i)
        if (open[i] && !spanning[(size_t) lab[i]]) ++R.isolatedVoxels;

    R.quasi2D = true;
    for (int z = 1; z < g.nz && R.quasi2D; ++z)
        for (int y = 0; y < g.ny && R.quasi2D; ++y)
            for (int x = 0; x < g.nx; ++x)
                if (g.at(x, y, z) != g.at(x, y, 0)) { R.quasi2D = false; break; }

    return R;
}

inline std::string inspectionText(const Grid &g, const Inspection &R)
{
    char b[512];
    std::string s;
    std::snprintf(b, sizeof(b), "  [GEOM] %d x %d x %d = %lu voxels\n",
                  g.nx, g.ny, g.nz, (unsigned long) g.size()); s += b;
    for (size_t k = 0; k < R.counts.size(); ++k) {
        const char *nm = R.counts[k].first == SOLID ? "solid" :
                         R.counts[k].first == WALL ? "bounce_back" :
                         R.counts[k].first == PORE ? "pore" : "microbe seed";
        std::snprintf(b, sizeof(b), "  [GEOM]   material %-3d %-13s %10ld  %6.2f%%\n",
                      R.counts[k].first, nm, R.counts[k].second,
                      100.0 * (double) R.counts[k].second / (double) g.size());
        s += b;
    }
    std::snprintf(b, sizeof(b), "  [GEOM] porosity %.4f\n", R.porosity); s += b;

    if (R.percolates) {
        std::snprintf(b, sizeof(b), "  [GEOM] percolates along x: yes, %d spanning cluster(s)\n",
                      R.spanningClusters);
        s += b;
    } else {
        s += "  [GEOM] percolates along x: NO\n";
        s += "  [GEOM] WARNING: there is no open path from the inlet face to the outlet face.\n";
        s += "  [GEOM]          A pressure drop across a sealed domain is not a well-posed\n";
        s += "  [GEOM]          problem and the flow solver will not converge to anything useful.\n";
    }
    if (R.isolatedVoxels > 0) {
        const double frac = 100.0 * (double) R.isolatedVoxels / (double) (R.openVoxels ? R.openVoxels : 1);
        std::snprintf(b, sizeof(b), "  [GEOM] isolated pore space: %ld voxel(s), %.1f%% of the open space\n",
                      R.isolatedVoxels, frac);
        s += b;
        if (frac > 5.0) {
            s += "  [GEOM] WARNING: solute in there can never reach the inlet, so it sits at its\n";
            s += "  [GEOM]          initial value forever and drags every domain-averaged number.\n";
        }
    }
    if (g.nz > 1 && R.quasi2D)
        s += "  [GEOM] every z-plane is identical, so the answer must not depend on z\n";
    return s;
}

}  // namespace complab_geom

#endif  // COMPLAB3D_GEOMETRY_HH
