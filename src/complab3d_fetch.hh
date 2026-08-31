/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_fetch.hh  --  GET A GENOME-SCALE MODEL WITHOUT LEAVING THE XML
 * ================================================================================================
 *
 *      <model_source>bigg:iJO1366</model_source>
 *      <model_cache>models</model_cache>
 *      <allow_download>true</allow_download>
 *
 *  THREE PLACES ARE TRIED, IN ORDER, AND ONLY THE LAST TOUCHES THE NETWORK:
 *
 *      1. bundled     models/iJO1366.xml.gz, shipped with the code. Decompressed once into the
 *                     cache. Several common models are bundled; see models/NOTICE.md.
 *      2. cache       models/iJO1366.xml, from a previous run or a previous download.
 *      3. download    only if <allow_download> is true AND the machine has network.
 *
 *  So on a cluster the normal case is that nothing is ever downloaded at all.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE FILE IS CHECKED AGAINST models/manifest.txt, AND THAT IS THE POINT
 *
 *  A model file being present is not the same as it being the RIGHT file. Two failures matter:
 *
 *    * a corrupted or truncated download, caught by the checksum;
 *    * a DIFFERENT REVISION of the same model, caught by the metabolite and reaction counts.
 *
 *  The second is the dangerous one. BiGG revises models. <exchange_reaction_indices> in CompLaB.xml
 *  is POSITIONAL, so if a revision inserts one reaction, every index after it now points somewhere
 *  else -- and the run proceeds happily, on the wrong exchanges, with no error at all. The manifest
 *  records the dimensions each model is expected to have, so that becomes a warning instead of a
 *  silently wrong answer.
 *
 *  The real fix for that class of bug is to name exchanges rather than number them; see
 *  resolveExchangeNames() in complab3d_sbml.hh.
 *
 *  ------------------------------------------------------------------------------------------------
 *  READ THIS BEFORE RELYING ON IT
 *
 *  COMPUTE NODES USUALLY HAVE NO OUTBOUND NETWORK. On Tahoma and most clusters, the login node can
 *  reach the internet and the compute nodes cannot. A job that tries to download at start-up will
 *  hang or fail, on every rank at once.
 *
 *  So the intended pattern is: fetch once, on the login node, into a cache directory, then run
 *  offline. The code is built for that --
 *
 *    * only rank 0 ever downloads; the others wait for the file to appear,
 *    * a cached file is used without any network access at all,
 *    * <allow_download>false</allow_download> makes a missing model a clean start-up error
 *      naming the URL, instead of a hang.
 *
 *  Setting allow_download true in a batch job with no network is the one way to use this badly,
 *  and the message says so.
 *
 *  ------------------------------------------------------------------------------------------------
 *  HOW IT DOWNLOADS
 *
 *  By running curl or wget, not by linking libcurl. That is a deliberate trade: linking libcurl
 *  would add a hard dependency to every build, including the overwhelming majority that never
 *  download anything, in exchange for a feature used once per model. Shelling out costs nothing
 *  when unused, and both tools are present on every cluster.
 *
 *  The download goes to a temporary name and is renamed only after it succeeds and the file looks
 *  like SBML. A half-downloaded model that is silently treated as complete would be far worse than
 *  no download at all.
 * ================================================================================================
 */
#ifndef COMPLAB3D_FETCH_HH
#define COMPLAB3D_FETCH_HH

#include <cstdio>
#include <cstdlib>
#include <string>

namespace complab_fetch {

inline bool fileExistsRaw(const std::string &p);

struct ManifestEntry;

struct Result {
    std::string path;          // where the model now is
    bool downloaded;           // true if it was fetched over the network this run
    bool fromBundle;           // true if it came from the models/ shipped with the code
    bool ok;
    bool dimensionMismatch;    // the file is a DIFFERENT REVISION than the manifest expects
    std::string message;
    Result();
};

/* ------------------------------------------------------------------------------------------------
 *  FNV-1a, 64 bit.
 *
 *  A corruption and wrong-file check, not a security measure -- so a short, dependency-free hash is
 *  the right tool. Linking a cryptographic hash for this would be a dependency bought for nothing.
 * ------------------------------------------------------------------------------------------------ */
inline unsigned long long fnv1a64(const unsigned char *p, size_t n)
{
    unsigned long long h = 14695981039346656037ULL;
    for (size_t i = 0; i < n; ++i) {
        h ^= (unsigned long long) p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

inline bool hashFile(const std::string &path, unsigned long long &out, size_t &bytes)
{
    std::FILE *f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    unsigned long long h = 14695981039346656037ULL;
    unsigned char buf[65536];
    size_t total = 0, n;
    while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0) {
        total += n;
        for (size_t i = 0; i < n; ++i) { h ^= (unsigned long long) buf[i]; h *= 1099511628211ULL; }
    }
    std::fclose(f);
    out = h;
    bytes = total;
    return true;
}

/* ------------------------------------------------------------------------------------------------
 *  One manifest row.
 * ------------------------------------------------------------------------------------------------ */
struct ManifestEntry {
    std::string name, objective, url;
    int rawSpecies, rawReactions;    // <species / <reaction tags in the file: a pre-parse check
    int nmet, nrxn;                  // what the PARSED model has: the exact check
    unsigned long long hash;
    size_t bytes;
    bool found;
    ManifestEntry() : rawSpecies(0), rawReactions(0), nmet(0), nrxn(0),
                      hash(0), bytes(0), found(false) {}
};

inline Result::Result() : downloaded(false), fromBundle(false), ok(false),
                          dimensionMismatch(false) {}

inline ManifestEntry readManifest(const std::string &dir, const std::string &name)
{
    ManifestEntry e;
    const std::string path = dir + "/manifest.txt";
    std::FILE *f = std::fopen(path.c_str(), "r");
    if (!f) return e;
    char line[1024];
    while (std::fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char nm[256] = {0}, ob[256] = {0}, ur[512] = {0};
        int rs = 0, rr = 0, nmet = 0, nrxn = 0;
        unsigned long long h = 0, by = 0;
        if (std::sscanf(line, "%255s %d %d %d %d %255s %llu %llu %511s",
                        nm, &rs, &rr, &nmet, &nrxn, ob, &h, &by, ur) >= 5) {
            if (name == nm) {
                e.name = nm;
                e.rawSpecies = rs; e.rawReactions = rr;
                e.nmet = nmet; e.nrxn = nrxn; e.objective = ob;
                e.hash = h; e.bytes = (size_t) by; e.url = ur; e.found = true;
                break;
            }
        }
    }
    std::fclose(f);
    return e;
}

/* Decompress bundled/<name>.xml.gz into cache/<name>.xml.
 *
 * By shelling out to gunzip rather than linking zlib, for the same reason the download shells out
 * to curl: a build that never uses this should not carry the dependency. gunzip is on every
 * cluster; if it somehow is not, the message says exactly what to run by hand. */
inline bool gunzipTo(const std::string &gz, const std::string &out, std::string &msg)
{
    const std::string tmp = out + ".part";
    std::string cmd = "gunzip -c '" + gz + "' > '" + tmp + "' 2>/dev/null";
    int rc = std::system(cmd.c_str());
    if (rc != 0) {
        cmd = "zcat '" + gz + "' > '" + tmp + "' 2>/dev/null";
        rc = std::system(cmd.c_str());
    }
    if (rc != 0 || !fileExistsRaw(tmp)) {
        std::remove(tmp.c_str());
        msg += "  [MODEL] could not decompress " + gz + ".\n"
               "  [MODEL] Neither gunzip nor zcat worked. Unpack it by hand:\n"
               "  [MODEL]     gunzip -c " + gz + " > " + out + "\n";
        return false;
    }
    if (std::rename(tmp.c_str(), out.c_str()) != 0) {
        std::remove(tmp.c_str());
        msg += "  [MODEL] unpacked " + gz + " but could not write " + out + "\n";
        return false;
    }
    return true;
}

inline bool fileExistsRaw(const std::string &p)
{
    std::FILE *f = std::fopen(p.c_str(), "rb");
    if (!f) return false;
    std::fclose(f);
    return true;
}

/* Does this look like an SBML model, rather than an HTML error page? Servers answer a missing
 * model with a 200 and a page of HTML surprisingly often, and a 40 kB "Not Found" page saved as
 * iAF987.xml would fail much later with a confusing parse error. */
inline bool fileExists(const std::string &p) { return fileExistsRaw(p); }

inline bool looksLikeSbml(const std::string &p)
{
    std::FILE *f = std::fopen(p.c_str(), "rb");
    if (!f) return false;
    char buf[4096];
    const size_t n = std::fread(buf, 1, sizeof(buf) - 1, f);
    std::fclose(f);
    buf[n] = '\0';
    const std::string head(buf, n);
    return head.find("<sbml") != std::string::npos || head.find("<Metabolic_Model") != std::string::npos;
}

/* Where a source name lives. Only BiGG is wired up; the shape is here for others. */
inline bool resolveUrl(const std::string &source, std::string &url, std::string &id,
                       std::string &err)
{
    const size_t colon = source.find(':');
    if (colon == std::string::npos) {
        err = "<model_source>" + source + "</model_source> has no provider prefix. "
              "Write it as bigg:iAF987.";
        return false;
    }
    const std::string prov = source.substr(0, colon);
    id = source.substr(colon + 1);
    if (id.empty()) { err = "<model_source>" + source + "</model_source> names no model"; return false; }

    if (prov == "bigg") {
        url = "http://bigg.ucsd.edu/static/models/" + id + ".xml";
        return true;
    }
    if (prov == "url") {              // an explicit URL, for anything else
        url = source.substr(colon + 1);
        const size_t slash = url.find_last_of('/');
        id = (slash == std::string::npos) ? "model" : url.substr(slash + 1);
        const size_t dot = id.find_last_of('.');
        if (dot != std::string::npos) id = id.substr(0, dot);
        return true;
    }
    err = "provider '" + prov + "' is not known. Use bigg: or url:.";
    return false;
}

/* ------------------------------------------------------------------------------------------------
 *  isMaster: only rank 0 should download. Other ranks call this too, with isMaster false, and get
 *  the resolved path without touching the network.
 * ------------------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------------------
 *  Is this the file the manifest says it is?
 *
 *  Two questions, and the second is the one that matters.
 *
 *    checksum   catches a truncated or corrupted copy. A mismatch here is not necessarily fatal --
 *               BiGG reformats files without changing the science -- so it is a note, not an error.
 *
 *    dimensions caught by counting <species> and <reaction> in the file. A DIFFERENT COUNT means a
 *               different revision of the model, and that is dangerous in a specific way:
 *               <exchange_reaction_indices> is positional. Insert one reaction upstream and every
 *               index after it silently points at the wrong thing, with no error anywhere. This is
 *               the check that turns that into a warning.
 * ------------------------------------------------------------------------------------------------ */
inline void countSbmlEntities(const std::string &path, int &nspecies, int &nreactions)
{
    /* Read the whole file rather than streaming it.
     *
     * The streaming version of this counted a tag twice whenever it straddled a buffer boundary,
     * which showed up as a phantom "different revision" warning on iJO1366. A model file is at most
     * a few tens of megabytes and the SBML parser is about to build a DOM of the same file anyway,
     * so chunking bought nothing and cost correctness.
     *
     * These are RAW TAG COUNTS. They are deliberately not the same as the model's metabolite count,
     * because boundaryCondition species are tags here and are not metabolites in the model. The
     * manifest stores raw counts to match, and the exact check happens after parsing. */
    nspecies = nreactions = 0;
    std::FILE *f = std::fopen(path.c_str(), "rb");
    if (!f) return;
    std::fseek(f, 0, SEEK_END);
    const long len = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    if (len <= 0) { std::fclose(f); return; }

    std::string all((size_t) len, '\0');
    const size_t got = std::fread(&all[0], 1, (size_t) len, f);
    std::fclose(f);
    all.resize(got);

    size_t p = 0;
    while ((p = all.find("<species ", p)) != std::string::npos) { ++nspecies; p += 9; }
    p = 0;
    while ((p = all.find("<reaction ", p)) != std::string::npos) { ++nreactions; p += 10; }
}

/* The exact check, once the model has actually been parsed. Cheap, and it compares like with like:
 * metabolites the LP will contain against metabolites the manifest recorded. */
inline std::string verifyParsedDimensions(const std::string &name, int nmet, int nrxn,
                                          int wantMet, int wantRxn)
{
    if (wantMet <= 0 || wantRxn <= 0) return "";
    if (nmet == wantMet && nrxn == wantRxn) return "";
    char b[512];
    std::snprintf(b, sizeof(b),
        "  [MODEL] WARNING: %s parsed as %d metabolites x %d reactions;\n"
        "  [MODEL]   the manifest records %d x %d. This is a DIFFERENT REVISION.\n"
        "  [MODEL]   <exchange_reaction_indices> is positional, so your indices may now point at\n"
        "  [MODEL]   different reactions and nothing would complain. Use\n"
        "  [MODEL]   <exchange_reaction_names> instead, or re-check every index.\n",
        name.c_str(), nmet, nrxn, wantMet, wantRxn);
    return b;
}

inline std::string verifyAgainstManifest(const std::string &path, const ManifestEntry &e)
{
    if (!e.found) {
        return "  [MODEL] not in models/manifest.txt, so its revision cannot be checked.\n"
               "  [MODEL] If you index exchanges positionally, confirm the indices yourself.\n";
    }
    std::string s;
    char b[512];

    unsigned long long h = 0;
    size_t bytes = 0;
    if (e.hash != 0 && hashFile(path, h, bytes)) {
        if (h != e.hash) {
            std::snprintf(b, sizeof(b),
                "  [MODEL] note: checksum differs from the manifest (%llu bytes here, %llu expected).\n"
                "  [MODEL]       Harmless if the file was merely reformatted; see the counts below.\n",
                (unsigned long long) bytes, (unsigned long long) e.bytes);
            s += b;
        }
    }

    if (e.rawSpecies > 0 && e.rawReactions > 0) {
        int ns = 0, nr = 0;
        countSbmlEntities(path, ns, nr);
        if (ns > 0 && nr > 0 && (ns != e.rawSpecies || nr != e.rawReactions)) {
            std::snprintf(b, sizeof(b),
                "  [MODEL] WARNING: this looks like a DIFFERENT REVISION of %s.\n"
                "  [MODEL]   manifest expects %d species and %d reaction tags; this file has %d and %d.\n"
                "  [MODEL]   <exchange_reaction_indices> is POSITIONAL, so your indices may now\n"
                "  [MODEL]   point at different reactions and the run would not complain.\n"
                "  [MODEL]   Use <exchange_reaction_names> instead, or re-check the indices.\n",
                e.name.c_str(), e.rawSpecies, e.rawReactions, ns, nr);
            s += b;
        } else if (ns > 0) {
            std::snprintf(b, sizeof(b),
                "  [MODEL] revision check: %d species, %d reaction tags, as the manifest expects\n", ns, nr);
            s += b;
        }
    }
    return s;
}

inline Result ensureModel(const std::string &source, const std::string &cacheDir,
                          bool allowDownload, bool isMaster,
                          const std::string &bundleDir = "models")
{
    Result R;
    std::string url, id, err;
    if (!resolveUrl(source, url, id, err)) { R.message = err; return R; }

    const std::string dir = cacheDir.empty() ? std::string(".") : cacheDir;
    R.path = dir + "/" + id + ".xml";
    const ManifestEntry entry = readManifest(bundleDir, id);
    if (entry.found && !entry.url.empty()) url = entry.url;

    /* ---- 1. already unpacked? ---- */
    if (fileExists(R.path)) {
        if (!looksLikeSbml(R.path)) {
            R.message = "  [MODEL] '" + R.path + "' exists but is not a model file. It is probably a\n"
                        "  [MODEL] failed download (servers answer a missing model with an HTML page\n"
                        "  [MODEL] surprisingly often). Delete it and try again.\n";
            return R;
        }
        R.ok = true;
        R.message = "  [MODEL] using " + R.path + "\n" + verifyAgainstManifest(R.path, entry);
        return R;
    }

    /* ---- 2. bundled with the code? ---- */
    const std::string gz = bundleDir + "/" + id + ".xml.gz";
    if (fileExists(gz)) {
        if (!isMaster) {
            R.message = "  [MODEL] waiting for rank 0 to unpack " + gz + "\n";
            R.ok = fileExists(R.path);
            return R;
        }
        std::string msg = "  [MODEL] unpacking bundled " + gz + "\n";
        if (gunzipTo(gz, R.path, msg)) {
            R.ok = true;
            R.fromBundle = true;
            msg += verifyAgainstManifest(R.path, entry);
        }
        R.message = msg;
        return R;
    }

    /* ---- 3. download ---- */
    if (!allowDownload) {
        R.message =
            "  [MODEL] '" + id + "' is not bundled, not cached, and <allow_download> is false.\n"
            "  [MODEL] Fetch it once on a machine with network access:\n"
            "  [MODEL]     mkdir -p " + dir + " && curl -L -o " + R.path + " " + url + "\n"
            "  [MODEL] then run again. Compute nodes usually cannot reach the internet, so doing\n"
            "  [MODEL] it this way round is the right thing on a cluster.\n";
        return R;
    }

    if (!isMaster) {
        R.message = "  [MODEL] waiting for rank 0 to provide " + R.path + "\n";
        R.ok = fileExists(R.path);
        return R;
    }

    const std::string tmp = R.path + ".part";
    std::string cmd = "curl -fsSL -o '" + tmp + "' '" + url + "' 2>/dev/null";
    int rc = std::system(cmd.c_str());
    if (rc != 0) {
        cmd = "wget -q -O '" + tmp + "' '" + url + "' 2>/dev/null";
        rc = std::system(cmd.c_str());
    }

    if (rc != 0 || !fileExists(tmp) || !looksLikeSbml(tmp)) {
        std::remove(tmp.c_str());
        R.message =
            "  [MODEL] could not download " + url + "\n"
            "  [MODEL] Most likely this machine has no outbound network, which is normal on a\n"
            "  [MODEL] compute node. Fetch it on the login node instead:\n"
            "  [MODEL]     mkdir -p " + dir + " && curl -L -o " + R.path + " " + url + "\n"
            "  [MODEL] and set <allow_download>false</allow_download> for the batch run.\n";
        return R;
    }

    if (std::rename(tmp.c_str(), R.path.c_str()) != 0) {
        std::remove(tmp.c_str());
        R.message = "  [MODEL] downloaded " + url + " but could not move it into place at "
                    + R.path + "\n";
        return R;
    }

    R.ok = true;
    R.downloaded = true;
    R.message = "  [MODEL] downloaded " + url + "\n  [MODEL] cached at " + R.path
              + " (delete it to force a refresh)\n"
              + verifyAgainstManifest(R.path, entry);
    return R;
}

}  // namespace complab_fetch

#endif  // COMPLAB3D_FETCH_HH
