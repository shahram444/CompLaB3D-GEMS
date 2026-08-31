/* Just enough of Palabos to compile the two new processors on their own.
 * Every name and signature here was copied from what the real headers use. */
#ifndef PALABOS_STUB_HH
#define PALABOS_STUB_HH
#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

typedef long long plint;
typedef double    T;      /* Palabos convention: complab.cpp typedefs T before any include */

namespace plb {

typedef ::plint plint;                 /* the real Palabos plint lives in namespace plb */
const double thrd = 1e-12;             /* complab3d_processors.hh's global threshold */

struct Box3D { plint x0,x1,y0,y1,z0,z1; };
struct Dot3D { plint x,y,z; };

template<typename T, int N> struct Array { T v[N]; T& operator[](int i){return v[i];} const T& operator[](int i) const {return v[i];} };

namespace util { inline plint roundToInt(double d){ return (plint) std::floor(d+0.5); } }
namespace global { struct Timer { void restart(){} double getTime(){return 0;} void stop(){} };
                   inline Timer& timer(const std::string&){ static Timer t; return t; } }
namespace modif { enum ModifT { nothing, staticVariables }; }
namespace BlockDomain { enum DomainT { bulk, bulkAndEnvelope }; }

template<typename T> struct D3Q7 { };

/* A cell that actually holds something, so a processor can be RUN and not merely compiled.
 * rho is what computeDensity() reports; p is the D3Q7 population the processor deposits into. */
template<typename T, template<typename U> class Descriptor>
struct Cell {
    T rho;
    Array<T,7> p;
    Cell() : rho(T()) { for (int i = 0; i < 7; ++i) p[i] = T(); }
    void getPopulations(Array<T,7>& g) const { g = p; }
    void setPopulations(const Array<T,7>& g) { p = g; }
    /* Real Palabos reports the zeroth moment, so anything written into the populations shows up
     * in the density.  The stub has to do the same or a processor that deposits into a lattice
     * and then reads it back looks like a no-op. */
    T computeDensity() const { T s = rho; for (int i = 0; i < 7; ++i) s += p[i]; return s; }
};

template<typename T, template<typename U> class Descriptor>
struct BlockLattice3D {
    Cell<T,Descriptor> c;
    Cell<T,Descriptor>& get(plint,plint,plint) { return c; }
    Dot3D getLocation() const { Dot3D d = {0,0,0}; return d; }
    T deposited() const { T s = T(); for (int i = 0; i < 7; ++i) s += c.p[i]; return s; }
};

template<typename T, template<typename U> class Descriptor>
inline Dot3D computeRelativeDisplacement(const BlockLattice3D<T,Descriptor>&,
                                         const BlockLattice3D<T,Descriptor>&)
{ Dot3D d = {0,0,0}; return d; }

template<typename T2> struct ScalarField3D {
    T2 v;
    T2& get(plint,plint,plint){ return v; }
};

/* the scalar-field-plus-lattice functional base initPhaseLattice3D derives from */
template<typename T1, template<typename U> class Descriptor, typename T2>
struct BoxProcessingFunctional3D_LS {
    virtual ~BoxProcessingFunctional3D_LS() {}
};

template<typename T, template<typename U> class Descriptor>
struct LatticeBoxProcessingFunctional3D {
    virtual ~LatticeBoxProcessingFunctional3D() {}
    virtual void process(Box3D, std::vector<BlockLattice3D<T,Descriptor>*>) = 0;
};

/* the XML reader and the log stream, only as far as the metabolic header uses them */
struct PlbIOException { const char* what() const { return "stub"; } };
struct XMLreader {
    XMLreader() {}
    explicit XMLreader(const std::string&) {}
    XMLreader& operator[](const std::string&) { return *this; }
    XMLreader& operator[](plint) { return *this; }
    template<typename V> void read(V&) const { throw PlbIOException(); }
};
struct PlbOut {
    template<typename V> PlbOut& operator<<(const V&) { return *this; }
    PlbOut& operator<<(std::ostream& (*)(std::ostream&)) { return *this; }
};
extern PlbOut pcout;
inline std::string createFileName(std::string s, plint, plint) { return s; }

}  // namespace plb
#endif
