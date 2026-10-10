// Stress test for EB generation with geometries prone to roundoff: surfaces
// through grid nodes, faces on grid planes, and nodes a few ulps off the
// surface.  Each case is generated from a seed and built from an implicit
// function or (3D only) a watertight STL file, then checked against the exact
// geometry.  Only geometries the EB supports are generated: unions of convex
// pieces a few cells thick and apart.  See the inputs file for the options.

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IndexSpace_STL.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_WriteEBSurface.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

namespace {

constexpr int max_pieces = 16;

// In 2D the pieces do not depend on z and are evaluated at z = 0.
constexpr Real huge_len = 1.e30_rt;

enum PieceKind : int { AxisBox = 0, RotBox, Sphere, Cylinder, Polyhedron };

// A convex piece.  The body is the union of the pieces.
//   AxisBox:    lo p[0..2], hi p[3..5]
//   RotBox:     center p[0..2], half sizes p[3..5], rows of rotation p[6..14]
//   Sphere:     center p[0..2], radius p[3]
//   Cylinder:   center p[0..2], radius p[3], axis direction p[4]
//   Polyhedron: planes [begin,end) of (n,c) with n.x <= c inside
struct Piece
{
    int kind = AxisBox;
    int begin = 0;
    int end = 0;
    Real p[15] = {};
};

struct Shape
{
    int n = 0;
    Piece piece[max_pieces];
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real box_sd (Real const q[3], Real const lo[3], Real const hi[3])
{
    Real dist_out2 = 0.0_rt;
    Real din = std::numeric_limits<Real>::max();
    for (int d = 0; d < 3; ++d) {
        Real const v = amrex::max(lo[d]-q[d], q[d]-hi[d]);
        if (v > 0.0_rt) { dist_out2 += v*v; }
        din = amrex::min(din, -v);
    }
    return (dist_out2 > 0.0_rt) ? -std::sqrt(dist_out2) : din;
}

// Signed distance to a piece, positive inside.  For a polyhedron the
// magnitude outside is a lower bound.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real piece_sd (Piece const& pc, Real const* planes, Real x, Real y, Real z)
{
    Real const q[3] = {x, y, z};
    Real const* p = pc.p;
    switch (pc.kind) {
    case AxisBox:
        return box_sd(q, p, p+3);
    case RotBox: {
        Real loc[3], lo[3], hi[3];
        for (int d = 0; d < 3; ++d) {
            loc[d] = p[6+3*d]*(x-p[0]) + p[7+3*d]*(y-p[1]) + p[8+3*d]*(z-p[2]);
            lo[d] = -p[3+d];
            hi[d] =  p[3+d];
        }
        return box_sd(loc, lo, hi);
    }
    case Sphere:
        return p[3] - std::sqrt((x-p[0])*(x-p[0]) + (y-p[1])*(y-p[1]) + (z-p[2])*(z-p[2]));
    case Cylinder: {
        int const dir = static_cast<int>(p[4]);
        Real r2 = 0.0_rt;
        for (int d = 0; d < 3; ++d) {
            if (d != dir) { r2 += (q[d]-p[d])*(q[d]-p[d]); }
        }
        return p[3] - std::sqrt(r2);
    }
    default: {
        Real f = std::numeric_limits<Real>::max();
        for (int i = pc.begin; i < pc.end; ++i) {
            Real const* pl = planes + 4*i;
            f = amrex::min(f, pl[3] - (pl[0]*x + pl[1]*y + pl[2]*z));
        }
        return f;
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real shape_sd (Shape const& s, Real const* planes, Real x, Real y, Real z)
{
    Real f = std::numeric_limits<Real>::lowest();
    for (int i = 0; i < s.n; ++i) {
        f = amrex::max(f, piece_sd(s.piece[i], planes, x, y, z));
    }
    return f;
}

struct ShapeIF : GPUable
{
    Shape shape;
    Real const* planes = nullptr;

    AMREX_GPU_HOST_DEVICE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
#if (AMREX_SPACEDIM == 3)
        return shape_sd(shape, planes, x, y, z);
#else
        return shape_sd(shape, planes, x, y, 0.0_rt);
#endif
    }

    Real operator() (RealArray const& p) const noexcept
    {
        return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
    }
};

// splitmix64
struct Rng
{
    std::uint64_t s;
    std::uint64_t next ()
    {
        std::uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    int uniform (int lo, int hi) { return lo + static_cast<int>(next() % std::uint64_t(hi-lo+1)); }
    double real () { return double(next() >> 11) * (1.0/9007199254740992.0); }
};

using V3 = std::array<Real,3>;

V3 operator* (Real s, V3 const& a) { return {s*a[0], s*a[1], s*a[2]}; }
Real dot (V3 const& a, V3 const& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
#if (AMREX_SPACEDIM == 3)
V3 operator+ (V3 const& a, V3 const& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
V3 operator- (V3 const& a, V3 const& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
V3 cross (V3 const& a, V3 const& b)
{
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
#endif
V3 normalize (V3 const& a) { return (1.0_rt/std::sqrt(dot(a,a))) * a; }

struct Triangle { V3 v[3]; };

struct Case
{
    std::string desc;
    bool stl = false;
    bool use_bvh = true;
    Box domain;
    RealBox rb;
    int max_grid_size = 32;
    int eb_max_grid_size = 64;
    int max_coarsening_level = 0;
    bool sharp = false; // has edges or corners
    Real vf_tol = 0.5_rt;
    Shape shape;
    std::vector<Real> planes;        // for Polyhedron pieces
    std::vector<Triangle> triangles; // STL
};

// Grid node coordinate, computed as the EB generation does; 0 for z in 2D
Real node (Geometry const& geom, int d, int i)
{
    return (d < AMREX_SPACEDIM) ? geom.ProbLo(d) + static_cast<Real>(i)*geom.CellSize(d)
                                : 0.0_rt;
}

// A coordinate on or near node i: exactly on it, a few ulps off, or at least
// 1/100 of a cell off.  Offsets in between only make cells that small cell
// removal discards by design.
Real near_node (Rng& rng, Geometry const& geom, int d, int i)
{
    if (d >= AMREX_SPACEDIM) { return 0.0_rt; }
    Real x = node(geom, d, i);
    Real const h = geom.CellSize(d);
    switch (rng.uniform(0, 8)) {
    case 0: case 1: case 2: return x;
    case 3: case 4: {
        int const k = rng.uniform(1, 3);
        Real const to = rng.uniform(0,1) ? std::numeric_limits<Real>::max()
                                         : std::numeric_limits<Real>::lowest();
        for (int n = 0; n < k; ++n) { x = std::nextafter(x, to); }
        return x;
    }
    case 5: return x + Real(0.01)*h;
    case 6: return x + Real(0.25)*h;
    case 7: return x + Real(0.5)*h;
    default: return x + Real(0.01 + 0.98*rng.real())*h;
    }
}

void add_box_triangles (std::vector<Triangle>& tris, V3 const& c, V3 const& h, Real const R[9])
{
    // corners in local coordinates, rotated back: x = c + R^T loc
    auto corner = [&] (int i, int j, int k) {
        Real const loc[3] = {(i ? h[0] : -h[0]), (j ? h[1] : -h[1]), (k ? h[2] : -h[2])};
        V3 x;
        for (int d = 0; d < 3; ++d) {
            x[d] = c[d] + R[d]*loc[0] + R[3+d]*loc[1] + R[6+d]*loc[2];
        }
        return x;
    };
    // faces as quads, counterclockwise seen from outside
    int const quads[6][4][3] = {{{0,0,0},{0,0,1},{0,1,1},{0,1,0}},
                                {{1,0,0},{1,1,0},{1,1,1},{1,0,1}},
                                {{0,0,0},{1,0,0},{1,0,1},{0,0,1}},
                                {{0,1,0},{0,1,1},{1,1,1},{1,1,0}},
                                {{0,0,0},{0,1,0},{1,1,0},{1,0,0}},
                                {{0,0,1},{1,0,1},{1,1,1},{0,1,1}}};
    for (auto const& q : quads) {
        V3 const a = corner(q[0][0],q[0][1],q[0][2]);
        V3 const b = corner(q[1][0],q[1][1],q[1][2]);
        V3 const cc = corner(q[2][0],q[2][1],q[2][2]);
        V3 const dd = corner(q[3][0],q[3][1],q[3][2]);
        tris.push_back({{a, b, cc}});
        tris.push_back({{a, cc, dd}});
    }
}

// Sphere approximated by a subdivided octahedron with the six axis vertices
// on the sphere; also returns the planes of its faces.  In 2D, a regular
// polygon with 4*2^level vertices including the four axis points.
void add_octasphere (std::vector<Triangle>& tris, std::vector<Real>& planes,
                     V3 const& c, Real r, int level)
{
#if (AMREX_SPACEDIM == 3)
    std::vector<std::array<V3,3>> faces;
    for (int sx = -1; sx <= 1; sx += 2) {
        for (int sy = -1; sy <= 1; sy += 2) {
            for (int sz = -1; sz <= 1; sz += 2) {
                V3 const a{Real(sx),0,0}, b{0,Real(sy),0}, cz{0,0,Real(sz)};
                if (sx*sy*sz > 0) { faces.push_back({a, b, cz}); }
                else              { faces.push_back({a, cz, b}); }
            }
        }
    }
    for (int l = 0; l < level; ++l) {
        std::vector<std::array<V3,3>> fine;
        for (auto const& f : faces) {
            V3 const m01 = normalize(f[0]+f[1]);
            V3 const m12 = normalize(f[1]+f[2]);
            V3 const m20 = normalize(f[2]+f[0]);
            fine.push_back({f[0], m01, m20});
            fine.push_back({m01, f[1], m12});
            fine.push_back({m20, m12, f[2]});
            fine.push_back({m01, m12, m20});
        }
        faces = std::move(fine);
    }
    for (auto const& f : faces) {
        Triangle t;
        for (int k = 0; k < 3; ++k) { t.v[k] = c + r*f[k]; }
        tris.push_back(t);
        V3 const n = normalize(cross(t.v[1]-t.v[0], t.v[2]-t.v[0]));
        planes.insert(planes.end(), {n[0], n[1], n[2], dot(n, t.v[0])});
    }
#else
    amrex::ignore_unused(tris);
    int const nv = 4 << level;
    for (int k = 0; k < nv; ++k) {
        double const t0 = 2.0*Math::pi<double>()*double(k)/double(nv);
        double const t1 = 2.0*Math::pi<double>()*double(k+1)/double(nv);
        // exact axis points
        auto vert = [&] (int kk, double t) -> V3 {
            int const q = kk % nv;
            if (q*4 == 0)    { return {c[0]+r, c[1], 0}; }
            if (q*4 == nv)   { return {c[0], c[1]+r, 0}; }
            if (q*4 == 2*nv) { return {c[0]-r, c[1], 0}; }
            if (q*4 == 3*nv) { return {c[0], c[1]-r, 0}; }
            return {c[0] + r*Real(std::cos(t)), c[1] + r*Real(std::sin(t)), 0};
        };
        V3 const a = vert(k, t0);
        V3 const b = vert(k+1, t1);
        V3 const n = normalize(V3{b[1]-a[1], a[0]-b[0], 0});
        planes.insert(planes.end(), {n[0], n[1], 0.0_rt, dot(n, a)});
    }
#endif
}

// Rotation about axis `axis` by the angle whose tangent is a/b
void rotation (int axis, int a, int b, Real R[9])
{
    Real const len = std::sqrt(Real(a*a + b*b));
    Real const cs = Real(b)/len;
    Real const sn = Real(a)/len;
    for (int i = 0; i < 9; ++i) { R[i] = 0.0_rt; }
    int const u = (axis+1)%3;
    int const v = (axis+2)%3;
    R[3*axis+axis] = 1.0_rt;
    R[3*u+u] =  cs; R[3*u+v] = sn;
    R[3*v+u] = -sn; R[3*v+v] = cs;
}

Case make_case (int icase, std::uint64_t seed)
{
    Rng rng{seed + 0x632be59bd9b4e019ULL * std::uint64_t(icase+1)};
    Case c;
    std::ostringstream desc;
    constexpr bool is3d = (AMREX_SPACEDIM == 3);

    // Grid: cell counts (log uniform), scale, origin and aspect ratio vary.
    double const nmin = is3d ? 20 : 32;
    double const nmax = is3d ? 64 : 512;
    auto ncell = [&] () {
        return int(std::lround(nmin * std::pow(nmax/nmin, rng.real())));
    };
    std::array<int,3> n{ncell(), ncell(), is3d ? ncell() : 1};
    // Feature sizes and counts grow with the domain
    int const nsmall = is3d ? std::min({n[0], n[1], n[2]}) : std::min(n[0], n[1]);
    int const grow = std::max(1, nsmall/32);
    c.domain = Box(IntVect(0), IntVect(AMREX_D_DECL(n[0]-1, n[1]-1, n[2]-1)));
    double const scales[] = {1.e-3, 0.37, 1.0, 4.2, 1.e3};
    double const s = scales[rng.uniform(0,4)];
    bool const aniso = rng.uniform(0,3) == 0;
    double const aspect[] = {1.0, 1.5, 0.75};
    RealArray lo, hi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        double const h = s / 32.0 * (aniso ? aspect[rng.uniform(0,2)] : 1.0);
        // up to 17000 cells from the origin, where roundoff in single
        // precision is still well below a cell
        double const origin[] = {0.0, -1.234*s, 0.3*s, -0.5*n[d]*h, -1234.5*h, 17000.25*h};
        lo[d] = Real(origin[rng.uniform(0,5)]);
        hi[d] = Real(double(lo[d]) + n[d]*h);
    }
    c.rb = RealBox(lo, hi);
    Geometry const geom(c.domain, c.rb, 0, {AMREX_D_DECL(0,0,0)});
    c.max_grid_size = 8 << rng.uniform(0,3);
    c.eb_max_grid_size = 16 << rng.uniform(0,2);
    c.max_coarsening_level = rng.uniform(0,1) ? 2 : 0;
    c.use_bvh = rng.uniform(0,1);
    desc << AMREX_SPACEDIM << "d n=" << n[0] << "x" << n[1];
    if (is3d) { desc << "x" << n[2]; }
    desc << " scale=" << s << (aniso ? " aniso" : "") << " mgs=" << c.max_grid_size
         << " ebmgs=" << c.eb_max_grid_size << " coarsen=" << c.max_coarsening_level;

    Real dxmin = std::numeric_limits<Real>::max();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { dxmin = std::min(dxmin, geom.CellSize(d)); }

    Shape& sh = c.shape;
    std::vector<RealBox> bounds; // of the pieces so far

    // At least 3 cells apart in some direction
    auto separated = [&] (RealBox const& rb) {
        for (auto const& o : bounds) {
            bool apart = false;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real const gap = amrex::max(rb.lo(d), o.lo(d)) - amrex::min(rb.hi(d), o.hi(d));
                apart = apart || gap >= 3.0_rt*geom.CellSize(d);
            }
            if (!apart) { return false; }
        }
        return true;
    };

    auto piece_bounds = [&] (Piece const& pc) {
        Real const* p = pc.p;
        RealArray blo, bhi;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real ext = 0.0_rt;
            switch (pc.kind) {
            case AxisBox: blo[d] = p[d]; bhi[d] = p[3+d]; continue;
            case RotBox: ext = std::abs(p[6+d])*p[3] + std::abs(p[9+d])*p[4]
                             + std::abs(p[12+d])*p[5]; break;
            case Cylinder: ext = (d == static_cast<int>(p[4])) ? huge_len : p[3]; break;
            default: ext = p[3];
            }
            blo[d] = p[d] - ext;
            bhi[d] = p[d] + ext;
        }
        return RealBox(blo, bhi);
    };

    // A center on a node, a cell center, or near a node
    auto center = [&] (std::array<int,3> const& iv) {
        V3 cc{0,0,0};
        int const k = rng.uniform(0, 2);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            cc[d] = (k == 0) ? node(geom, d, iv[d])
                : ((k == 1) ? node(geom, d, iv[d]) + 0.5_rt*geom.CellSize(d)
                            : near_node(rng, geom, d, iv[d]));
        }
        return cc;
    };

    // Convex pieces of the given kinds placed at random, at least 3 cells apart
    auto add_pieces = [&] (std::vector<int> const& kinds, int npieces) {
        for (int b = 0; b < npieces; ++b) {
            for (int tries = 0; tries < 50; ++tries) {
                Piece pc;
                pc.kind = kinds[rng.uniform(0, int(kinds.size())-1)];
                std::array<int,3> iv{rng.uniform(-2, n[0]+1), rng.uniform(-2, n[1]+1),
                                     rng.uniform(-2, n[2]+1)};
                V3 const cc = center(iv);
                std::vector<Real> planes;
                std::vector<Triangle> tris;
                if (pc.kind == RotBox) {
                    int const ta = rng.uniform(1, 3);
                    int const tb = rng.uniform(0, 3);
                    Real const len = std::sqrt(Real(ta*ta + tb*tb));
                    for (int d = 0; d < 3; ++d) {
                        // a multiple of dx/len keeps the faces on nodes for square cells
                        pc.p[3+d] = Real(rng.uniform(int(3*len)+1, int(Real(6*grow)*len)+1)) * dxmin / len;
                    }
                    rotation(is3d ? rng.uniform(0, 2) : 2, ta, tb, pc.p+6);
                    if (!is3d) { pc.p[5] = huge_len; }
                    add_box_triangles(tris, cc, V3{pc.p[3],pc.p[4],pc.p[5]}, pc.p+6);
                } else if (pc.kind == Sphere || pc.kind == Cylinder) {
                    // radius through nodes for a center on a node
                    int const triples[][3] = {{4,0,0},{3,4,0},{2,3,6},{4,4,2},{5,0,0},{3,3,3}};
                    auto const& t = triples[rng.uniform(0,5)];
                    pc.p[3] = std::sqrt(Real(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]))
                        * Real(rng.uniform(1, grow)) * dxmin;
                    pc.p[4] = Real(is3d ? rng.uniform(0, 2) : 2);
                    if (!is3d) { pc.kind = Cylinder; }
                } else if (pc.kind == Polyhedron) {
                    pc.p[3] = Real(rng.uniform(4, 7*grow)) * dxmin;
                    add_octasphere(tris, planes, cc, pc.p[3], rng.uniform(1, 2));
                }
                for (int d = 0; d < 3; ++d) { pc.p[d] = cc[d]; }
                RealBox const rb = piece_bounds(pc);
                if (!separated(rb)) { continue; }
                bounds.push_back(rb);
                if (pc.kind == Polyhedron) {
                    pc.begin = int(c.planes.size()/4);
                    c.planes.insert(c.planes.end(), planes.begin(), planes.end());
                    pc.end = int(c.planes.size()/4);
                }
                if (c.stl) { c.triangles.insert(c.triangles.end(), tris.begin(), tris.end()); }
                sh.piece[sh.n++] = pc;
                break;
            }
        }
    };

    int const family = rng.uniform(0, 6);
    if (family == 0 || family == 1) {
        // Axis-aligned boxes on a lattice of node-related coordinates
        c.stl = is3d && (family == 1);
        c.sharp = true;
        desc << (c.stl ? " stl" : " if") << " grid-boxes";
        std::array<std::vector<Real>,3> lat;
        for (int d = 0; d < 3; ++d) {
            if (!is3d && d == 2) {
                lat[d] = {-huge_len, huge_len};
                continue;
            }
            int i = rng.uniform(-4, 4);
            int const m = rng.uniform(1, is3d ? 3 : std::min(4, 1+grow));
            for (int k = 0; k < m; ++k) {
                int const len = rng.uniform(4, 8*grow);
                int const gap = rng.uniform(4, 6*grow);
                lat[d].push_back(near_node(rng, geom, d, i));
                lat[d].push_back(near_node(rng, geom, d, i+len));
                i += len + gap;
            }
        }
        for (std::size_t a = 0; a < lat[0].size(); a += 2) {
            for (std::size_t b = 0; b < lat[1].size(); b += 2) {
                for (std::size_t e = 0; e < lat[2].size(); e += 2) {
                    if (sh.n == max_pieces) { continue; }
                    Piece& pc = sh.piece[sh.n++];
                    pc.kind = AxisBox;
                    pc.p[0] = lat[0][a]; pc.p[1] = lat[1][b]; pc.p[2] = lat[2][e];
                    pc.p[3] = lat[0][a+1]; pc.p[4] = lat[1][b+1]; pc.p[5] = lat[2][e+1];
                    if (c.stl) {
                        // the corners exactly, not center +- half size
                        V3 const cc{0.5_rt*(pc.p[0]+pc.p[3]), 0.5_rt*(pc.p[1]+pc.p[4]),
                                    0.5_rt*(pc.p[2]+pc.p[5])};
                        Real const Rid[9] = {1,0,0, 0,1,0, 0,0,1};
                        std::vector<Triangle> t;
                        add_box_triangles(t, cc, V3{1,1,1}, Rid);
                        for (auto& tr : t) {
                            for (auto& v : tr.v) {
                                for (int d = 0; d < 3; ++d) {
                                    v[d] = (v[d] < cc[d]) ? pc.p[d] : pc.p[3+d];
                                }
                            }
                            c.triangles.push_back(tr);
                        }
                    }
                }
            }
        }
    } else if (family == 2) {
        c.stl = is3d && rng.uniform(0, 1);
        c.sharp = true;
        desc << (c.stl ? " stl" : " if") << " rotated-boxes";
        add_pieces({RotBox}, rng.uniform(1, std::min(max_pieces, 4*grow)));
    } else if (family == 3) {
        desc << " if spheres-cylinders";
        add_pieces({Sphere, Cylinder}, rng.uniform(1, std::min(max_pieces, 4*grow)));
    } else if (family == 4) {
        c.stl = is3d && rng.uniform(0, 1);
        c.sharp = true;
        desc << (c.stl ? " stl" : " if") << " octaspheres";
        add_pieces({Polyhedron}, rng.uniform(1, std::min(max_pieces, 4*grow)));
    } else if (family == 5) {
        c.sharp = true;
        desc << " if mixed";
        add_pieces({RotBox, Sphere, Cylinder, Polyhedron},
                   rng.uniform(2, std::min(max_pieces, 6*grow)));
    } else {
        // A convex polytope of half spaces through nodes with rational
        // normals.  Normals at most 90 degrees apart keep the body edges
        // from being sharper than 90 degrees.
        c.sharp = true;
        desc << " if polytope";
        int const np = rng.uniform(1, 3);
        std::array<int,3> const iv{rng.uniform(8, n[0]-8), rng.uniform(8, n[1]-8),
                                   is3d ? rng.uniform(8, n[2]-8) : 0};
        std::vector<V3> normals;
        Piece& pc = sh.piece[sh.n++];
        pc.kind = Polyhedron;
        pc.begin = int(c.planes.size()/4);
        for (int b = 0; b < np; ++b) {
            V3 nv{};
            bool ok = false;
            for (int tries = 0; tries < 100 && !ok; ++tries) {
                nv = V3{Real(rng.uniform(-2,2)), Real(rng.uniform(-2,2)),
                        is3d ? Real(rng.uniform(-2,2)) : 0.0_rt};
                if (dot(nv,nv) == 0) { continue; }
                nv = normalize(nv);
                ok = true;
                for (auto const& m : normals) { ok = ok && (dot(nv, m) >= 0); }
            }
            if (!ok) { continue; }
            normals.push_back(nv);
            V3 pt;
            for (int d = 0; d < 3; ++d) { pt[d] = node(geom, d, iv[d] + rng.uniform(-2,2)); }
            c.planes.insert(c.planes.end(), {nv[0], nv[1], nv[2], dot(nv, pt)});
        }
        pc.end = int(c.planes.size()/4);
    }
    desc << " pieces=" << sh.n;
    if (c.stl) { desc << " triangles=" << c.triangles.size() << " bvh=" << c.use_bvh; }
    c.desc = desc.str();
    return c;
}

void write_stl (std::string const& fname, std::vector<Triangle> const& tris)
{
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream ofs(fname);
        ofs << std::setprecision(std::numeric_limits<Real>::max_digits10);
        ofs << "solid eb_stress\n";
        for (auto const& t : tris) {
            ofs << "facet normal 0 0 0\n outer loop\n";
            for (auto const& v : t.v) {
                ofs << "  vertex " << v[0] << " " << v[1] << " " << v[2] << "\n";
            }
            ofs << " endloop\nendfacet\n";
        }
        ofs << "endsolid eb_stress\n";
    }
    ParallelDescriptor::Barrier();
}

struct Result
{
    Long ncut = 0;
    Long nbad_nodes = 0;
    Long nbad_data = 0;
    Real max_dist = 0;
    Real max_vf_err = 0;
};

Result check (Case const& c, Geometry const& geom, BoxArray const& ba,
              DistributionMapping const& dm, Real const* planes, int verbose,
              int write_surface)
{
    auto factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
    auto const& flags = factory->getMultiEBCellFlagFab();
    auto const& vfrac = factory->getVolFrac();
    auto const& levset = factory->getLevelSet();

    GpuArray<Real,3> dx{0,0,0};
    GpuArray<Real,3> problo{0,0,0};
    Real dxmin = std::numeric_limits<Real>::max();
    Real dxmax = 0.0_rt;
    Real diag2 = 0.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = geom.CellSize(d);
        problo[d] = geom.ProbLo(d);
        dxmin = amrex::min(dxmin, dx[d]);
        dxmax = amrex::max(dxmax, dx[d]);
        diag2 += dx[d]*dx[d];
    }
    Real const diag = std::sqrt(diag2);
    Shape const shape = c.shape;
    Real const vf_tol = c.vf_tol;
    bool const sharp = c.sharp;

    auto const fa = flags.const_arrays();
    auto const vf_a = vfrac.const_arrays();
    auto const bc_a = factory->getBndryCent().const_arrays();
    auto const ba_a = factory->getBndryArea().const_arrays();
    auto const bn_a = factory->getBndryNormal().const_arrays();
    auto const ls_a = levset.const_arrays();
    Array<MultiFab,AMREX_SPACEDIM> ap;
    Array<MultiArray4<Real const>,AMREX_SPACEDIM> ap_a;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        ap[d] = factory->getAreaFrac()[d]->ToMultiFab(1.0_rt, 0.0_rt);
        ap_a[d] = ap[d].const_arrays();
    }

    Result res;

    // Per cell: volume fraction against the exact one (sampled near the
    // surface), EB centroid distance to the surface, and data sanity.
    auto r = ParReduce(TypeList<ReduceOpMax,ReduceOpMax,ReduceOpSum,ReduceOpSum>{},
                       TypeList<Real,Real,Long,Long>{}, flags, IntVect(0),
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> GpuTuple<Real,Real,Long,Long>
    {
        int const ijk[3] = {i, j, k};
        Real clo[3] = {0, 0, 0};
        Real cc[3] = {0, 0, 0};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            clo[d] = problo[d] + Real(ijk[d])*dx[d];
            cc[d] = clo[d] + 0.5_rt*dx[d];
        }
        Real const sdc = shape_sd(shape, planes, cc[0], cc[1], cc[2]);
        Real exact;
        bool const far = std::abs(sdc) > 0.51_rt*diag;
        if (far) {
            exact = (sdc > 0.0_rt) ? 0.0_rt : 1.0_rt;
        } else {
            constexpr int ns = 8;
            constexpr int nsz = (AMREX_SPACEDIM == 3) ? ns : 1;
            int nfluid = 0;
            for (int kk = 0; kk < nsz; ++kk) {
            for (int jj = 0; jj < ns; ++jj) {
            for (int ii = 0; ii < ns; ++ii) {
                Real const x = clo[0] + (Real(ii)+0.5_rt)/Real(ns)*dx[0];
                Real const y = clo[1] + (Real(jj)+0.5_rt)/Real(ns)*dx[1];
                Real const z = clo[2] + (Real(kk)+0.5_rt)/Real(ns)*dx[2];
                nfluid += (shape_sd(shape, planes, x, y, z) > 0.0_rt) ? 0 : 1;
            }}}
            exact = Real(nfluid) / Real(ns*ns*nsz);
        }
        Real const v = vf_a[b](i,j,k);
        // With sharp edges, a cell that looks cut twice is covered and its
        // nodes are set to zero.
        bool covered_cut = false;
        if (sharp && v == 0.0_rt && !far) {
            for (int kk = 0; kk < AMREX_SPACEDIM-1; ++kk) {
            for (int jj = 0; jj < 2; ++jj) {
            for (int ii = 0; ii < 2; ++ii) {
                covered_cut = covered_cut || ls_a[b](i+ii,j+jj,k+kk) == 0.0_rt;
            }}}
        }
        Real const vf_err = covered_cut ? 0.0_rt : std::abs(v - exact);

        auto const flag = fa[b](i,j,k);
        // fractions may be off by roundoff
        constexpr Real ftol = 16.0_rt*std::numeric_limits<Real>::epsilon();
        // away from the surface the cell must be regular or covered
        bool bad = (far && v != exact) || !amrex::isfinite(v) || v < -ftol || v > 1.0_rt+ftol;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            IntVect ivp(AMREX_D_DECL(i,j,k));
            ivp[d] += 1;
            Real const am = ap_a[d][b](i,j,k);
            Real const app = ap_a[d][b](ivp);
            bad = bad || !amrex::isfinite(am) || am < -ftol || am > 1.0_rt+ftol
                      || !amrex::isfinite(app) || app < -ftol || app > 1.0_rt+ftol;
        }
        if (flag.isRegular()) {
            bad = bad || v != 1.0_rt;
        } else if (flag.isCovered()) {
            bad = bad || v != 0.0_rt;
        }

        Real dist = 0.0_rt;
        Long cut = 0;
        if (flag.isSingleValued()) {
            cut = 1;
            Real q[3] = {0, 0, 0};
            Real bn2 = 0.0_rt;
            bool finite = true;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real const bc = bc_a[b](i,j,k,d);
                bn2 += bn_a[b](i,j,k,d)*bn_a[b](i,j,k,d);
                finite = finite && amrex::isfinite(bc);
                bad = bad || !amrex::isfinite(bc) || std::abs(bc) > 0.5_rt + 1.e-4_rt;
                q[d] = clo[d] + (0.5_rt+bc)*dx[d];
            }
            Real const area = ba_a[b](i,j,k);
            bad = bad || !amrex::isfinite(area) || area <= 0.0_rt
                      || std::abs(bn2 - 1.0_rt) > 1.e-3_rt;
            dist = finite ? std::abs(shape_sd(shape, planes, q[0], q[1], q[2])) / dxmax
                          : std::numeric_limits<Real>::max();
        }
        if (verbose > 1 && (bad || vf_err > vf_tol || dist > 2.0_rt)) {
            AMREX_DEVICE_PRINTF("  bad cell (%d,%d,%d) volfrac %g exact %g dist %g data %d\n",
                                i, j, k, double(v), double(exact), double(dist), int(bad));
        }
        return {vf_err, dist, Long(bad), cut};
    });
    res.max_vf_err = amrex::get<0>(r);
    res.max_dist = amrex::get<1>(r);
    res.nbad_data = amrex::get<2>(r);
    res.ncut = amrex::get<3>(r);

    // Nodes clearly off the surface must be on the correct side.  Clearly
    // means more than 0.1 dx and well beyond roundoff in the coordinates.
    // For shapes with sharp edges, a cell that looks cut twice is covered,
    // which moves the surface onto its nodes (zero level set).
    Real xmax = 0.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        xmax = amrex::max(xmax, std::abs(geom.ProbLo(d)), std::abs(geom.ProbHi(d)));
    }
    Real const off = amrex::max(0.1_rt*dxmin, 64.0_rt*std::numeric_limits<Real>::epsilon()*xmax);
    Real const zero_ok = sharp ? amrex::max(diag, off) : off;
    res.nbad_nodes = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Long>{}, levset, IntVect(0),
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> GpuTuple<Long>
    {
        int const ijk[3] = {i, j, k};
        Real q[3] = {0, 0, 0};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { q[d] = problo[d] + Real(ijk[d])*dx[d]; }
        Real const sd = shape_sd(shape, planes, q[0], q[1], q[2]);
        Real const ls = ls_a[b](i,j,k);
        bool const bad = std::abs(sd) > off && (sd > 0.0_rt) != (ls >= 0.0_rt)
            && !(ls == 0.0_rt && std::abs(sd) <= zero_ok);
        if (bad && verbose > 1) {
            AMREX_DEVICE_PRINTF("  bad node (%d,%d,%d) levelset %g distance %g\n",
                                i, j, k, double(ls), double(sd));
        }
        return {Long(bad)};
    });

    ParallelDescriptor::ReduceRealMax({res.max_vf_err, res.max_dist});
    ParallelDescriptor::ReduceLongSum({res.nbad_nodes, res.nbad_data, res.ncut});

#if (AMREX_SPACEDIM == 3)
    if (write_surface) {
        WriteEBSurface(ba, dm, geom, factory.get());
    }
#else
    amrex::ignore_unused(ba, dm, write_surface);
#endif
    return res;
}

// Everything needed to rerun a case, for the failure log
std::string case_record (Case const& c, Geometry const& geom, Long seed, int icase,
                         std::string const& cmdline)
{
    std::ostringstream os;
    os << std::setprecision(std::numeric_limits<Real>::max_digits10);
    os << "rerun: " << cmdline << " seed=" << seed << " case=" << icase
       << "  (same AMReX version)\n";
    os << "  amrex " << amrex::Version() << ", " << AMREX_SPACEDIM << "D, Real is "
       << (sizeof(Real) == sizeof(double) ? "double" : "float")
       << ", MPI ranks " << ParallelDescriptor::NProcs()
#if defined(AMREX_USE_CUDA)
       << ", CUDA"
#elif defined(AMREX_USE_HIP)
       << ", HIP"
#elif defined(AMREX_USE_SYCL)
       << ", SYCL"
#endif
       << "\n";
    os << "  seed " << seed << " case " << icase << ": " << c.desc << "\n";
    os << "  domain " << c.domain << " prob_lo";
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { os << " " << geom.ProbLo(d); }
    os << " dx";
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { os << " " << geom.CellSize(d); }
    os << "\n";
    for (int ip = 0; ip < c.shape.n; ++ip) {
        auto const& pc = c.shape.piece[ip];
        os << "  piece " << ip << " kind " << pc.kind << " planes [" << pc.begin << ","
           << pc.end << ") p";
        for (Real v : pc.p) { os << " " << v; }
        os << "\n";
    }
    for (std::size_t ip = 0; ip+3 < c.planes.size(); ip += 4) {
        os << "  plane " << ip/4 << " n " << c.planes[ip] << " " << c.planes[ip+1] << " "
           << c.planes[ip+2] << " c " << c.planes[ip+3] << "\n";
    }
    return os.str();
}

std::string failure_log;  // appended to, and closed, at each failure
std::string current_case; // record of the case being run

// Called by amrex::Abort: log the case that was running, then abort.
void abort_handler (const char* msg)
{
    amrex::SetErrorHandler(nullptr);
    if (!current_case.empty()) {
        std::ofstream ofs(failure_log, std::ios::app);
        ofs << "ABORTED (rank " << ParallelDescriptor::MyProc() << "): "
            << (msg ? msg : "") << "\n" << current_case << "\n";
    }
    amrex::Abort(msg ? msg : "EB stress test aborted");
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        int ncases = 100;
        pp.query("ncases", ncases);
        int first_case = 0;
        pp.query("first_case", first_case);
        int single_case = -1;
        pp.query("case", single_case);
        // random unless given; printed so that a failure can be rerun
        Long seed = -1;
        pp.query("seed", seed);
        if (seed < 0) {
            if (ParallelDescriptor::IOProcessor()) {
                std::random_device rd;
                seed = Long((std::uint64_t(rd()) << 31) ^ std::uint64_t(rd())) & 0x7fffffffffffLL;
            }
            ParallelDescriptor::Bcast(&seed, 1, ParallelDescriptor::IOProcessorNumber());
        }
        // 0: only failures and the summary, 1: a line per case, 2: bad cells
        int verbose = 1;
        pp.query("verbose", verbose);
        failure_log = "eb_stress_failures.log";
        pp.query("failure_log", failure_log);
        // Largest volume fraction error of a cut cell.  With one cut per cell
        // based on node values, an edge or corner just inside a cell can cost
        // most of it; smooth surfaces cost far less.
        Real max_vf_err = 0.4_rt;
        pp.query("max_volfrac_error", max_vf_err);
        Real max_vf_err_sharp = 0.9_rt;
        pp.query("max_volfrac_error_sharp", max_vf_err_sharp);
        int write_surface = 0;
        pp.query("write_surface", write_surface);
        {
            // Sharp edges crossing faces obliquely can make cells that look
            // cut twice; cover them unless told otherwise.
            ParmParse ppeb("eb2");
            bool cover_multiple_cuts = true;
            ppeb.queryAdd("cover_multiple_cuts", cover_multiple_cuts);
        }
        // If positive, run cases until this many seconds have passed
        // instead of ncases cases.
        double run_time = 0;
        pp.query("run_time", run_time);
        // If positive, print progress every this many seconds
        double progress_interval = 0;
        pp.query("progress_interval", progress_interval);
        if (single_case >= 0) { run_time = 0; }

        int const cbegin = (single_case >= 0) ? single_case : first_case;
        int const cend = (single_case >= 0) ? single_case+1 : first_case+ncases;
        if (run_time > 0) {
            amrex::Print() << "EB stress test seed = " << seed << ", cases from " << cbegin
                           << " for " << run_time << " seconds\n";
        } else {
            amrex::Print() << "EB stress test seed = " << seed << ", cases " << cbegin << " to "
                           << cend-1 << "\n";
        }
        double const start_time = amrex::second();
        double next_progress = progress_interval;

        std::string cmdline = argv[0];
        for (int i = 1; i < argc; ++i) { cmdline += std::string(" ") + argv[i]; }

        amrex::SetErrorHandler(abort_handler);

        std::vector<int> failed;
        int icase = cbegin;
        for (; ; ++icase) {
            if (run_time > 0) {
                int stop = (amrex::second() - start_time >= double(run_time)) ? 1 : 0;
                ParallelDescriptor::Bcast(&stop, 1, ParallelDescriptor::IOProcessorNumber());
                if (stop) { break; }
            } else if (icase >= cend) {
                break;
            }
            Case c = make_case(icase, std::uint64_t(seed));
            c.vf_tol = c.sharp ? max_vf_err_sharp : max_vf_err;
            Geometry const geom(c.domain, c.rb, 0, {AMREX_D_DECL(0,0,0)});
            BoxArray ba(c.domain);
            ba.maxSize(c.max_grid_size);
            DistributionMapping const dm(ba);
            current_case = case_record(c, geom, seed, icase, cmdline);
            if (verbose > 2) { amrex::Print() << current_case; }

            Gpu::DeviceVector<Real> planes(c.planes.size());
            Gpu::copyAsync(Gpu::hostToDevice, c.planes.begin(), c.planes.end(), planes.begin());
            Gpu::streamSynchronize();
            Real const* pplanes = planes.data();

            int const saved_max_grid_size = EB2::max_grid_size;
            EB2::max_grid_size = c.eb_max_grid_size;
            std::string stl_file;
            if (c.stl) {
                stl_file = "eb_stress_" + std::to_string(seed) + "_" + std::to_string(icase) + ".stl";
                write_stl(stl_file, c.triangles);
                EB2::IndexSpace::push(std::make_unique<EB2::IndexSpaceSTL>
                    (stl_file, 1.0_rt, Array<Real,3>{0.0_rt,0.0_rt,0.0_rt}, 0, geom, 0,
                     c.max_coarsening_level, 4, true, true, 0, c.use_bvh, false));
            } else {
                ShapeIF f;
                f.shape = c.shape;
                f.planes = pplanes;
                EB2::Build(EB2::makeShop(f), geom, 0, c.max_coarsening_level);
            }
            EB2::max_grid_size = saved_max_grid_size;

            Result const r = check(c, geom, ba, dm, pplanes, verbose, write_surface);
            EB2::IndexSpace::pop();

            bool const ok = r.nbad_nodes == 0 && r.nbad_data == 0 && r.max_dist <= 2.0_rt
                && r.max_vf_err <= c.vf_tol;
            std::ostringstream line;
            line << "case " << icase << " " << (ok ? "ok  " : "FAIL") << " [" << c.desc
                 << "] cut=" << r.ncut << " vf_err=" << std::setprecision(3) << r.max_vf_err
                 << " dist=" << r.max_dist << " bad_nodes=" << r.nbad_nodes
                 << " bad_data=" << r.nbad_data;
            if (verbose > 0 || !ok) { amrex::Print() << line.str() << "\n"; }

            // Log a failure right away and keep its STL file; nothing is
            // left behind for a case that passes.
            ParallelDescriptor::Barrier();
            if (ParallelDescriptor::IOProcessor()) {
                if (!ok) {
                    std::ofstream ofs(failure_log, std::ios::app);
                    ofs << "FAILED: " << line.str() << "\n" << current_case;
                    if (c.stl) { ofs << "  stl file " << stl_file << "\n"; }
                    ofs << "\n";
                } else if (c.stl) {
                    std::remove(stl_file.c_str());
                }
            }
            if (!ok) { failed.push_back(icase); }
            current_case.clear();

            if (progress_interval > 0) {
                double const elapsed = amrex::second() - start_time;
                if (elapsed >= next_progress) {
                    int const ndone = icase - cbegin + 1;
                    amrex::Print() << "progress: " << ndone << " cases (" << cbegin << " to "
                                   << icase << "), " << failed.size() << " failed, "
                                   << std::fixed << std::setprecision(0) << elapsed << " s, "
                                   << std::setprecision(1) << double(ndone)/elapsed
                                   << " cases/s\n" << std::defaultfloat;
                    while (next_progress <= elapsed) { next_progress += progress_interval; }
                }
            }
        }

        amrex::SetErrorHandler(nullptr);
        int const nrun = icase - cbegin;
        amrex::Print() << "EB stress test: " << (nrun-int(failed.size())) << " of " << nrun
                       << " cases (" << cbegin << " to " << icase-1 << ") passed with seed = "
                       << seed << "\n";
        if (!failed.empty()) {
            amrex::Print() << "Failed cases (rerun with seed=" << seed << " case=N):";
            for (int f : failed) { amrex::Print() << " " << f; }
            amrex::Print() << "\nDetails in " << failure_log << "\n";
            amrex::Abort("EB stress test failed");
        }
    }
    amrex::Finalize();
}
