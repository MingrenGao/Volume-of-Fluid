// curvature_convergence_paraboloid_fixed.cpp

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <array>
#include <stdexcept>
#include <limits>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ====================== Helpers ======================
static inline double clamp01(double x) { return std::min(1.0, std::max(0.0, x)); }
static std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

// ====================== Simple 2D grid wrapper ======================
struct Grid2D {
    int nx, ny;
    double dx, dy;
    std::vector<double> f;   // size = nx * ny

    Grid2D() : nx(0), ny(0), dx(0.0), dy(0.0) {}
    Grid2D(int nx_, int ny_, double dx_, double dy_)
        : nx(nx_), ny(ny_), dx(dx_), dy(dy_), f(nx_ * ny_, 0.0) {}

    inline int idx(int i, int j) const { return j * nx + i; }
    double& operator()(int i, int j)       { return f[idx(i, j)]; }
    double  operator()(int i, int j) const { return f[idx(i, j)]; }
};

// ====================== Radial kernel K(r) =========================
// K(r) = (1 - r^2)^4, r < 1; otherwise 0
inline double kernel(double r) {
    if (r >= 1.0) return 0.0;
    double t  = 1.0 - r * r;
    double t2 = t * t;
    return t2 * t2;
}

// ====================== Generic convolution: g_tilde = K * g =======
Grid2D convolveField(const Grid2D& g, double h) {
    Grid2D gt(g.nx, g.ny, g.dx, g.dy);

    int rad_x = static_cast<int>(std::ceil(h / g.dx));
    int rad_y = static_cast<int>(std::ceil(h / g.dy));

    for (int j = 0; j < g.ny; ++j) {
        double yj = (j + 0.5) * g.dy;
        for (int i = 0; i < g.nx; ++i) {
            double xi = (i + 0.5) * g.dx;

            double sum_wf = 0.0;
            double sum_w  = 0.0;

            for (int jj = std::max(0, j - rad_y); jj <= std::min(g.ny - 1, j + rad_y); ++jj) {
                double yb = (jj + 0.5) * g.dy;
                for (int ii = std::max(0, i - rad_x); ii <= std::min(g.nx - 1, i + rad_x); ++ii) {
                    double xb = (ii + 0.5) * g.dx;

                    double dx = xi - xb;
                    double dy = yj - yb;
                    double r  = std::sqrt(dx * dx + dy * dy);

                    double w = kernel(r / h);
                    if (w == 0.0) continue;

                    sum_wf += w * g(ii, jj);
                    sum_w  += w;
                }
            }

            if (sum_w > 0.0) gt(i, j) = sum_wf / sum_w;
            else             gt(i, j) = g(i, j);
        }
    }

    return gt;
}
inline Grid2D convolveVOF(const Grid2D& g, double h) { return convolveField(g, h); }

// ====================== Config struct & parser =====================
struct Config {
    // domain
    double prob_lo_x = 0.0, prob_lo_y = 0.0;
    double prob_hi_x = 1.0, prob_hi_y = 1.0;

    int nx = 128, ny = 128;
    std::vector<int> n_list;

    // initial geometry: circle / cosine
    std::string ic_shape = "circle";

    // circle parameters
    double R  = 0.25, cx = 0.5, cy = 0.5;

    // cosine parameters
    double cos_A = 3.0, cos_L = 8.0;
    int    cos_n = 1;

    // smoothing radius: h = h_factor * dx
    double h_factor = 4.0;

    // local recursive refinement depth for VOF initialization
    int vof_refine_max_level = 3;

    // curvature method:
    // CV / RDF_raw / RDF_smooth / PARABOLOID / HF / VOF / PC / PV
    std::string curv_method = "CV";

    // rdf.grad_stencil = 5pt / 9pt
    std::string rdf_grad_stencil = "9pt";

    // RDF parameters
    double rdf_A         = 25.0;
    double rdf_B         = 0.01;
    double rdf_rad_cells = 3;

    // Paraboloid/parabola fit parameters
    int parab_stencil = 1;
    int parab_min_segments = 3;

    // HF-equivalent local column half-length in cells (2*m+1) used by PARABOLOID method
    int hf_m = 3;

    // error band threshold for "interface cell"
    double eps_if = 1e-5;

    // ---- NEW: VOF init method ----
    // "recursive" (your original recursive PLIC-based init) or "exact_poly" (circle-only exact polygon clipping)
    std::string vof_init_method = "recursive";

    // ---- NEW: exact circle polygon resolution ----
    int hf_Ntheta = 720;

    // ---- NEW: PPM dump ----
    int dump_ppm = 0;

    // ---- VFM options ----
    int    vfm_stencil = 7;                 // odd, e.g. 5/7/9
    double vfm_weight_radius_factor = 3.5;  // weightRadiusD = factor * dx
    int    vfm_maxIters = 80;
    double vfm_tolRel   = 1e-12;
    double vfm_lambda0  = 1e-3;
    double vfm_lambdaUp = 10.0;
    double vfm_lambdaDown = 0.3;
    double vfm_fdEps    = 1e-6;
    // ---- PC (PLIC-centroidal) options ----
    int    pc_stencil = 7;        // odd, >=3
    double pc_d_over_dx = 3.0;    // dWeight = pc_d_over_dx * dx

};

bool load_config(const std::string& fname, Config& cfg) {
    std::ifstream fin(fname);
    if (!fin) {
        std::cerr << "Warning: cannot open input file '" << fname
                  << "', using default parameters.\n";
        return false;
    }

    std::string line;
    while (std::getline(fin, line)) {
        auto pos_comment = line.find('#');
        if (pos_comment != std::string::npos) line = line.substr(0, pos_comment);
        line = trim(line);
        if (line.empty()) continue;

        auto pos_eq = line.find('=');
        if (pos_eq == std::string::npos) continue;

        std::string key = trim(line.substr(0, pos_eq));
        std::string val = trim(line.substr(pos_eq + 1));
        if (key.empty() || val.empty()) continue;

        std::istringstream iss(val);

        if (key == "geometry.prob_lo") {
            double xlo, ylo; if (iss >> xlo >> ylo) { cfg.prob_lo_x = xlo; cfg.prob_lo_y = ylo; }
        } else if (key == "geometry.prob_hi") {
            double xhi, yhi; if (iss >> xhi >> yhi) { cfg.prob_hi_x = xhi; cfg.prob_hi_y = yhi; }
        } else if (key == "amr.n_cell") {
            int nx, ny; if (iss >> nx >> ny) { cfg.nx = nx; cfg.ny = ny; }
        } else if (key == "amr.n_list") {
            cfg.n_list.clear(); int N; while (iss >> N) cfg.n_list.push_back(N);
        } else if (key == "ic.shape") {
            cfg.ic_shape = val;
        } else if (key == "circle.R") {
            double R; if (iss >> R) cfg.R = R;
        } else if (key == "circle.center") {
            double cx, cy; if (iss >> cx >> cy) { cfg.cx = cx; cfg.cy = cy; }
        } else if (key == "cosine.A") {
            double A; if (iss >> A) cfg.cos_A = A;
        } else if (key == "cosine.L") {
            double L; if (iss >> L) cfg.cos_L = L;
        } else if (key == "cosine.n_periods") {
            int n; if (iss >> n) cfg.cos_n = n;
        } else if (key == "curv.h_factor") {
            double hf; if (iss >> hf) cfg.h_factor = hf;
        } else if (key == "vof.refine_max_level") {
            int lev; if (iss >> lev) cfg.vof_refine_max_level = lev;
        } else if (key == "curv.method") {
            cfg.curv_method = val;
        } else if (key == "rdf.grad_stencil") {
            cfg.rdf_grad_stencil = val;
        } else if (key == "rdf.A") {
            double A; if (iss >> A) cfg.rdf_A = A;
        } else if (key == "rdf.B") {
            double B; if (iss >> B) cfg.rdf_B = B;
        } else if (key == "rdf.rad_cells") {
            double r; if (iss >> r) cfg.rdf_rad_cells = r;
        } else if (key == "parab.stencil") {
            int s; if (iss >> s) cfg.parab_stencil = std::max(1, s);
        } else if (key == "parab.min_segments") {
            int m; if (iss >> m) cfg.parab_min_segments = std::max(1, m);
        } else if (key == "hf.m") {
            int m; if (iss >> m) cfg.hf_m = std::max(1, m);
        } else if (key == "error.eps_if") {
            double e; if (iss >> e) cfg.eps_if = e;
        } else if (key == "vof.init_method") {
            cfg.vof_init_method = val; // "recursive" or "exact_poly"
        } else if (key == "hf.Ntheta") {
            int nt; if (iss >> nt) cfg.hf_Ntheta = std::max(32, nt);
        } else if (key == "dump.ppm") {
            int d; if (iss >> d) cfg.dump_ppm = (d != 0);
        } else if (key == "vfm.stencil") {
            int s; if (iss >> s) cfg.vfm_stencil = std::max(3, (s % 2 == 0 ? s+1 : s));
        } else if (key == "vfm.weight_radius_factor") {
            double f; if (iss >> f) cfg.vfm_weight_radius_factor = f;
        } else if (key == "vfm.maxIters") {
            int it; if (iss >> it) cfg.vfm_maxIters = std::max(1, it);
        } else if (key == "vfm.tolRel") {
            double t; if (iss >> t) cfg.vfm_tolRel = t;
        } else if (key == "vfm.lambda0") {
            double t; if (iss >> t) cfg.vfm_lambda0 = t;
        } else if (key == "vfm.lambdaUp") {
            double t; if (iss >> t) cfg.vfm_lambdaUp = t;
        } else if (key == "vfm.lambdaDown") {
            double t; if (iss >> t) cfg.vfm_lambdaDown = t;
        } else if (key == "vfm.fdEps") {
            double t; if (iss >> t) cfg.vfm_fdEps = t;
        } else if (key == "pc.stencil") {
            int s; if (iss >> s) cfg.pc_stencil = std::max(3, (s % 2 == 0 ? s+1 : s));
        } else if (key == "pc.d_over_dx") {
            double t; if (iss >> t) cfg.pc_d_over_dx = t;
        }

    }
    return true;
}

// ====================== Geometry "inside" tests ====================
inline bool in_circle(double x, double y, const Config& cfg) {
    double dx = x - cfg.cx, dy = y - cfg.cy;
    return (dx*dx + dy*dy) <= cfg.R * cfg.R;
}
inline bool in_cosine_region(double x, double y, const Config& cfg) {
    const double two_pi = 2.0 * M_PI;
    double alpha = two_pi * cfg.cos_n / cfg.cos_L;
    double y_if  = cfg.cos_A - std::cos(alpha * x);
    return (y <= y_if);
}

// ====================== Signed level-set and gradient ==============
inline double phi_circle(double x, double y, const Config& cfg) {
    double dx = x - cfg.cx, dy = y - cfg.cy;
    double r  = std::sqrt(dx*dx + dy*dy);
    return r - cfg.R;
}
inline void grad_phi_circle(double x, double y, const Config& cfg, double& gx, double& gy) {
    double dx = x - cfg.cx, dy = y - cfg.cy;
    double r2 = dx*dx + dy*dy;
    double r  = std::sqrt(r2);
    if (r < 1e-14) { gx = 0.0; gy = 0.0; }
    else          { gx = dx / r; gy = dy / r; }
}
inline double phi_cosine(double x, double y, const Config& cfg) {
    const double two_pi = 2.0 * M_PI;
    double alpha = two_pi * cfg.cos_n / cfg.cos_L;
    double y_if  = cfg.cos_A - std::cos(alpha * x);
    return y - y_if;
}
inline void grad_phi_cosine(double x, double y, const Config& cfg, double& gx, double& gy) {
    const double two_pi = 2.0 * M_PI;
    double alpha = two_pi * cfg.cos_n / cfg.cos_L;
    gx = -alpha * std::sin(alpha * x);
    gy = 1.0;
    (void)y;
}

// ====================== Small polygon utility for PLIC =============
struct Point2D { double x, y; };

// Clip polygon with half-plane: n·x + d <= 0
std::vector<Point2D> clipPolygonWithHalfplane(const std::vector<Point2D>& poly,
                                              double nx, double ny, double d)
{
    std::vector<Point2D> out;
    if (poly.empty()) return out;

    auto eval = [&](const Point2D& p) { return nx * p.x + ny * p.y + d; };

    for (size_t i = 0; i < poly.size(); ++i) {
        const Point2D& P = poly[i];
        const Point2D& Q = poly[(i + 1) % poly.size()];

        double sP = eval(P);
        double sQ = eval(Q);

        bool insideP = (sP <= 0.0);
        bool insideQ = (sQ <= 0.0);

        if (insideP && insideQ) {
            out.push_back(Q);
        } else if (insideP && !insideQ) {
            double denom = (nx * (Q.x - P.x) + ny * (Q.y - P.y));
            if (std::fabs(denom) > 1e-14) {
                double t = - (nx * P.x + ny * P.y + d) / denom;
                Point2D I{ P.x + t * (Q.x - P.x), P.y + t * (Q.y - P.y) };
                out.push_back(I);
            }
        } else if (!insideP && insideQ) {
            double denom = (nx * (Q.x - P.x) + ny * (Q.y - P.y));
            if (std::fabs(denom) > 1e-14) {
                double t = - (nx * P.x + ny * P.y + d) / denom;
                Point2D I{ P.x + t * (Q.x - P.x), P.y + t * (Q.y - P.y) };
                out.push_back(I);
            }
            out.push_back(Q);
        }
    }

    return out;
}

double polygonArea(const std::vector<Point2D>& poly)
{
    size_t n = poly.size();
    if (n < 3) return 0.0;
    double area = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const Point2D& P = poly[i];
        const Point2D& Q = poly[(i + 1) % n];
        area += P.x * Q.y - Q.x * P.y;
    }
    return 0.5 * std::fabs(area);
}

// ====================== PLIC volume fraction in a cell =============
double volume_fraction_PLIC(double xc, double yc,
                            double dx, double dy,
                            const Config& cfg,
                            const std::string& shape,
                            int inside_count_corner)
{
    if (inside_count_corner == 4) return 1.0;
    if (inside_count_corner == 0) return 0.0;

    double phi_c = 0.0;
    double gx = 0.0, gy = 0.0;

    if (shape == "circle") {
        phi_c = phi_circle(xc, yc, cfg);
        grad_phi_circle(xc, yc, cfg, gx, gy);
    } else if (shape == "cosine") {
        phi_c = phi_cosine(xc, yc, cfg);
        grad_phi_cosine(xc, yc, cfg, gx, gy);
    } else {
        return inside_count_corner / 4.0;
    }

    double grad_norm = std::sqrt(gx*gx + gy*gy);
    if (grad_norm < 1e-14) return inside_count_corner / 4.0;

    double nx = gx / grad_norm;
    double ny = gy / grad_norm;

    // φ(x) ≈ φc + ∇φ·(x-xc) = 0 -> n·x + (φc - n·xc) = 0
    double d = phi_c - (nx * xc + ny * yc);

    double xL = xc - 0.5 * dx, xR = xc + 0.5 * dx;
    double yB = yc - 0.5 * dy, yT = yc + 0.5 * dy;

    std::vector<Point2D> rect = {{xL,yB},{xR,yB},{xR,yT},{xL,yT}};
    std::vector<Point2D> clipped = clipPolygonWithHalfplane(rect, nx, ny, d);

    double area_inside = polygonArea(clipped);
    double cell_area   = dx * dy;
    if (cell_area <= 0.0) return 0.0;

    double vf = area_inside / cell_area;
    return clamp01(vf);
}

// ====================== Recursive local refinement for VOF =========
template <typename InsideFunc>
double compute_volume_fraction_recursive(
    double xc, double yc,
    double dx, double dy,
    const Config& cfg,
    InsideFunc inside_func,
    int level,
    int maxLevel,
    const std::string& shape)
{
    double xL = xc - 0.5 * dx, xR = xc + 0.5 * dx;
    double yB = yc - 0.5 * dy, yT = yc + 0.5 * dy;

    bool c1 = inside_func(xL, yB, cfg);
    bool c2 = inside_func(xR, yB, cfg);
    bool c3 = inside_func(xL, yT, cfg);
    bool c4 = inside_func(xR, yT, cfg);

    int inside_count = (int)c1 + (int)c2 + (int)c3 + (int)c4;

    if (inside_count == 4) return 1.0;
    if (inside_count == 0) return 0.0;

    if (level >= maxLevel) {
        return volume_fraction_PLIC(xc, yc, dx, dy, cfg, shape, inside_count);
    }

    double dx2 = 0.5 * dx, dy2 = 0.5 * dy;

    double xc1 = xc - 0.25 * dx, yc1 = yc - 0.25 * dy;
    double xc2 = xc + 0.25 * dx, yc2 = yc - 0.25 * dy;
    double xc3 = xc - 0.25 * dx, yc3 = yc + 0.25 * dy;
    double xc4 = xc + 0.25 * dx, yc4 = yc + 0.25 * dy;

    double vf1 = compute_volume_fraction_recursive(xc1, yc1, dx2, dy2, cfg, inside_func, level+1, maxLevel, shape);
    double vf2 = compute_volume_fraction_recursive(xc2, yc2, dx2, dy2, cfg, inside_func, level+1, maxLevel, shape);
    double vf3 = compute_volume_fraction_recursive(xc3, yc3, dx2, dy2, cfg, inside_func, level+1, maxLevel, shape);
    double vf4 = compute_volume_fraction_recursive(xc4, yc4, dx2, dy2, cfg, inside_func, level+1, maxLevel, shape);

    return 0.25 * (vf1 + vf2 + vf3 + vf4);
}

// ====================== Initialize VOF (original) ==================
Grid2D initializeCircleVOF_recursive(int nx, int ny, const Config& cfg)
{
    double Lx = cfg.prob_hi_x - cfg.prob_lo_x;
    double Ly = cfg.prob_hi_y - cfg.prob_lo_y;
    Grid2D g(nx, ny, Lx / nx, Ly / ny);

    for (int j = 0; j < ny; ++j) {
        double y = cfg.prob_lo_y + (j + 0.5) * g.dy;
        for (int i = 0; i < nx; ++i) {
            double x = cfg.prob_lo_x + (i + 0.5) * g.dx;

            double vf = compute_volume_fraction_recursive(
                x, y, g.dx, g.dy, cfg,
                in_circle, 0, cfg.vof_refine_max_level, "circle"
            );

            g(i, j) = vf;
        }
    }
    return g;
}
Grid2D initializeCosineVOF_recursive(int nx, int ny, const Config& cfg)
{
    double Lx = cfg.prob_hi_x - cfg.prob_lo_x;
    double Ly = cfg.prob_hi_y - cfg.prob_lo_y;
    Grid2D g(nx, ny, Lx / nx, Ly / ny);

    for (int j = 0; j < ny; ++j) {
        double y = cfg.prob_lo_y + (j + 0.5) * g.dy;
        for (int i = 0; i < nx; ++i) {
            double x = cfg.prob_lo_x + (i + 0.5) * g.dx;

            double vf = compute_volume_fraction_recursive(
                x, y, g.dx, g.dy, cfg,
                in_cosine_region, 0, cfg.vof_refine_max_level, "cosine"
            );

            g(i, j) = vf;
        }
    }
    return g;
}


// ===================================================================
struct Pt { double x, y; };
using Poly = std::vector<Pt>;

static double poly_area(const Poly &p) {
    int n = (int)p.size();
    if (n < 3) return 0.0;
    double A = 0.0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        A += p[i].x * p[j].y - p[j].x * p[i].y;
    }
    return std::fabs(A) * 0.5;
}

// Sutherland–Hodgman clip against a directed edge a->b, keeping left side
static Poly clip_poly_with_edge(const Poly &subject, const Pt &a, const Pt &b) {
    Poly out;
    auto is_left = [&](const Pt &p) {
        double cross = (b.x - a.x)*(p.y - a.y) - (b.y - a.y)*(p.x - a.x);
        return cross >= -1e-12;
    };

    int n = (int)subject.size();
    if (n == 0) return out;

    for (int i = 0; i < n; ++i) {
        Pt cur  = subject[i];
        Pt prev = subject[(i + n - 1) % n];

        bool cur_in  = is_left(cur);
        bool prev_in = is_left(prev);

        auto intersect = [&](const Pt& P, const Pt& Q)->Pt {
            double A1x = Q.x - P.x, A1y = Q.y - P.y;
            double A2x = b.x - a.x, A2y = b.y - a.y;
            double denom = A1x*A2y - A1y*A2x;
            if (std::fabs(denom) < 1e-18) return Q;
            double s = ((a.x - P.x)*A2y - (a.y - P.y)*A2x) / denom;
            return { P.x + s*A1x, P.y + s*A1y };
        };

        if (cur_in) {
            if (!prev_in) out.push_back(intersect(prev, cur));
            out.push_back(cur);
        } else if (prev_in) {
            out.push_back(intersect(prev, cur));
        }
    }
    return out;
}

static Poly clip_poly_with_rect(const Poly &poly,
                                double xmin, double xmax,
                                double ymin, double ymax) {
    Poly res = poly;
    Pt e1{xmin,ymin}, e2{xmax,ymin}, e3{xmax,ymax}, e4{xmin,ymax};
    res = clip_poly_with_edge(res, e1, e2);
    res = clip_poly_with_edge(res, e2, e3);
    res = clip_poly_with_edge(res, e3, e4);
    res = clip_poly_with_edge(res, e4, e1);
    return res;
}

static Poly build_circle_poly(double cx, double cy, double R, int Ntheta) {
    Poly p; p.reserve(Ntheta);
    for (int k = 0; k < Ntheta; ++k) {
        double th = 2.0 * M_PI * k / Ntheta;
        p.push_back({ cx + R*std::cos(th), cy + R*std::sin(th) });
    }
    return p;
}

// Exact VOF init for circle: clip circle polygon against each cell rectangle
static Grid2D initializeCircleVOF_exact_poly(int nx, int ny, const Config& cfg) {
    double Lx = cfg.prob_hi_x - cfg.prob_lo_x;
    double Ly = cfg.prob_hi_y - cfg.prob_lo_y;
    Grid2D g(nx, ny, Lx / nx, Ly / ny);

    const double dx = g.dx;
    const double dy = g.dy;
    const double cell_area = dx * dy;

    Poly circle_poly = build_circle_poly(cfg.cx, cfg.cy, cfg.R, cfg.hf_Ntheta);

    for (int j = 0; j < ny; ++j) {
        double ymin = cfg.prob_lo_y + j * dy;
        double ymax = ymin + dy;
        for (int i = 0; i < nx; ++i) {
            double xmin = cfg.prob_lo_x + i * dx;
            double xmax = xmin + dx;

            Poly clipped = clip_poly_with_rect(circle_poly, xmin, xmax, ymin, ymax);
            double vf = (cell_area > 0.0) ? (poly_area(clipped) / cell_area) : 0.0;
            g(i, j) = clamp01(vf);
        }
    }
    return g;
}

// ====================== CV curvature from smoothed VOF =============
Grid2D computeCurvatureCV(const Grid2D& ft) {
    Grid2D kappa(ft.nx, ft.ny, ft.dx, ft.dy);

    const int nx = ft.nx;
    const int ny = ft.ny;
    const double dx = ft.dx;
    const double dy = ft.dy;

    Grid2D nx_field(nx, ny, dx, dy);
    Grid2D ny_field(nx, ny, dx, dy);

    const double eps = 1e-12;

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            double fx = (ft(i + 1, j) - ft(i - 1, j)) / (2.0 * dx);
            double fy = (ft(i, j + 1) - ft(i, j - 1)) / (2.0 * dy);

            double mag = std::sqrt(fx * fx + fy * fy);
            if (mag < eps) {
                nx_field(i, j) = 0.0;
                ny_field(i, j) = 0.0;
            } else {
                nx_field(i, j) = -fx / mag;
                ny_field(i, j) = -fy / mag;
            }
        }
    }

    for (int i = 0; i < nx; ++i) {
        nx_field(i, 0)      = nx_field(i, 1);
        ny_field(i, 0)      = ny_field(i, 1);
        nx_field(i, ny - 1) = nx_field(i, ny - 2);
        ny_field(i, ny - 1) = ny_field(i, ny - 2);
    }
    for (int j = 0; j < ny; ++j) {
        nx_field(0,      j) = nx_field(1,      j);
        ny_field(0,      j) = ny_field(1,      j);
        nx_field(nx - 1, j) = nx_field(nx - 2, j);
        ny_field(nx - 1, j) = ny_field(nx - 2, j);
    }

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            double dnx_dx = (nx_field(i + 1, j) - nx_field(i - 1, j)) / (2.0 * dx);
            double dny_dy = (ny_field(i, j + 1) - ny_field(i, j - 1)) / (2.0 * dy);
            kappa(i, j) = dnx_dx + dny_dy;
        }
    }

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double ft_val = ft(i, j);
            if (ft_val <= 1e-3 || ft_val >= 1.0 - 1e-3) kappa(i, j) = 0.0;
        }
    }

    return kappa;
}

// ====================== RDF helpers ================================
struct InterfaceCell {
    int i, j;
    double xc, yc;     // physical coordinates (cell center)
    double fg;         // VOF
    double nx, ny;     // outward unit normal (from grad f)
    double d_line;     // line: n·x + d = 0, inside: <=0
    double xint, yint; // a point on reconstructed interface line (projection)
};

// 5-point gradient
inline void gradient_5pt(const Grid2D& vof, int i, int j, double dx, double dy,
                         double& fx, double& fy)
{
    fx = (vof(i + 1, j)   - vof(i - 1, j))   / (2.0 * dx);
    fy = (vof(i,     j+1) - vof(i,     j-1)) / (2.0 * dy);
}

// 9-point gradient (Sobel-like)
inline void gradient_9pt(const Grid2D& vof, int i, int j, double dx, double dy,
                         double& fx, double& fy)
{
    double f_ip1_jp1 = vof(i+1, j+1);
    double f_ip1_j   = vof(i+1, j);
    double f_ip1_jm1 = vof(i+1, j-1);
    double f_im1_jp1 = vof(i-1, j+1);
    double f_im1_j   = vof(i-1, j);
    double f_im1_jm1 = vof(i-1, j-1);
    double f_i_jp1   = vof(i,   j+1);
    double f_i_jm1   = vof(i,   j-1);

    fx = ( f_ip1_jp1 + 2.0*f_ip1_j + f_ip1_jm1
         - f_im1_jp1 - 2.0*f_im1_j - f_im1_jm1 ) / (8.0 * dx);

    fy = ( f_ip1_jp1 + 2.0*f_i_jp1 + f_im1_jp1
         - f_ip1_jm1 - 2.0*f_i_jm1 - f_im1_jm1 ) / (8.0 * dy);
}

// Compute area fraction for rectangle cut by n·x + d <= 0
static double vf_from_line(double xc, double yc, double dx, double dy,
                           double nx, double ny, double d_line)
{
    double xL = xc - 0.5 * dx, xR = xc + 0.5 * dx;
    double yB = yc - 0.5 * dy, yT = yc + 0.5 * dy;

    std::vector<Point2D> rect = {{xL,yB},{xR,yB},{xR,yT},{xL,yT}};
    auto clipped = clipPolygonWithHalfplane(rect, nx, ny, d_line);
    double a = polygonArea(clipped);
    return a / (dx * dy);
}

// Solve d_line so that area fraction equals f_target (bisection)
static double solve_d_for_vf(double xc, double yc, double dx, double dy,
                            double nx, double ny, double f_target)
{
    double L = 0.5 * (std::fabs(nx)*dx + std::fabs(ny)*dy);

    double base = -(nx*xc + ny*yc);
    double lo = base - L;
    double hi = base + L;

    double f_lo = vf_from_line(xc,yc,dx,dy,nx,ny,lo);
    double f_hi = vf_from_line(xc,yc,dx,dy,nx,ny,hi);

    if (f_lo < f_hi) { std::swap(lo, hi); std::swap(f_lo, f_hi); }

    f_target = std::min(f_lo, std::max(f_hi, f_target));

    for (int it = 0; it < 60; ++it) {
        double mid = 0.5*(lo+hi);
        double f_mid = vf_from_line(xc,yc,dx,dy,nx,ny,mid);
        if (f_mid > f_target) lo = mid;
        else                  hi = mid;
    }
    return 0.5*(lo+hi);
}

// Build interface cell list for RDF (with PLIC line reconstruction)
std::vector<InterfaceCell> buildInterfaceCells(const Grid2D& vof,
                                               const Config& cfg)
{
    std::vector<InterfaceCell> cells;

    const int nx = vof.nx;
    const int ny = vof.ny;
    const double dx = vof.dx;
    const double dy = vof.dy;

    const double B  = cfg.rdf_B;
    const double grad_min = B / std::min(dx, dy);

    bool use9 = (cfg.rdf_grad_stencil == "9pt" ||
                 cfg.rdf_grad_stencil == "nine" ||
                 cfg.rdf_grad_stencil == "NINE");

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            double f = vof(i, j);
            if (f <= 1e-12 || f >= 1.0 - 1e-12) continue;

            double fx, fy;
            if (use9) gradient_9pt(vof, i, j, dx, dy, fx, fy);
            else      gradient_5pt(vof, i, j, dx, dy, fx, fy);

            double grad = std::sqrt(fx * fx + fy * fy);
            if (grad < grad_min) continue;

            InterfaceCell c;
            c.i  = i;
            c.j  = j;
            c.xc = cfg.prob_lo_x + (i + 0.5) * dx;
            c.yc = cfg.prob_lo_y + (j + 0.5) * dy;
            c.fg = f;

            c.nx = -fx / grad;
            c.ny = -fy / grad;

            c.d_line = solve_d_for_vf(c.xc, c.yc, dx, dy, c.nx, c.ny, c.fg);

            double s0 = c.nx*c.xc + c.ny*c.yc + c.d_line;
            c.xint = c.xc - s0 * c.nx;
            c.yint = c.yc - s0 * c.ny;

            cells.push_back(c);
        }
    }

    return cells;
}

// ====================== RDF: build φ from VOF ======================
Grid2D buildRDF(const Grid2D& vof, const Config& cfg)
{
    Grid2D phi(vof.nx, vof.ny, vof.dx, vof.dy);

    const int nx = vof.nx;
    const int ny = vof.ny;
    const double dx = vof.dx;
    const double dy = vof.dy;

    std::vector<InterfaceCell> iface = buildInterfaceCells(vof, cfg);
    if (iface.empty()) {
        std::cerr << "Warning: no interface cells found for RDF.\n";
        return phi;
    }

    const double Aexp = cfg.rdf_A;
    const int rad_i = std::max(1, (int)std::ceil(cfg.rdf_rad_cells));
    const int rad_j = std::max(1, (int)std::ceil(cfg.rdf_rad_cells));
    const double Rphys = cfg.rdf_rad_cells * std::min(dx, dy);

    for (int j = 0; j < ny; ++j) {
        double yij = cfg.prob_lo_y + (j + 0.5) * dy;
        for (int i = 0; i < nx; ++i) {
            double xij = cfg.prob_lo_x + (i + 0.5) * dx;

            double sum_WD = 0.0;
            double sum_W  = 0.0;

            for (const auto& g : iface) {
                if (std::abs(g.i - i) > rad_i) continue;
                if (std::abs(g.j - j) > rad_j) continue;

                double Dgij = g.nx * xij + g.ny * yij + g.d_line;
                double r = std::hypot(xij - g.xint, yij - g.yint);
                if (r < 1e-14) continue;

                double cos_theta = std::fabs(Dgij) / r;
                cos_theta = std::min(1.0, std::max(0.0, cos_theta));
                double w_dir = std::pow(cos_theta, Aexp);

                double q = r / (Rphys + 1e-30);
                double w_r = kernel(q);
                if (w_r == 0.0) continue;

                double Wgij = g.fg * (1.0 - g.fg) * w_dir * w_r;

                sum_WD += Wgij * Dgij;
                sum_W  += Wgij;
            }

            phi(i, j) = (sum_W > 0.0) ? (sum_WD / sum_W) : 0.0;
        }
    }

    return phi;
}

// ====================== RDF curvature from φ =======================
Grid2D computeCurvatureRDF(const Grid2D& phi)
{
    Grid2D kappa(phi.nx, phi.ny, phi.dx, phi.dy);

    const int nx = phi.nx;
    const int ny = phi.ny;
    const double dx = phi.dx;
    const double dy = phi.dy;

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            double phix = (phi(i + 1, j) - phi(i - 1, j)) / (2.0 * dx);
            double phiy = (phi(i, j + 1) - phi(i, j - 1)) / (2.0 * dy);

            double phixx = (phi(i + 1, j) - 2.0 * phi(i, j) + phi(i - 1, j)) / (dx * dx);
            double phiyy = (phi(i, j + 1) - 2.0 * phi(i, j) + phi(i, j - 1)) / (dy * dy);
            double phixy = (phi(i + 1, j + 1) - phi(i + 1, j - 1)
                          - phi(i - 1, j + 1) + phi(i - 1, j - 1)) / (4.0 * dx * dy);

            double nx_ = phix;
            double ny_ = phiy;

            double grad2 = nx_ * nx_ + ny_ * ny_;
            if (grad2 < 1e-20) { kappa(i, j) = 0.0; continue; }

            double num = phixx * ny_ * ny_ + phiyy * nx_ * nx_
                       - 2.0 * nx_ * ny_ * phixy;
            double den = std::pow(grad2, 1.5);

            kappa(i, j) = num / den;
        }
    }

    return kappa;
}

// ================================================================
//     Parabolic
// ================================================================

// line-rectangle intersection
static std::vector<Point2D> intersectLineWithRect(double nx,double ny,double d,
                                                  double xL,double xR,double yB,double yT)
{
    std::vector<Point2D> pts;
    auto add_unique = [&](double x,double y){
        for(const auto& p:pts){
            if(std::hypot(p.x-x,p.y-y)<1e-12) return;
        }
        pts.push_back({x,y});
    };

    if (std::fabs(ny) > 1e-14) {
        double y = -(nx*xL + d)/ny;
        if (y >= yB-1e-14 && y <= yT+1e-14) add_unique(xL,y);
        y = -(nx*xR + d)/ny;
        if (y >= yB-1e-14 && y <= yT+1e-14) add_unique(xR,y);
    }
    if (std::fabs(nx) > 1e-14) {
        double x = -(ny*yB + d)/nx;
        if (x >= xL-1e-14 && x <= xR+1e-14) add_unique(x,yB);
        x = -(ny*yT + d)/nx;
        if (x >= xL-1e-14 && x <= xR+1e-14) add_unique(x,yT);
    }

    if (pts.size() > 2) {
        double best = -1.0;
        int a=0,b=1;
        for (int i=0;i<(int)pts.size();++i){
            for(int j=i+1;j<(int)pts.size();++j){
                double dist=std::hypot(pts[i].x-pts[j].x, pts[i].y-pts[j].y);
                if(dist>best){best=dist;a=i;b=j;}
            }
        }
        return {pts[a], pts[b]};
    }

    return pts;
}

// column height used by paraboloid method
static double columnHeightHF(const Grid2D& vof,
                             int i, int j,
                             bool column_in_y,
                             int half_span)
{
    const int nx = vof.nx, ny = vof.ny;
    const double dx = vof.dx, dy = vof.dy;
    (void)dx;

    double H = 0.0;
    if (column_in_y) {
        for (int k = -half_span; k <= half_span; ++k) {
            int jj = j + k;
            if (jj < 0 || jj >= ny) continue;
            H += clamp01(vof(i, jj)) * dy;
        }
    } else {
        for (int k = -half_span; k <= half_span; ++k) {
            int ii = i + k;
            if (ii < 0 || ii >= nx) continue;
            H += clamp01(vof(ii, j)) * vof.dx;
        }
    }
    return H;
}

Grid2D computeCurvatureParabolicReconstructionPaperCartesian(const Grid2D& vof,
                                                             const Config& cfg)
{
    Grid2D kappa(vof.nx, vof.ny, vof.dx, vof.dy);

    std::vector<InterfaceCell> iface = buildInterfaceCells(vof, cfg);
    if (iface.empty()) return kappa;

    const int nx = vof.nx;
    const int ny = vof.ny;
    const double dx = vof.dx;
    const double dy = vof.dy;

    const int half_span = cfg.hf_m;

    for (const auto& c0 : iface) {
        const int i0 = c0.i;
        const int j0 = c0.j;

        bool tangent_in_x = (std::fabs(c0.ny) >= std::fabs(c0.nx));

        const int im1 = i0 - 1, ip1 = i0 + 1;
        const int jm1 = j0 - 1, jp1 = j0 + 1;

        if (tangent_in_x) {
            if (im1 < 0 || ip1 >= nx) { kappa(i0, j0) = 0.0; continue; }
        } else {
            if (jm1 < 0 || jp1 >= ny) { kappa(i0, j0) = 0.0; continue; }
        }

        double w = 0.0;
        double Hm1 = 0.0, H0 = 0.0, Hp1 = 0.0;

        if (tangent_in_x) {
            w   = dx;
            Hm1 = columnHeightHF(vof, im1, j0, true,  half_span);
            H0  = columnHeightHF(vof, i0,  j0, true,  half_span);
            Hp1 = columnHeightHF(vof, ip1, j0, true,  half_span);
        } else {
            w   = dy;
            Hm1 = columnHeightHF(vof, i0, jm1, false, half_span);
            H0  = columnHeightHF(vof, i0, j0,  false, half_span);
            Hp1 = columnHeightHF(vof, i0, jp1, false, half_span);
        }

        const double a = (Hm1 - 2.0 * H0 + Hp1) / (2.0 * w * w);
        const double b = (Hp1 - Hm1) / (2.0 * w);

        const double denom = std::pow(1.0 + b * b, 1.5);
        if (denom < 1e-30) { kappa(i0, j0) = 0.0; continue; }

        const double k = (2.0 * a) / denom;

        // sign convention kept
        kappa(i0, j0) = -k;
    }

    return kappa;
}

// ===================================================================
//  HF (Cummins et al. 2005) on Cartesian grid, using local height functions
//  - Local stencil: 3 columns/rows around interface cell, each height is a column-sum over (2*m+1) cells
//  - Derivatives by central differences:
//      Yx  = (Y_{i+1}-Y_{i-1})/(2 dx)
//      Yxx = (Y_{i+1}-2Y_i+Y_{i-1})/(dx^2)
//  - Curvature (Eq. 15):
//      kappa = sign(n_y) * Yxx / (1 + Yx^2)^(3/2)   when using y = Y(x)
//    For x = X(y) form, symmetric:
//      kappa = -sign(n_x) * Xyy / (1 + Xy^2)^(3/2)
// ===================================================================

static inline double global_sign_from_circle_outward(const InterfaceCell& c,
                                                     const Config& cfg)
{
    
    double rx = c.xc - cfg.cx;
    double ry = c.yc - cfg.cy;

    
    double dot = c.nx * rx + c.ny * ry;
    return (dot >= 0.0) ? 1.0 : -1.0;
}


static inline double global_sign_fallback(const InterfaceCell& c)
{
    
    (void)c;
    return 1.0;
}


static inline double sgn(double a) { return (a >= 0.0) ? 1.0 : -1.0; }

static inline bool height_in_cell_Y(double Y, double yB, double yT) {
    return (Y >= yB) && (Y <= yT);
}

static inline bool height_in_cell_X(double X, double xL, double xR) {
    return (X >= xL) && (X <= xR);
}


static Grid2D computeCurvatureHF_CosineFullColumnY(const Grid2D& vof, const Config& cfg)
{
    Grid2D kappa(vof.nx, vof.ny, vof.dx, vof.dy);

    const int nx = vof.nx;
    const int ny = vof.ny;
    const double dx = vof.dx;
    const double dy = vof.dy;
    const double eps = cfg.eps_if;

    
    std::vector<double> Y(nx, 0.0);
    for (int i = 0; i < nx; ++i) {
        double H = 0.0; 
        for (int j = 0; j < ny; ++j) {
            H += clamp01(vof(i, j)) * dy;
        }
        Y[i] = cfg.prob_lo_y + H; 
    }

   
    std::vector<double> K(nx, 0.0);

    auto safe_denom = [](double Yx) {
        double v = 1.0 + Yx * Yx;
        return std::pow(v, 1.5);
    };

    if (nx >= 4) {
        // i = 0 (one-sided)
        {
            double Yx  = (-3.0*Y[0] + 4.0*Y[1] - 1.0*Y[2]) / (2.0*dx);
            double Yxx = ( 2.0*Y[0] - 5.0*Y[1] + 4.0*Y[2] - 1.0*Y[3]) / (dx*dx);
            double denom = safe_denom(Yx);
            K[0] = (denom > 1e-30) ? (-Yxx / denom) : 0.0;
        }
        // i = nx-1 (one-sided)
        {
            int n = nx - 1;
            double Yx  = ( 3.0*Y[n] - 4.0*Y[n-1] + 1.0*Y[n-2]) / (2.0*dx);
            double Yxx = ( 2.0*Y[n] - 5.0*Y[n-1] + 4.0*Y[n-2] - 1.0*Y[n-3]) / (dx*dx);
            double denom = safe_denom(Yx);
            K[n] = (denom > 1e-30) ? (-Yxx / denom) : 0.0;
        }
        // interior (central)
        for (int i = 1; i <= nx - 2; ++i) {
            double Yx  = (Y[i+1] - Y[i-1]) / (2.0*dx);
            double Yxx = (Y[i+1] - 2.0*Y[i] + Y[i-1]) / (dx*dx);
            double denom = safe_denom(Yx);
            K[i] = (denom > 1e-30) ? (-Yxx / denom) : 0.0;
        }
    } else {
        // too small, do nothing (won't happen in your convergence list)
        return kappa;
    }

    // 3) write kappa only on interface cells
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double f = vof(i, j);
            if (f > eps && f < 1.0 - eps) kappa(i, j) = K[i];
            else                          kappa(i, j) = 0.0;
        }
    }

    return kappa;
}

static Grid2D computeCurvatureHF(const Grid2D& vof, const Config& cfg)
{
    // cosine: keep your stable full-column HF
    if (cfg.ic_shape == "cosine") {
        return computeCurvatureHF_CosineFullColumnY(vof, cfg);
    }

    Grid2D kappa(vof.nx, vof.ny, vof.dx, vof.dy);

    const int nx = vof.nx;
    const int ny = vof.ny;
    const double dx = vof.dx;
    const double dy = vof.dy;
    const int m = cfg.hf_m;
    const double eps = cfg.eps_if;

  
    const bool use9 = (cfg.rdf_grad_stencil == "9pt" ||
                       cfg.rdf_grad_stencil == "nine" ||
                       cfg.rdf_grad_stencil == "NINE");


    const double grad_min = 1e-14;

  
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {

            const double f = vof(i, j);
            if (!(f > eps && f < 1.0 - eps)) continue;

         
            double fx=0.0, fy=0.0;
            if (use9) gradient_9pt(vof, i, j, dx, dy, fx, fy);
            else      gradient_5pt(vof, i, j, dx, dy, fx, fy);

            const double g = std::sqrt(fx*fx + fy*fy);
            if (g < grad_min) continue;

            const double nx0 = -fx / g;
            const double ny0 = -fy / g;

        
            const bool use_Y_of_x = (std::fabs(ny0) >= std::fabs(nx0));

            double Ym1=0.0, Y0=0.0, Yp1=0.0;
            double Xm1=0.0, X0=0.0, Xp1=0.0;
            double k_local = 0.0;

            if (use_Y_of_x) {
                if (i - 1 < 0 || i + 1 >= nx) continue;

                Ym1 = columnHeightHF(vof, i - 1, j, /*column_in_y=*/true,  m);
                Y0  = columnHeightHF(vof, i,     j, /*column_in_y=*/true,  m);
                Yp1 = columnHeightHF(vof, i + 1, j, /*column_in_y=*/true,  m);

                const double Yx  = (Yp1 - Ym1) / (2.0 * dx);
                const double Yxx = (Yp1 - 2.0 * Y0 + Ym1) / (dx * dx);

                const double denom = std::pow(1.0 + Yx*Yx, 1.5);
                if (denom < 1e-30) continue;

                k_local = Yxx / denom;

            
                if (cfg.ic_shape == "circle") {
               
                    InterfaceCell tmp{};
                    tmp.xc = cfg.prob_lo_x + (i + 0.5) * dx;
                    tmp.yc = cfg.prob_lo_y + (j + 0.5) * dy;
                    tmp.nx = nx0;
                    tmp.ny = ny0;
                    const double sgn_global = global_sign_from_circle_outward(tmp, cfg);
                    kappa(i, j) = sgn_global * std::fabs(k_local);
                } else {
                    kappa(i, j) = sgn(ny0) * k_local;
                }

            } else {
                if (j - 1 < 0 || j + 1 >= ny) continue;

                Xm1 = columnHeightHF(vof, i, j - 1, /*column_in_y=*/false, m);
                X0  = columnHeightHF(vof, i, j,     /*column_in_y=*/false, m);
                Xp1 = columnHeightHF(vof, i, j + 1, /*column_in_y=*/false, m);

                const double Xy  = (Xp1 - Xm1) / (2.0 * dy);
                const double Xyy = (Xp1 - 2.0 * X0 + Xm1) / (dy * dy);

                const double denom = std::pow(1.0 + Xy*Xy, 1.5);
                if (denom < 1e-30) continue;

                k_local = Xyy / denom;

                // Cummins sign for x=X(y)
                if (cfg.ic_shape == "circle") {
                    InterfaceCell tmp{};
                    tmp.xc = cfg.prob_lo_x + (i + 0.5) * dx;
                    tmp.yc = cfg.prob_lo_y + (j + 0.5) * dy;
                    tmp.nx = nx0;
                    tmp.ny = ny0;
                    const double sgn_global = global_sign_from_circle_outward(tmp, cfg);
                    kappa(i, j) = sgn_global * std::fabs(k_local);
                } else {
                    kappa(i, j) = -sgn(nx0) * k_local;
                }
            }
        }
    }

    return kappa;
}


// ===================================================================
//  VFM (Volume Fraction Matching) curvature: parabola fit in local (t,n) frame
// ===================================================================
namespace vfm {

struct Vec2 {
    double x=0, y=0;
    Vec2()=default;
    Vec2(double X,double Y):x(X),y(Y){}
    Vec2 operator+(const Vec2& o) const { return {x+o.x, y+o.y}; }
    Vec2 operator-(const Vec2& o) const { return {x-o.x, y-o.y}; }
    Vec2 operator*(double s) const { return {x*s, y*s}; }
};

static inline double dot(const Vec2& a,const Vec2& b){ return a.x*b.x + a.y*b.y; }
static inline double norm(const Vec2& a){ return std::sqrt(dot(a,a)); }

static inline Vec2 normalize(const Vec2& a){
    const double n = norm(a);
    if(n < 1e-14) return {0,1};
    return a*(1.0/n);
}


struct Frame2 {
    Vec2 t; // tangent
    Vec2 n; // normal (outward)
    Vec2 toLocal(const Vec2& v) const { return { dot(t,v), dot(n,v) }; }
};

static inline Frame2 makeFrameFromNormal(const Vec2& nHat){
    Vec2 n = normalize(nHat);
    Vec2 t = Vec2{-n.y, n.x};
    return {t, n};
}


static double polygonAreaLocal(const std::vector<Vec2>& P){
    double A=0.0;
    for(size_t i=0;i<P.size();i++){
        const Vec2& a = P[i];
        const Vec2& b = P[(i+1)%P.size()];
        A += a.x*b.y - a.y*b.x;
    }
    return 0.5*std::fabs(A);
}

static bool yRangeAtX(const std::vector<Vec2>& P, double x, double& yMin, double& yMax){
    std::vector<double> ys;
    ys.reserve(8);

    for(size_t i=0;i<P.size();i++){
        Vec2 a = P[i];
        Vec2 b = P[(i+1)%P.size()];

        const double x1=a.x, x2=b.x;
        if((x < std::min(x1,x2)) || (x > std::max(x1,x2))) continue;

        if(std::fabs(x2-x1) < 1e-14){
            if(std::fabs(x - x1) < 1e-14){
                ys.push_back(a.y);
                ys.push_back(b.y);
            }
            continue;
        }

        const double t = (x - x1)/(x2-x1);
        if(t < -1e-12 || t > 1.0+1e-12) continue;
        const double y = a.y + t*(b.y - a.y);
        ys.push_back(y);
    }

    if(ys.size() < 2) return false;
    std::sort(ys.begin(), ys.end());
    yMin = ys.front();
    yMax = ys.back();
    return (yMax > yMin);
}

static inline double wendlandWeight(double r, double d){
    if(r > d) return 0.0;
    const double q = r/d;
    const double t = 1.0 - q;
    return (1.0 + 4.0*q) * t*t*t*t;
}


static inline double f_parab(double x, const std::array<double,3>& a){
    return a[0] + a[1]*x + a[2]*x*x;
}


static inline double curvatureFromParabolaAtOrigin(const std::array<double,3>& a){
    const double fp = a[1];
    const double fpp = 2.0*a[2];
    return fpp / std::pow(1.0 + fp*fp, 1.5);
}

static constexpr int NG = 16;
static constexpr double xg16[NG] = {
  -0.989400934991649932596154173450,
  -0.944575023073232576077988415535,
  -0.865631202387831743880467897712,
  -0.755404408355003033895101194847,
  -0.617876244402643748446671764049,
  -0.458016777657227386342419442984,
  -0.281603550779258913230460501460,
  -0.095012509837637440185319335425,
   0.095012509837637440185319335425,
   0.281603550779258913230460501460,
   0.458016777657227386342419442984,
   0.617876244402643748446671764049,
   0.755404408355003033895101194847,
   0.865631202387831743880467897712,
   0.944575023073232576077988415535,
   0.989400934991649932596154173450
};
static constexpr double wg16[NG] = {
  0.027152459411754094851780572456,
  0.062253523938647892862843836994,
  0.095158511682492784809925107602,
  0.124628971255533872052476282192,
  0.149595988816576732081501730547,
  0.169156519395002538189312079030,
  0.182603415044923588866763667969,
  0.189450610455068496285396723208,
  0.189450610455068496285396723208,
  0.182603415044923588866763667969,
  0.169156519395002538189312079030,
  0.149595988816576732081501730547,
  0.124628971255533872052476282192,
  0.095158511682492784809925107602,
  0.062253523938647892862843836994,
  0.027152459411754094851780572456
};

struct Cell2D {
    double xL,xR,yB,yT;
    double alpha;
    Vec2 center;
};

static inline double alphaTildeCell(const Cell2D& c,
                                   const Vec2& x0World,
                                   const Frame2& F,
                                   const std::array<double,3>& a)
{
    std::vector<Vec2> Pw = {
        {c.xL, c.yB}, {c.xR, c.yB}, {c.xR, c.yT}, {c.xL, c.yT}
    };

    std::vector<Vec2> P;
    P.reserve(4);
    for(const auto& w : Pw){
        Vec2 v = w - x0World;
        Vec2 l = F.toLocal(v);
        P.push_back(l);
    }

    const double polyA = polygonAreaLocal(P);
    if(polyA < 1e-20) return 0.0;

    double xMin = +std::numeric_limits<double>::infinity();
    double xMax = -std::numeric_limits<double>::infinity();
    for(const auto& p : P){
        xMin = std::min(xMin, p.x);
        xMax = std::max(xMax, p.x);
    }
    if(!(xMax > xMin)) return 0.0;

    const double Jx = 0.5*(xMax-xMin);
    const double cx = 0.5*(xMax+xMin);

    double areaUnder = 0.0;
    for(int i=0;i<NG;i++){
        const double x = cx + Jx*xg16[i];
        const double w = wg16[i];

        double yLo=0, yHi=0;
        if(!yRangeAtX(P, x, yLo, yHi)) continue;

        const double fy = f_parab(x, a);

        const double len = std::clamp(fy - yLo, 0.0, yHi - yLo);
        areaUnder += w * len;
    }
    areaUnder *= Jx;

    double at = areaUnder / polyA;
    return std::clamp(at, 0.0, 1.0);
}

static std::array<double,3> solve3(std::array<std::array<double,3>,3> A,
                                  std::array<double,3> b)
{
    for(int k=0;k<3;k++){
        int piv=k;
        double best=std::fabs(A[k][k]);
        for(int r=k+1;r<3;r++){
            if(std::fabs(A[r][k]) > best){ best=std::fabs(A[r][k]); piv=r; }
        }
        if(best < 1e-18) throw std::runtime_error("Singular 3x3 in LM.");
        if(piv!=k){ std::swap(A[piv],A[k]); std::swap(b[piv],b[k]); }

        double inv = 1.0/A[k][k];
        for(int j=k;j<3;j++) A[k][j] *= inv;
        b[k] *= inv;

        for(int i=k+1;i<3;i++){
            double f = A[i][k];
            for(int j=k;j<3;j++) A[i][j] -= f*A[k][j];
            b[i] -= f*b[k];
        }
    }

    std::array<double,3> x{};
    for(int i=2;i>=0;i--){
        double s=b[i];
        for(int j=i+1;j<3;j++) s -= A[i][j]*x[j];
        x[i]=s;
    }
    return x;
}

struct VFM2DOptions {
    int maxIters = 80;
    double tolRel = 1e-12;
    double lambda0 = 1e-3;
    double lambdaUp = 10.0;
    double lambdaDown = 0.3;
    double fdEps = 1e-6;
    double weightRadiusD = 0.02;
};

struct VFM2DResult {
    std::array<double,3> a;
    double curvature;
    double cost;
};

static VFM2DResult fitVFM2D(const std::vector<Cell2D>& stencil,
                            const Vec2& x0World,
                            const Vec2& nHatWorld,
                            const std::array<double,3>& aInit,
                            const VFM2DOptions& opt)
{
    Frame2 F = makeFrameFromNormal(nHatWorld);

    auto costFn = [&](const std::array<double,3>& a)->double{
        double J=0.0;
        for(const auto& c : stencil){
            double r = norm(c.center - x0World);
            double w = wendlandWeight(r, opt.weightRadiusD);
            if(w==0.0) continue;
            const double at = alphaTildeCell(c, x0World, F, a);
            const double d  = (at - c.alpha);
            J += w*d*d;
        }
        return J;
    };

    std::array<double,3> a = aInit;
    double lambda = opt.lambda0;
    double cost = costFn(a);

    for(int iter=0; iter<opt.maxIters; iter++){
        std::vector<std::array<double,3>> Jrows;
        std::vector<double> rvec;
        Jrows.reserve(stencil.size());
        rvec.reserve(stencil.size());

        for(const auto& c : stencil){
            double r = norm(c.center - x0World);
            double w = wendlandWeight(r, opt.weightRadiusD);
            if(w==0.0) continue;

            const double sw = std::sqrt(w);
            const double base = alphaTildeCell(c, x0World, F, a);
            rvec.push_back(sw*(base - c.alpha));

            std::array<double,3> Ji{};
            for(int j=0;j<3;j++){
                auto ap=a, am=a;
                double h = opt.fdEps*(1.0 + std::fabs(a[j]));
                ap[j]+=h; am[j]-=h;
                const double fp = alphaTildeCell(c, x0World, F, ap);
                const double fm = alphaTildeCell(c, x0World, F, am);
                Ji[j] = sw*(fp - fm)/(2.0*h);
            }
            Jrows.push_back(Ji);
        }

        std::array<std::array<double,3>,3> A{};
        std::array<double,3> g{}; // J^T r
        for(size_t i=0;i<rvec.size();i++){
            const auto& Ji = Jrows[i];
            const double ri = rvec[i];
            for(int p=0;p<3;p++){
                g[p] += Ji[p]*ri;
                for(int q=0;q<3;q++) A[p][q] += Ji[p]*Ji[q];
            }
        }

        for(int d=0; d<3; d++) A[d][d] *= (1.0 + lambda);

        std::array<double,3> rhs{-g[0],-g[1],-g[2]};
        std::array<double,3> delta = solve3(A, rhs);

        auto cand = a;
        for(int j=0;j<3;j++) cand[j] += delta[j];
        const double costCand = costFn(cand);

        if(costCand < cost){
            const double relImprove = (cost - costCand)/(cost + 1e-30);
            a = cand;
            cost = costCand;
            lambda *= opt.lambdaDown;
            if(relImprove < opt.tolRel) break;
        } else {
            lambda *= opt.lambdaUp;
        }
    }

    return {a, curvatureFromParabolaAtOrigin(a), cost};
}

} // namespace vfm

static inline vfm::Vec2 analyticNormalAt(const Config& cfg, double x, double y)
{
    double gx=0.0, gy=0.0;
    if (cfg.ic_shape == "circle") grad_phi_circle(x,y,cfg,gx,gy);
    else                          grad_phi_cosine(x,y,cfg,gx,gy);
    double n = std::sqrt(gx*gx + gy*gy);
    if (n < 1e-14) return vfm::Vec2{0,1};
    return vfm::Vec2{gx/n, gy/n}; // outward (phi increasing)
}

static Grid2D computeCurvatureVFM(const Grid2D& vof, const Config& cfg)
{
    Grid2D kappa(vof.nx, vof.ny, vof.dx, vof.dy);
    auto iface = buildInterfaceCells(vof, cfg);
    if (iface.empty()) return kappa;

    const int S = std::max(3, (cfg.vfm_stencil % 2 == 0 ? cfg.vfm_stencil+1 : cfg.vfm_stencil));
    const int hw = (S-1)/2;

    vfm::VFM2DOptions opt;
    opt.maxIters = cfg.vfm_maxIters;
    opt.tolRel   = cfg.vfm_tolRel;
    opt.lambda0  = cfg.vfm_lambda0;
    opt.lambdaUp = cfg.vfm_lambdaUp;
    opt.lambdaDown = cfg.vfm_lambdaDown;
    opt.fdEps    = cfg.vfm_fdEps;
    opt.weightRadiusD = cfg.vfm_weight_radius_factor * vof.dx;

    const std::array<double,3> aInit{0,0,0};

    for (const auto& c : iface) {
        if (!(c.fg > cfg.eps_if && c.fg < 1.0 - cfg.eps_if)) continue;

        const int i0 = c.i;
        const int j0 = c.j;

        vfm::Vec2 x0{c.xc, c.yc};
        vfm::Vec2 nHat = analyticNormalAt(cfg, c.xc, c.yc);

        std::vector<vfm::Cell2D> st;
        st.reserve(S*S);

        for(int dj=-hw; dj<=hw; dj++){
            for(int di=-hw; di<=hw; di++){
                int i = i0 + di;
                int j = j0 + dj;
                if(i<0||i>=vof.nx||j<0||j>=vof.ny) continue;

                double xL = cfg.prob_lo_x + i * vof.dx;
                double xR = xL + vof.dx;
                double yB = cfg.prob_lo_y + j * vof.dy;
                double yT = yB + vof.dy;
                double xc = xL + 0.5*vof.dx;
                double yc = yB + 0.5*vof.dy;

                vfm::Cell2D cc;
                cc.xL=xL; cc.xR=xR; cc.yB=yB; cc.yT=yT;
                cc.alpha = clamp01(vof(i,j));
                cc.center = vfm::Vec2{xc,yc};
                st.push_back(cc);
            }
        }

        if ((int)st.size() < 6) { kappa(i0,j0) = 0.0; continue; }

        auto res = vfm::fitVFM2D(st, x0, nHat, aInit, opt);
        kappa(i0,j0) = -res.curvature;
    }

    return kappa;
}


// ===================================================================
//  PC (PLIC-centroidal) curvature: parabola fit in local (t,n) frame
//  using segment centroids + weights:
//    w = Wendland(r,dWeight) * (segment_length * max(0, dot(n0, ni)))
// ===================================================================

struct PLICSegmentPC {
    int i, j;          // cell indices
    double xc, yc;     // centroid point of PLIC segment (here: midpoint)
    double nx, ny;     // unit normal (outward)
    double length;     // segment length
};
///----PV Implementation---------------//
struct PLICSegmentPV {
    int i, j;          // cell indices
    double xc, yc;     // segment midpoint
    double nx, ny;     // unit normal (outward)
    double length;     // segment length
    Point2D p1, p2;    // endpoints in global coordinates (REQUIRED by PV)
};
///----PV Implementation---------------//

// Wendland: (1+4q)(1-q)^4, q=r/d in [0,1]
static inline double wendland_PC(double r, double d)
{
    if (r <= 0.0) return 1.0;
    if (r >= d)   return 0.0;
    double q = 1.0 - r / d;
    double t = r / d;
    return (1.0 + 4.0 * t) * q * q * q * q;
}

// 3x3 solver (Gaussian elimination), returns false if singular
static inline bool solve3x3_PC(const double A_[3][3], const double b_[3], double x[3])
{
    double A[3][3];
    double b[3];
    for (int i=0;i<3;i++){
        b[i]=b_[i];
        for (int j=0;j<3;j++) A[i][j]=A_[i][j];
    }

    for (int k=0;k<3;k++){
        int piv=k;
        double best=std::fabs(A[k][k]);
        for (int i=k+1;i<3;i++){
            double v=std::fabs(A[i][k]);
            if(v>best){best=v;piv=i;}
        }
        if(best<1e-14) return false;
        if(piv!=k){
            for (int j=0;j<3;j++) std::swap(A[k][j],A[piv][j]);
            std::swap(b[k],b[piv]);
        }

        double diag=A[k][k];
        for(int j=k;j<3;j++) A[k][j]/=diag;
        b[k]/=diag;

        for(int i=0;i<3;i++){
            if(i==k) continue;
            double f=A[i][k];
            for(int j=k;j<3;j++) A[i][j]-=f*A[k][j];
            b[i]-=f*b[k];
        }
    }
    for(int i=0;i<3;i++) x[i]=b[i];
    return true;
}
///----PV Implementation---------------//
static std::vector<PLICSegmentPV> buildPLICSegmentsPV(const Grid2D& vof, const Config& cfg)
{
    std::vector<PLICSegmentPV> segs;
    auto iface = buildInterfaceCells(vof, cfg);
    if (iface.empty()) return segs;

    segs.reserve(iface.size());

    for (const auto& c : iface) {
        if (!(c.fg > cfg.eps_if && c.fg < 1.0 - cfg.eps_if)) continue;

        const double xL = cfg.prob_lo_x + c.i * vof.dx;
        const double xR = xL + vof.dx;
        const double yB = cfg.prob_lo_y + c.j * vof.dy;
        const double yT = yB + vof.dy;

        auto pts = intersectLineWithRect(c.nx, c.ny, c.d_line, xL, xR, yB, yT);
        if (pts.size() != 2) continue;

        const double dxs = pts[1].x - pts[0].x;
        const double dys = pts[1].y - pts[0].y;
        const double len = std::hypot(dxs, dys);
        if (len < 1e-14) continue;

        PLICSegmentPV s;
        s.i = c.i; s.j = c.j;
        s.p1 = pts[0];
        s.p2 = pts[1];
        s.xc = 0.5*(pts[0].x + pts[1].x);
        s.yc = 0.5*(pts[0].y + pts[1].y);

        // normal already unit from buildInterfaceCells, but normalize defensively
        double nlen = std::hypot(c.nx, c.ny);
        if (nlen < 1e-14) continue;
        s.nx = c.nx / nlen;
        s.ny = c.ny / nlen;

        s.length = len;
        segs.push_back(s);
    }

    return segs;
}
///----PV Implementation---------------//


// Build PLIC segments from your existing interface-cell reconstruction
static std::vector<PLICSegmentPC> buildPLICSegmentsPC(const Grid2D& vof, const Config& cfg)
{
    std::vector<PLICSegmentPC> segs;
    auto iface = buildInterfaceCells(vof, cfg);
    if (iface.empty()) return segs;

    segs.reserve(iface.size());

    for (const auto& c : iface) {
        if (!(c.fg > cfg.eps_if && c.fg < 1.0 - cfg.eps_if)) continue;

        double xL = cfg.prob_lo_x + c.i * vof.dx;
        double xR = xL + vof.dx;
        double yB = cfg.prob_lo_y + c.j * vof.dy;
        double yT = yB + vof.dy;

        auto pts = intersectLineWithRect(c.nx, c.ny, c.d_line, xL, xR, yB, yT);
        if (pts.size() != 2) continue;

        double dxs = pts[1].x - pts[0].x;
        double dys = pts[1].y - pts[0].y;
        double len = std::hypot(dxs, dys);
        if (len < 1e-14) continue;

        PLICSegmentPC s;
        s.i = c.i; s.j = c.j;
        s.xc = 0.5*(pts[0].x + pts[1].x);
        s.yc = 0.5*(pts[0].y + pts[1].y);

        // normal already unit from buildInterfaceCells (nx=-fx/|grad| etc)
        double nlen = std::hypot(c.nx, c.ny);
        if (nlen < 1e-14) continue;
        s.nx = c.nx / nlen;
        s.ny = c.ny / nlen;

        s.length = len;
        segs.push_back(s);
    }
    return segs;
}

///----PV Implementation---------------//

static bool curvature_PV_for_target(
    const std::vector<PLICSegmentPV>& segs,
    int ic, int jc,
    int stencilWidth,   // odd
    double dWeight,     // = cfg.pc_d_over_dx * dx   (reuse your PC knob)
    double& kappa_out,
    double& x0, double& y0,
    double& nx0, double& ny0)
{
    const PLICSegmentPV* tgt = nullptr;
    for (const auto& s : segs) {
        if (s.i == ic && s.j == jc) { tgt = &s; break; }
    }
    if (!tgt) return false;

    // target frame
    nx0 = tgt->nx; ny0 = tgt->ny;
    double nlen = std::hypot(nx0, ny0);
    if (nlen < 1e-14) return false;
    nx0 /= nlen; ny0 /= nlen;

    const double tx0 = -ny0;
    const double ty0 =  nx0;

    x0 = tgt->xc;
    y0 = tgt->yc;

    const int half = stencilWidth / 2;

    // Normal equations (PV)
    double A00=0, A01=0, A02=0;
    double A11=0, A12=0, A22=0;
    double B0 =0, B1 =0, B2 =0;

    int nUsed = 0;

    auto to_local = [&](double X, double Y, double& xp, double& zp) {
        const double rx = X - x0;
        const double ry = Y - y0;
        xp = rx*tx0 + ry*ty0;   // tangent
        zp = rx*nx0 + ry*ny0;   // normal
    };

    for (const auto& s : segs) {
        const int di = std::abs(s.i - ic);
        const int dj = std::abs(s.j - jc);
        if (di > half || dj > half) continue;

        double x1,z1,x2,z2;
        to_local(s.p1.x, s.p1.y, x1, z1);
        to_local(s.p2.x, s.p2.y, x2, z2);

        // Need a meaningful projected interval in x'
        if (std::fabs(x2 - x1) < 1e-14) continue;

        // Ensure ordered [xA,xB]
        double xA = x1, xB = x2;
        double zz1 = z1, zz2 = z2;
        if (xB < xA) { std::swap(xA, xB); std::swap(zz1, zz2); std::swap(x1, x2); std::swap(z1, z2); }

        const double Lproj = xB - xA;
        if (Lproj < 1e-14) continue;

        // Neighbor line g(x') through endpoints: z = b0 + b1 x
        const double b1 = (z2 - z1) / (x2 - x1);
        const double b0 = z1 - b1*x1;

        // Radial weight based on neighbor midpoint (in target frame)
        double xcm=0, zcm=0;
        to_local(s.xc, s.yc, xcm, zcm);

        const double r  = std::hypot(xcm, zcm);
        const double wR = wendland_PC(r, dWeight); // reuse your existing Wendland
        if (wR == 0.0) continue;

        // PV moments over [xA,xB]
        const double s0 = (xB - xA);
        const double s1 = 0.5 * (xB*xB - xA*xA);
        const double s2 = (1.0/3.0) * (xB*xB*xB - xA*xA*xA);

        // integral of neighbor line over [xA,xB]
        const double gInt = b0*s0 + b1*s1;

        // accumulate normal equations
        A00 += wR * s0*s0;
        A01 += wR * s0*s1;
        A02 += wR * s0*s2;
        A11 += wR * s1*s1;
        A12 += wR * s1*s2;
        A22 += wR * s2*s2;

        B0  += wR * s0 * gInt;
        B1  += wR * s1 * gInt;
        B2  += wR * s2 * gInt;

        ++nUsed;
    }

    if (nUsed < 3) return false;

    const double M[3][3] = {
        {A00, A01, A02},
        {A01, A11, A12},
        {A02, A12, A22}
    };
    const double bb[3] = {B0, B1, B2};
    double a[3];

    if (!solve3x3_PC(M, bb, a)) return false;

    const double fp  = a[1];
    const double fpp = 2.0 * a[2];

    const double denom = std::pow(1.0 + fp*fp, 1.5);
    if (denom < 1e-30) return false;

    kappa_out = -fpp / denom;
    return true;
}


///----PV Implementation---------------//

// Fit curvature at a target segment (ic,jc) using weighted parabola z = a0 + a1 x + a2 x^2
static bool curvature_PC_for_target(
    const std::vector<PLICSegmentPC>& segs,
    int ic, int jc,
    int stencilWidth,     // odd
    double dWeight,       // = cfg.pc_d_over_dx * dx
    double& kappa_out, double& x0, double& y0,
    double& nx0, double& ny0)
{
    const PLICSegmentPC* tgt = nullptr;
    for (const auto& s : segs) {
        if (s.i == ic && s.j == jc) { tgt = &s; break; }
    }
    if (!tgt) return false;

    // target frame
    nx0 = tgt->nx; ny0 = tgt->ny;
    double nlen = std::hypot(nx0, ny0);
    if (nlen < 1e-14) return false;
    nx0 /= nlen; ny0 /= nlen;

    const double tx0 = -ny0;
    const double ty0 =  nx0;

    const int half = stencilWidth / 2;
    x0 = tgt->xc;
    y0 = tgt->yc;

    // Normal equations for weighted LS on basis [1, x, x^2]
    double S00=0, S01=0, S02=0;
    double S11=0, S12=0, S22=0;
    double B0 =0, B1 =0, B2 =0;

    int nUsed = 0;

    for (const auto& s : segs) {
        int di = std::abs(s.i - ic);
        int dj = std::abs(s.j - jc);
        if (di > half || dj > half) continue;

        double rx = s.xc - x0;
        double ry = s.yc - y0;

        // local coords: x along tangent, z along normal
        double x = rx*tx0 + ry*ty0;
        double z = rx*nx0 + ry*ny0;

        double r = std::hypot(x, z);
        double wR = wendland_PC(r, dWeight);
        if (wR == 0.0) continue;

        // neighbor normal alignment + segment length weight
        double nxd = s.nx;
        double nyd = s.ny;
        double nld = std::hypot(nxd, nyd);     // FIXED BUG: was nyd*nyd mistyped in your snippet
        if (nld < 1e-14) continue;
        nxd /= nld; nyd /= nld;

        double dotnn = nx0*nxd + ny0*nyd;
        if (dotnn < 0.0) dotnn = 0.0;

        double wA = s.length * dotnn;
        if (wA <= 0.0) continue;

        double w = wR * wA;

        const double x2 = x*x;
        const double x3 = x2*x;
        const double x4 = x2*x2;

        S00 += w;
        S01 += w * x;
        S02 += w * x2;
        S11 += w * x2;
        S12 += w * x3;
        S22 += w * x4;

        B0  += w * z;
        B1  += w * x * z;
        B2  += w * x2 * z;

        ++nUsed;
    }

    if (nUsed < 3) return false;

    double A[3][3] = {
        { S00, S01, S02 },
        { S01, S11, S12 },
        { S02, S12, S22 }
    };
    double b[3] = { B0, B1, B2 };
    double a[3];

    if (!solve3x3_PC(A, b, a)) return false;

    // curvature at x=0: kappa_local = -f''/(1+f'^2)^(3/2), f'(0)=a1, f''(0)=2a2
    const double fp  = a[1];
    const double fpp = 2.0 * a[2];
    const double denom = std::pow(1.0 + fp*fp, 1.5);
    if (denom < 1e-30) return false;

    kappa_out = -fpp / denom;
    return true;
}

///----PV Implementation---------------//

static Grid2D computeCurvaturePV(const Grid2D& vof, const Config& cfg)
{
    Grid2D kappa(vof.nx, vof.ny, vof.dx, vof.dy);

    auto segs = buildPLICSegmentsPV(vof, cfg);
    if (segs.empty()) return kappa;

    const int S = std::max(3, (cfg.pc_stencil % 2 == 0 ? cfg.pc_stencil+1 : cfg.pc_stencil));
    const double dWeight = cfg.pc_d_over_dx * vof.dx;

    for (const auto& s : segs) {
        double k=0, x0=0, y0=0, nx0=0, ny0=0;
        if (!curvature_PV_for_target(segs, s.i, s.j, S, dWeight, k, x0, y0, nx0, ny0))
            continue;

        if (cfg.ic_shape == "circle") {
            InterfaceCell tmp{};
            tmp.xc = x0; tmp.yc = y0;
            tmp.nx = nx0; tmp.ny = ny0;
            const double sgn_global = global_sign_from_circle_outward(tmp, cfg);
            kappa(s.i, s.j) = sgn_global * std::fabs(k);
        } else {
            kappa(s.i, s.j) = k;
        }
    }

    return kappa;
}


///----PV Implementation---------------//

static Grid2D computeCurvaturePC(const Grid2D& vof, const Config& cfg)
{
    Grid2D kappa(vof.nx, vof.ny, vof.dx, vof.dy);

    auto segs = buildPLICSegmentsPC(vof, cfg);
    if (segs.empty()) return kappa;

    const int S = std::max(3, (cfg.pc_stencil % 2 == 0 ? cfg.pc_stencil+1 : cfg.pc_stencil));
    const double dWeight = cfg.pc_d_over_dx * vof.dx;

    for (const auto& s : segs) {
        double k=0, x0=0, y0=0, nx0=0, ny0=0;
        if (!curvature_PC_for_target(segs, s.i, s.j, S, dWeight, k, x0, y0, nx0, ny0))
            continue;

        if (cfg.ic_shape == "circle") {
            InterfaceCell tmp{};
            tmp.xc = x0; tmp.yc = y0;
            tmp.nx = nx0; tmp.ny = ny0;
            const double sgn_global = global_sign_from_circle_outward(tmp, cfg);
            kappa(s.i, s.j) = sgn_global * std::fabs(k);
        } else {
            kappa(s.i, s.j) = k;
        }
    }

    return kappa;
}

// ====================== Error utilities (CONSISTENT) ===============
struct ErrNorms {
    double L2   = 0.0;
    double Linf = 0.0;
    int    N    = 0;
};

static ErrNorms ErrorParaboloidConsistentCircle(const Grid2D& kappa_cell,
                                                const Grid2D& vof,
                                                const Config& cfg)
{
    ErrNorms E;
    const double k_exact = 1.0 / cfg.R;
    const double eps_if  = cfg.eps_if;

    auto iface = buildInterfaceCells(vof, cfg);
    for (const auto& c : iface) {
        if (!(c.fg > eps_if && c.fg < 1.0 - eps_if)) continue;
        double kI = kappa_cell(c.i, c.j);
        double e  = kI - k_exact;
        E.Linf = std::max(E.Linf, std::fabs(e));
        E.L2   += e * e;
        E.N    += 1;
    }
    if (E.N > 0) E.L2 = std::sqrt(E.L2 / E.N);
    return E;
}

static ErrNorms ErrorParaboloidConsistentCosine(const Grid2D& kappa_cell,
                                                const Grid2D& vof,
                                                const Config& cfg)
{
    ErrNorms E;
    const double eps_if  = cfg.eps_if;
    const double alpha = 2.0 * M_PI * cfg.cos_n / cfg.cos_L;

    auto iface = buildInterfaceCells(vof, cfg);
    for (const auto& c : iface) {
        if (!(c.fg > eps_if && c.fg < 1.0 - eps_if)) continue;

        double xI = c.xint;
        double phi = alpha * xI;
        double yx  = alpha * std::sin(phi);
        double yxx = alpha * alpha * std::cos(phi);
        double k_exact = - yxx / std::pow(1.0 + yx*yx, 1.5);

        double kI = kappa_cell(c.i, c.j);
        double e  = kI - k_exact;

        E.Linf = std::max(E.Linf, std::fabs(e));
        E.L2   += e * e;
        E.N    += 1;
    }
    if (E.N > 0) E.L2 = std::sqrt(E.L2 / E.N);
    return E;
}

static void write_ppm(const std::string &filename, int W, int H,
                      const std::vector<unsigned char> &rgb) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        std::cerr << "ERROR: Could not open " << filename << " for writing!\n";
        return;
    }
    ofs << "P6\n" << W << " " << H << "\n255\n";
    ofs.write(reinterpret_cast<const char*>(rgb.data()), (std::streamsize)rgb.size());
    ofs.close();
}

static void dump_alpha_and_curvature_ppm(const Grid2D& alpha,
                                        const Grid2D& kappa,
                                        const Config& cfg)
{
    const int nx = alpha.nx, ny = alpha.ny;
    std::vector<unsigned char> alpha_rgb(nx * ny * 3);
    std::vector<unsigned char> curv_rgb(nx * ny * 3);

    // compute min/max over interface cells
    double kmin = 1e300, kmax = -1e300;
    for (int j=0;j<ny;++j) for(int i=0;i<nx;++i){
        double f = alpha(i,j);
        if (f > cfg.eps_if && f < 1.0 - cfg.eps_if) {
            double v = kappa(i,j);
            if (std::isfinite(v)) { kmin = std::min(kmin,v); kmax = std::max(kmax,v); }
        }
    }
    if (!(kmin < kmax)) { kmin = -1.0; kmax = 1.0; }
    double maxabs = std::max(std::fabs(kmin), std::fabs(kmax));
    if (maxabs < 1e-12) maxabs = 1e-12;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int id = (j * nx + i) * 3;

            unsigned char g = (unsigned char)(255.0 * clamp01(alpha(i,j)));
            alpha_rgb[id+0] = g;
            alpha_rgb[id+1] = g;
            alpha_rgb[id+2] = g;

            double f = alpha(i,j);
            bool iface = (f > cfg.eps_if && f < 1.0 - cfg.eps_if);

            if (iface) {
                double v = kappa(i,j);
                double t = (v + maxabs) / (2.0 * maxabs);
                t = clamp01(t);

                unsigned char r, gg, b;
                if (t < 0.5) {
                    double u = t / 0.5;
                    r  = (unsigned char)(255 * u);
                    gg = (unsigned char)(255 * u);
                    b  = 255;
                } else {
                    double u = (t - 0.5) / 0.5;
                    r  = 255;
                    gg = (unsigned char)(255 * (1 - u));
                    b  = (unsigned char)(255 * (1 - u));
                }
                curv_rgb[id+0] = r;
                curv_rgb[id+1] = gg;
                curv_rgb[id+2] = b;
            } else {
                unsigned char bg = 120;
                curv_rgb[id+0] = bg;
                curv_rgb[id+1] = bg;
                curv_rgb[id+2] = bg;
            }
        }
    }

    write_ppm("alpha.ppm", nx, ny, alpha_rgb);
    write_ppm("curvature.ppm", nx, ny, curv_rgb);
    std::cout << "alpha.ppm / curvature.ppm written (nx=" << nx << ", ny=" << ny << ")\n";
}

// ====================== main() =====================================
int main(int argc, char** argv) {
    std::string input_file = "curv.inp";
    if (argc > 1) input_file = argv[1];

    Config cfg;
    load_config(input_file, cfg);

    std::vector<int> res_list = cfg.n_list.empty() ? std::vector<int>{cfg.nx} : cfg.n_list;

    std::cout << "=== Curvature convergence test ===\n";
    std::cout << "  ic.shape      = " << cfg.ic_shape << "\n";
    std::cout << "  domain x: [" << cfg.prob_lo_x << ", " << cfg.prob_hi_x << "]\n";
    std::cout << "  domain y: [" << cfg.prob_lo_y << ", " << cfg.prob_hi_y << "]\n";

    if (cfg.ic_shape == "circle") {
        std::cout << "  circle R      = " << cfg.R
                << " , center = (" << cfg.cx << ", " << cfg.cy << ")\n";
    } else if (cfg.ic_shape == "cosine") {
        std::cout << "  cosine A      = " << cfg.cos_A
                << " , L = " << cfg.cos_L
                << " , n = " << cfg.cos_n << "\n";
    }

    std::cout << "  vof.init_method = " << cfg.vof_init_method << "\n";
    std::cout << "  curv.method   = " << cfg.curv_method << "\n";
    std::cout << "  h_factor      = " << cfg.h_factor << " (h = h_factor * dx)\n";

    std::cout << "  rdf.grad_stencil = " << cfg.rdf_grad_stencil << "\n";
    std::cout << "  rdf.A         = " << cfg.rdf_A << "\n";
    std::cout << "  rdf.B         = " << cfg.rdf_B << "\n";
    std::cout << "  rdf.rad_cells = " << cfg.rdf_rad_cells << "\n";

    std::cout << "  hf.m          = " << cfg.hf_m << "\n";
    std::cout << "  error.eps_if  = " << cfg.eps_if << "\n";
    std::cout << "  dump.ppm      = " << cfg.dump_ppm << "\n";

    if (cfg.curv_method == "VFM" || cfg.curv_method == "vfm") {
        std::cout << "  vfm.stencil              = " << cfg.vfm_stencil << "\n";
        std::cout << "  vfm.weight_radius_factor = " << cfg.vfm_weight_radius_factor << "\n";
        std::cout << "  vfm.maxIters             = " << cfg.vfm_maxIters << "\n";
        std::cout << "  vfm.tolRel               = " << cfg.vfm_tolRel << "\n";
        std::cout << "  vfm.lambda0              = " << cfg.vfm_lambda0 << "\n";
        std::cout << "  vfm.lambdaUp             = " << cfg.vfm_lambdaUp << "\n";
        std::cout << "  vfm.lambdaDown           = " << cfg.vfm_lambdaDown << "\n";
        std::cout << "  vfm.fdEps                = " << cfg.vfm_fdEps << "\n";
    }
    if (cfg.curv_method == "PC") {
        std::cout << "  pc.stencil      = " << cfg.pc_stencil << "\n";
        std::cout << "  pc.d_over_dx    = " << cfg.pc_d_over_dx << "\n";
    }
    ///----PV Implementation---------------//
    if (cfg.curv_method == "PV") {
    std::cout << "  pc.stencil      = " << cfg.pc_stencil << "\n";
    std::cout << "  pc.d_over_dx    = " << cfg.pc_d_over_dx << "\n";
    }
    ///----PV Implementation---------------//


    std::cout << "==================================\n\n";


    std::cout << std::scientific << std::setprecision(6);

    std::vector<double> dx_list;
    std::vector<double> err_list;

    for (size_t idxRes = 0; idxRes < res_list.size(); ++idxRes) {
        int N  = res_list[idxRes];
        int nx = N;
        int ny = N;

        double Lx = cfg.prob_hi_x - cfg.prob_lo_x;

        // ---- VOF init ----
        Grid2D vof;
        if (cfg.ic_shape == "circle") {
            if (cfg.vof_init_method == "exact_poly") vof = initializeCircleVOF_exact_poly(nx, ny, cfg);
            else                                     vof = initializeCircleVOF_recursive(nx, ny, cfg);
        } else if (cfg.ic_shape == "cosine") {
            // exact_poly only implemented for circle; cosine uses original recursive init
            vof = initializeCosineVOF_recursive(nx, ny, cfg);
        } else {
            std::cerr << "Unknown ic.shape = " << cfg.ic_shape << ", fallback to circle.\n";
            vof = (cfg.vof_init_method == "exact_poly")
                    ? initializeCircleVOF_exact_poly(nx, ny, cfg)
                    : initializeCircleVOF_recursive(nx, ny, cfg);
        }

        double dx_cell = Lx / nx;
        double h_conv  = cfg.h_factor * dx_cell;

        Grid2D kappa;
        Grid2D phi_for_dump;
        bool dump_phi = false;

        if (cfg.curv_method == "CV") {
            Grid2D f_tilde = convolveVOF(vof, h_conv);
            kappa = computeCurvatureCV(f_tilde);
        } else if (cfg.curv_method == "RDF_raw") {
            Grid2D phi = buildRDF(vof, cfg);
            kappa      = computeCurvatureRDF(phi);
            phi_for_dump = phi;
            dump_phi = true;
        } else if (cfg.curv_method == "RDF_smooth") {
            Grid2D phi       = buildRDF(vof, cfg);
            Grid2D phi_smoid = convolveField(phi, h_conv);
            kappa            = computeCurvatureRDF(phi_smoid);
            phi_for_dump     = phi;
            dump_phi = true;
        } else if (cfg.curv_method == "PARABOLOID") {
            kappa = computeCurvatureParabolicReconstructionPaperCartesian(vof, cfg);
        } else if (cfg.curv_method == "HF") {
            kappa = computeCurvatureHF(vof, cfg);
        } else if (cfg.curv_method == "VFM" || cfg.curv_method == "vfm") {
            kappa = computeCurvatureVFM(vof, cfg);
        } else if (cfg.curv_method == "PC") {
            kappa = computeCurvaturePC(vof, cfg);

        ///----PV Implementation---------------//
        }else if (cfg.curv_method == "PV") {
            kappa = computeCurvaturePV(vof, cfg); 
        ///----PV Implementation---------------//       
        } else {
            std::cerr << "Unknown curv.method = " << cfg.curv_method << ", fallback to CV.\n";
            Grid2D f_tilde = convolveVOF(vof, h_conv);
            kappa = computeCurvatureCV(f_tilde);
        }

        ErrNorms E;
        if (cfg.ic_shape == "circle") E = ErrorParaboloidConsistentCircle(kappa, vof, cfg);
        else                          E = ErrorParaboloidConsistentCosine(kappa, vof, cfg);

        double L2_err = E.L2;

        std::cout << "N = " << std::setw(6) << N
                  << "  dx = " << std::setw(12) << dx_cell
                  << "  L2_error = " << std::setw(12) << L2_err
                  << "  Linf = " << std::setw(12) << E.Linf
                  << "  Niface = " << E.N
                  << "\n";

        dx_list.push_back(dx_cell);
        err_list.push_back(L2_err);

        if (idxRes == res_list.size() - 1) {
            // dump kappa field
            {
                std::ofstream fout("kappa_field.dat");
                if (fout) {
                    fout << "# i j x y kappa\n";
                    fout << std::scientific << std::setprecision(8);
                    for (int j = 0; j < ny; ++j) {
                        double y = cfg.prob_lo_y + (j + 0.5) * kappa.dy;
                        for (int i = 0; i < nx; ++i) {
                            double x = cfg.prob_lo_x + (i + 0.5) * kappa.dx;
                            fout << i << "  " << j << "  "
                                 << x << "  " << y << "  "
                                 << kappa(i,j) << "\n";
                        }
                    }
                    std::cout << "kappa_field.dat written for N = " << N << "\n";
                }
            }

            // dump phi if requested
            if (dump_phi) {
                std::ofstream fphi("phi_field.dat");
                if (fphi) {
                    fphi << "# i j x y phi\n";
                    fphi << std::scientific << std::setprecision(8);
                    for (int j = 0; j < ny; ++j) {
                        double y = cfg.prob_lo_y + (j + 0.5) * phi_for_dump.dy;
                        for (int i = 0; i < nx; ++i) {
                            double x = cfg.prob_lo_x + (i + 0.5) * phi_for_dump.dx;
                            fphi << i << "  " << j << "  "
                                 << x << "  " << y << "  "
                                 << phi_for_dump(i,j) << "\n";
                        }
                    }
                    std::cout << "phi_field.dat written for N = " << N << "\n";
                }
            }

            // dump interface segments (your existing feature)
            {
                auto iface = buildInterfaceCells(vof, cfg);
                std::ofstream fint("interface_segments.dat");
                if (fint) {
                    fint << "# i j x1 y1 x2 y2 xint yint fg nx ny d\n";
                    fint << std::scientific << std::setprecision(10);

                    for (const auto& c : iface) {
                        if (!(c.fg > cfg.eps_if && c.fg < 1.0 - cfg.eps_if)) continue;

                        double xL = cfg.prob_lo_x + c.i * vof.dx;
                        double xR = xL + vof.dx;
                        double yB = cfg.prob_lo_y + c.j * vof.dy;
                        double yT = yB + vof.dy;

                        auto pts = intersectLineWithRect(c.nx, c.ny, c.d_line, xL, xR, yB, yT);
                        if (pts.size() != 2) continue;

                        fint << c.i << " " << c.j << " "
                             << pts[0].x << " " << pts[0].y << " "
                             << pts[1].x << " " << pts[1].y << " "
                             << c.xint << " " << c.yint << " "
                             << c.fg << " "
                             << c.nx << " " << c.ny << " "
                             << c.d_line << "\n";
                    }
                    std::cout << "interface_segments.dat written for N = " << N << "\n";
                }
            }

            // ---- PPM dump (from your demo) ----
            if (cfg.dump_ppm) {
                dump_alpha_and_curvature_ppm(vof, kappa, cfg);
            }
        }
    }

    if (!dx_list.empty()) {
        int ref_idx = (int)dx_list.size() - 1;
        double dx_ref  = dx_list[ref_idx];
        double err_ref = err_list[ref_idx];
        double C = (dx_ref > 0.0) ? err_ref / (dx_ref * dx_ref) : 0.0;

        std::cout << "\n# Reference second-order line (anchored at finest grid)\n";
        std::cout << "#   dx           L2_error       L2_second_order\n";
        for (size_t k = 0; k < dx_list.size(); ++k) {
            double dx  = dx_list[k];
            double err = err_list[k];
            double err2 = C * dx * dx;
            std::cout << std::setw(12) << dx << "  "
                      << std::setw(12) << err << "  "
                      << std::setw(12) << err2 << "\n";
        }
    }

    return 0;
}
