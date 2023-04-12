#include <AMReX.H>
#include <AMReX_Parser.H>
#include <AMReX_IParser.H>
#include <map>

using namespace amrex;

static int max_stack_size = 0;
static int test_number = 0;

template <typename F>
int test1 (std::string const& f,
           std::map<std::string,Real> const& constants,
           Vector<std::string> const& variables,
           F && fb, Array<Real,1> const& lo, Array<Real,1> const& hi,
           int N, Real reltol, Real abstol)
{
    amrex::Print() << test_number++ << ". Testing \"" << f << "\"   ";

    Parser parser(f);
    for (auto const& kv : constants) {
        parser.setConstant(kv.first, kv.second);
    }
    parser.registerVariables(variables);
    auto const exe = parser.compile<1>();
    max_stack_size = std::max(max_stack_size, parser.maxStackSize());

    GpuArray<Real,1> dx{(hi[0]-lo[0]) / (N-1)};

    int nfail = 0;
    Real max_relerror = 0.;
    for (int i = 0; i < N; ++i) {
        Real x = lo[0] + i*dx[0];
        Real result = exe(x);
        Real benchmark = fb(x);
        Real abserror = std::abs(result-benchmark);
        Real relerror = abserror / (1.e-50 + std::max(std::abs(result),std::abs(benchmark)));
        if (abserror > abstol && relerror > reltol) {
            amrex::Print() << "\n    f(" << x << ") = " << result << ", "
                           << benchmark;
            max_relerror = std::max(max_relerror, relerror);
            ++nfail;
        }
    }
    if (nfail > 0) {
        amrex::Print() << "\n    failed " << nfail << " times.  Max rel. error: "
                       << max_relerror << "\n";
        return 1;
    } else {
        amrex::Print() << "    pass\n";
        return 0;
    }
}

template <typename F>
int test3 (std::string const& f,
           std::map<std::string,Real> const& constants,
           Vector<std::string> const& variables,
           F && fb, Array<Real,3> const& lo, Array<Real,3> const& hi,
           int N, Real reltol, Real abstol)
{
    amrex::Print() << test_number++ << ". Testing \"" << f << "\"   ";

    Parser parser(f);
    for (auto const& kv : constants) {
        parser.setConstant(kv.first, kv.second);
    }
    parser.registerVariables(variables);
    auto const exe = parser.compile<3>();
    max_stack_size = std::max(max_stack_size, parser.maxStackSize());

    GpuArray<Real,3> dx{(hi[0]-lo[0]) / (N-1),
                        (hi[1]-lo[1]) / (N-1),
                        (hi[2]-lo[2]) / (N-1)};
    int nfail = 0;
    for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
    for (int k = 0; k < N; ++k) {
        Real x = lo[0] + i*dx[0];
        Real y = lo[1] + j*dx[1];
        Real z = lo[2] + k*dx[2];
        Real result = exe(x,y,z);
        Real benchmark = fb(x,y,z);
        Real abserror = std::abs(result-benchmark);
        Real relerror = abserror / (1.e-50 + std::max(std::abs(result),std::abs(benchmark)));
        if (abserror > abstol && relerror > reltol) {
            amrex::Print() << "    f(" << x << "," << y << "," << z << ") = " << result << ", "
                           << benchmark << "\n";
            ++nfail;
        }
    }}}
    if (nfail > 0) {
        amrex::Print() << "    failed " << nfail << " times\n";
        return 1;
    } else {
        amrex::Print() << "    pass\n";
        return 0;
    }
}

template <typename F>
int test4 (std::string const& f,
           std::map<std::string,Real> const& constants,
           Vector<std::string> const& variables,
           F && fb, Array<Real,4> const& lo, Array<Real,4> const& hi,
           int N, Real reltol, Real abstol)
{
    amrex::Print() << test_number++ << ". Testing \"" << f << "\"   ";

    Parser parser(f);
    for (auto const& kv : constants) {
        parser.setConstant(kv.first, kv.second);
    }
    parser.registerVariables(variables);
    auto const exe = parser.compile<4>();
    max_stack_size = std::max(max_stack_size, parser.maxStackSize());

    GpuArray<Real,4> dx{(hi[0]-lo[0]) / (N-1),
                        (hi[1]-lo[1]) / (N-1),
                        (hi[2]-lo[2]) / (N-1),
                        (hi[3]-lo[3]) / (N-1)};
    int nfail = 0;
    for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
    for (int k = 0; k < N; ++k) {
    for (int m = 0; m < N; ++m) {
        Real x = lo[0] + i*dx[0];
        Real y = lo[1] + j*dx[1];
        Real z = lo[2] + k*dx[2];
        Real t = lo[3] + m*dx[3];
        Real result = exe(x,y,z,t);
        Real benchmark = fb(x,y,z,t);
        Real abserror = std::abs(result-benchmark);
        Real relerror = abserror / (1.e-50 + std::max(std::abs(result),std::abs(benchmark)));
        if (abserror > abstol && relerror > reltol) {
            amrex::Print() << "    f(" << x << "," << y << "," << z << "," << t << ") = " << result << ", "
                           << benchmark << "\n";
            ++nfail;
        }
    }}}}
    if (nfail > 0) {
        amrex::Print() << "    failed " << nfail << " times\n";
        return 1;
    } else {
        amrex::Print() << "    pass\n";
        return 0;
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        amrex::Print() << "\n";
        int nerror = 0;
        nerror += test3("if( ((z-zc)*(z-zc)+(y-yc)*(y-yc)+(x-xc)*(x-xc))^(0.5) < (r_star-dR), 0.0, if(((z-zc)*(z-zc)+(y-yc)*(y-yc)+(x-xc)*(x-xc))^(0.5) <= r_star, dens, 0.0))",
                        {{"xc", 0.1}, {"yc", -1.0}, {"zc", 0.2}, {"r_star", 0.73}, {"dR", 0.57}, {"dens", 12.}},
                        {"x","y","z"},
                        [=] (Real x, Real y, Real z) -> Real {
                            Real xc=0.1, yc=-1.0, zc=0.2, r_star=0.73, dR=0.57, dens=12.;
                            Real r = std::sqrt((z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc));
                            if (r >= r_star-dR && r <= r_star) {
                                return dens;
                            } else {
                                return 0.0;
                            }
                        },
                        {-1., -1., -1.0}, {1.0, 1.0, 1.0}, 100,
                        1.e-12, 1.e-15);

        nerror += test3("r=sqrt((z-zc)*(z-zc)+(y-yc)*(y-yc)+(x-xc)*(x-xc)); if(r < (r_star-dR), 0.0, if(r <= r_star, dens, 0.0))",
                        {{"xc", 0.1}, {"yc", -1.0}, {"zc", 0.2}, {"r_star", 0.73}, {"dR", 0.57}, {"dens", 12.}},
                        {"x","y","z"},
                        [=] (Real x, Real y, Real z) -> Real {
                            Real xc=0.1, yc=-1.0, zc=0.2, r_star=0.73, dR=0.57, dens=12.;
                            Real r = std::sqrt((z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc));
                            if (r >= r_star-dR && r <= r_star) {
                                return dens;
                            } else {
                                return 0.0;
                            }
                        },
                        {-1., -1., -1.0}, {1.0, 1.0, 1.0}, 100,
                        1.e-12, 1.e-15);

        nerror += test3("r2=(z-zc)*(z-zc)+(y-yc)*(y-yc)+(x-xc)*(x-xc); r=sqrt(r2); if(r < (r_star-dR), 0.0, if(r <= r_star, dens, 0.0))",
                        {{"xc", 0.1}, {"yc", -1.0}, {"zc", 0.2}, {"r_star", 0.73}, {"dR", 0.57}, {"dens", 12.}},
                        {"x","y","z"},
                        [=] (Real x, Real y, Real z) -> Real {
                            Real xc=0.1, yc=-1.0, zc=0.2, r_star=0.73, dR=0.57, dens=12.;
                            Real r = std::sqrt((z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc));
                            if (r >= r_star-dR && r <= r_star) {
                                return dens;
                            } else {
                                return 0.0;
                            }
                        },
                        {-1., -1., -1.0}, {1.0, 1.0, 1.0}, 100,
                        1.e-12, 1.e-15);

        nerror += test3("( ((( (z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc) )^(0.5))<=r_star) * ((( (z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc) )^(0.5))>=(r_star-dR)) )*dens",
                        {{"xc", 0.1}, {"yc", -1.0}, {"zc", 0.2}, {"r_star", 0.73}, {"dR", 0.57}, {"dens", 12.}},
                        {"x","y","z"},
                        [=] (Real x, Real y, Real z) -> Real {
                            Real xc=0.1, yc=-1.0, zc=0.2, r_star=0.73, dR=0.57, dens=12.;
                            Real r = std::sqrt((z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc));
                            if (r >= r_star-dR && r <= r_star) {
                                return dens;
                            } else {
                                return 0.0;
                            }
                        },
                        {-1., -1., -1.0}, {1.0, 1.0, 1.0}, 100,
                        1.e-12, 1.e-15);

        nerror += test4("( (( (z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc) )^(0.5))>=r_star)*(-( (t<to)*(t/to)*omega + (t>=to)*omega )*(((x-xc)*(x-xc) + (y-yc)*(y-yc))^(0.5))/((1.0-( ( (t<to)*(t/to)*omega + (t>=to)*omega)  *(((x-xc)*(x-xc) + (y-yc)*(y-yc))^(0.5))/c)^2)^(0.5)) * (y-yc)/(((x-xc)*(x-xc) + (y-yc)*(y-yc))^(0.5)))",
                        {{"xc", 0.1}, {"yc", -1.0}, {"zc", 0.2}, {"to", 3.}, {"omega", 0.33}, {"c", 30.}, {"r_star", 0.75}},
                        {"x","y","z","t"},
                        [=] (Real x, Real y, Real z, Real t) -> Real {
                            Real xc=0.1, yc=-1.0, zc=0.2, to=3., omega=0.33, c=30., r_star=0.75;
                            Real r = std::sqrt((z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc));
                            if (r >= r_star) {
                                Real tomega = (t>=to) ? omega : omega*(t/to);
                                Real r2d = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
                                return -tomega * r / std::sqrt(1.0-(tomega*r2d/c)*(tomega*r2d/c)) * (y-yc)/r;
                            } else {
                                return 0.0;
                            }
                        },
                        {-1., -1., -1.0, 0.0}, {1.0, 1.0, 1.0, 10}, 30,
                        1.e-12, 1.e-15);

        nerror += test4("r=sqrt((z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc)); tomega=if(t>=to, omega, omega*(t/to)); r2d=sqrt((y-yc)*(y-yc) + (x-xc)*(x-xc)); (r>=r_star)*(-tomega*r/(1.0-((tomega*r2d/c)^2))^0.5 * (y-yc)/r)",
                        {{"xc", 0.1}, {"yc", -1.0}, {"zc", 0.2}, {"to", 3.}, {"omega", 0.33}, {"c", 30.}, {"r_star", 0.75}},
                        {"x","y","z","t"},
                        [=] (Real x, Real y, Real z, Real t) -> Real {
                            Real xc=0.1, yc=-1.0, zc=0.2, to=3., omega=0.33, c=30., r_star=0.75;
                            Real r = std::sqrt((z-zc)*(z-zc) + (y-yc)*(y-yc) + (x-xc)*(x-xc));
                            if (r >= r_star) {
                                Real tomega = (t>=to) ? omega : omega*(t/to);
                                Real r2d = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
                                return -tomega * r / std::sqrt(1.0-(tomega*r2d/c)*(tomega*r2d/c)) * (y-yc)/r;
                            } else {
                                return 0.0;
                            }
                        },
                        {-1., -1., -1.0, 0.0}, {1.0, 1.0, 1.0, 10}, 30,
                        1.e-12, 1.e-15);

        nerror += test3("cos(m * pi / Lx * (x - Lx / 2)) * cos(n * pi / Ly * (y - Ly / 2)) * sin(p * pi / Lz * (z - Lz / 2))*mu_0*(x>-0.5)*(x<0.5)*(y>-0.5)*(y<0.5)*(z>-0.5)*(z<0.5)",
                        {{"m", 0.0}, {"n", 1.0}, {"pi", 3.14}, {"p", 1.0}, {"Lx", 1.}, {"Ly", 1.}, {"Lz", 1.}, {"mu_0", 1.27e-6}},
                        {"x","y","z"},
                        [=] (Real x, Real y, Real z) -> Real {
                            Real m=0.0,n=1.0,pi=3.14,p=1.0,Lx=1.,Ly=1.,Lz=1.,mu_0=1.27e-6;
                            if ((x>-0.5) && (x<0.5) && (y>-0.5) && (y<0.5) && (z>-0.5) &&(z<0.5)) {
                                return std::cos(m * pi / Lx * (x - Lx / 2)) * std::cos(n * pi / Ly * (y - Ly / 2)) * std::sin(p * pi / Lz * (z - Lz / 2))*mu_0;
                            } else {
                                return 0.0;
                            }
                        },
                        {-0.8, -0.8, -0.8}, {0.8, 0.8, 0.8}, 100,
                        1.e-12, 1.e-15);

        nerror += test3("if ((x>-0.5) and (x<0.5) and (y>-0.5) and (y<0.5) and (z>-0.5) and (z<0.5), cos(m * pi / Lx * (x - Lx / 2)) * cos(n * pi / Ly * (y - Ly / 2)) * sin(p * pi / Lz * (z - Lz / 2))*mu_0*(x>-0.5)*(x<0.5)*(y>-0.5)*(y<0.5)*(z>-0.5)*(z<0.5), 0)",
                        {{"m", 0.0}, {"n", 1.0}, {"pi", 3.14}, {"p", 1.0}, {"Lx", 1.}, {"Ly", 1.}, {"Lz", 1.}, {"mu_0", 1.27e-6}},
                        {"x","y","z"},
                        [=] (Real x, Real y, Real z) -> Real {
                            Real m=0.0,n=1.0,pi=3.14,p=1.0,Lx=1.,Ly=1.,Lz=1.,mu_0=1.27e-6;
                            if ((x>-0.5) && (x<0.5) && (y>-0.5) && (y<0.5) && (z>-0.5) &&(z<0.5)) {
                                return std::cos(m * pi / Lx * (x - Lx / 2)) * std::cos(n * pi / Ly * (y - Ly / 2)) * std::sin(p * pi / Lz * (z - Lz / 2))*mu_0;
                            } else {
                                return 0.0;
                            }
                        },
                        {-0.8, -0.8, -0.8}, {0.8, 0.8, 0.8}, 100,
                        1.e-12, 1.e-15);

        nerror += test3("2.*sqrt(2.)+sqrt(-log(x))*cos(2*pi*z)",
                        {{"pi", 3.14}},
                        {"x","y","z"},
                        [=] (Real x, Real, Real z) -> Real {
                            Real pi = 3.14;
                            return 2.*std::sqrt(2.)+std::sqrt(-std::log(x))*std::cos(2*pi*z);
                        },
                        {0.5, 0.8, 0.3}, {16, 16, 16}, 100,
                        1.e-12, 1.e-15);

        nerror += test1("nc*n0*(if(abs(z)<=r0, 1.0, if(abs(z)<r0+Lcut, exp((-abs(z)+r0)/L), 0.0)))",
                        {{"nc",1.742e27},{"n0",30.},{"r0",2.5e-6},{"Lcut",2.e-6},{"L",0.05e-6}},
                        {"z"},
                        [=] (Real z) -> Real {
                            Real nc=1.742e27, n0=30., r0=2.5e-6, Lcut=2.e-6, L=0.05e-6;
                            if (std::abs(z) <= r0) {
                                return nc*n0;
                            } else if (std::abs(z) < r0+Lcut) {
                                return nc*n0*std::exp((-std::abs(z)+r0)/L);
                            } else {
                                return 0.0;
                            }
                        },
                        {-5.e-6}, {25.e-6}, 10000,
                        1.e-12, 1.e-15);

        nerror += test1("(z<lramp)*0.5*(1-cos(pi*z/lramp))*dens+(z>lramp)*dens",
                        {{"lramp",8.e-3},{"pi",3.14},{"dens",1.e23}},
                        {"z"},
                        [=] (Real z) -> Real {
                            Real lramp=8.e-3, pi=3.14, dens=1.e23;
                            if (z < lramp) {
                                return 0.5*(1-std::cos(pi*z/lramp))*dens;
                            } else {
                                return dens;
                            }
                        },
                        {-149.e-6}, {1.e-6}, 1000,
                        1.e-12, 1.e-15);

        nerror += test1("if(z<lramp, 0.5*(1-cos(pi*z/lramp))*dens, dens)",
                        {{"lramp",8.e-3},{"pi",3.14},{"dens",1.e23}},
                        {"z"},
                        [=] (Real z) -> Real {
                            Real lramp=8.e-3, pi=3.14, dens=1.e23;
                            if (z < lramp) {
                                return 0.5*(1-std::cos(pi*z/lramp))*dens;
                            } else {
                                return dens;
                            }
                        },
                        {-149.e-6}, {1.e-6}, 1000,
                        1.e-12, 1.e-15);

        nerror += test1("if(z<zp, nc*exp((z-zc)/lgrad), if(z<=zp2, 2.*nc, nc*exp(-(z-zc2)/lgrad)))",
                        {{"zc",20.e-6},{"zp",20.05545177444479562e-6},{"nc",1.74e27},{"lgrad",0.08e-6},{"zp2",24.e-6},{"zc2",24.05545177444479562e6}},
                        {"z"},
                        [=] (Real z) -> Real {
                            Real zc=20.e-6, zp=20.05545177444479562e-6, nc=1.74e27, lgrad=0.08e-6, zp2=24.e-6, zc2=24.05545177444479562e6;
                            if (z < zp) {
                                return nc*std::exp((z-zc)/lgrad);
                            } else if (z <= zp2) {
                                return 2.*nc;
                            } else {
                                return nc*exp(-(z-zc2)/lgrad);
                            }
                        },
                        {0.}, {100.e-6}, 1000,
                        1.e-12, 1.e-15);

        nerror += test3("epsilon/kp*2*x/w0**2*exp(-(x**2+y**2)/w0**2)*sin(k0*z)",
                        {{"epsilon",0.01},{"kp",3.5},{"w0",5.e-6},{"k0",3.e5}},
                        {"x","y","z"},
                        [=] (Real x, Real y, Real z) -> Real {
                            Real epsilon=0.01, kp=3.5, w0=5.e-6, k0=3.e5;
                            return epsilon/kp*2*x/(w0*w0)*std::exp(-(x*x+y*y)/(w0*w0))*sin(k0*z);
                        },
                        {0.e-6, 0.0, -20.e-6}, {20.e-6, 1.e-10, 20.e-6}, 100,
                        1.e-12, 1.e-15);

        amrex::Print() << "\nMax stack size is " << max_stack_size << "\n";
        if (nerror > 0) {
            amrex::Print() << nerror << " tests failed\n";
            amrex::Abort();
        } else {
            amrex::Print() << "All tests passed\n";
        }
        amrex::Print() << "\n";
    }

    {
        int count = 0;
        int x = 11;
        {
            auto f = [&] (std::string const& s) -> int
            {
                amrex::Print() << count++ << ". Testing \"" << s << "\"\n";
                IParser iparser(s);
                iparser.registerVariables({"x"});
                auto exe = iparser.compileHost<1>();
                return exe(x);
            };
            AMREX_ALWAYS_ASSERT(f("2*(x/3)") == (2*(x/3)));
            AMREX_ALWAYS_ASSERT(f("2*(3/x)") == (2*(3/x)));
            AMREX_ALWAYS_ASSERT(f("2/(x/3)") == (2/(x/3)));
            AMREX_ALWAYS_ASSERT(f("2/(22/x)") == (2/(22/x)));
            AMREX_ALWAYS_ASSERT(f("x/13*5") == ((x/13)*5));
            AMREX_ALWAYS_ASSERT(f("13/x*5") == ((13/x)*5));
            AMREX_ALWAYS_ASSERT(f("x/13/5") == ((x/13)/5));
            AMREX_ALWAYS_ASSERT(f("13/x/5") == ((13/x)/5));

            auto g = [&] (std::string const& s, std::string const& c, int cv) -> int
            {
                amrex::Print() << count++ << ". Testing \"" << s << "\"\n";
                IParser iparser(s);
                iparser.registerVariables({"x"});
                iparser.setConstant(c, cv);
                auto exe = iparser.compileHost<1>();
                return exe(x);
            };
            AMREX_ALWAYS_ASSERT(g("a*(x/3)", "a", 2) == (2*(x/3)));
            AMREX_ALWAYS_ASSERT(g("a*(3/x)", "a", 2) == (2*(3/x)));
            AMREX_ALWAYS_ASSERT(g("a/(x/3)", "a", 2) == (2/(x/3)));
            AMREX_ALWAYS_ASSERT(g("a/(22/x)", "a", 2) == (2/(22/x)));
            AMREX_ALWAYS_ASSERT(g("x/b*5", "b", 13) == ((x/13)*5));
            AMREX_ALWAYS_ASSERT(g("b/x*5", "b", 13) == ((13/x)*5));
            AMREX_ALWAYS_ASSERT(g("x/b/5", "b", 13) == ((x/13)/5));
            AMREX_ALWAYS_ASSERT(g("b/x/5", "b", 13) == ((13/x)/5));

            auto h = [&] (std::string const& s) -> int
            {
                amrex::Print() << count++ << ". Testing \"" << s << "\"\n";
                IParser iparser(s);
                auto exe = iparser.compileHost<0>();
                return exe();
            };
            AMREX_ALWAYS_ASSERT(h("2**10") == 1024);
            AMREX_ALWAYS_ASSERT(h("3^-1") == 1/3);
            AMREX_ALWAYS_ASSERT(h("5^0") == 1);
            AMREX_ALWAYS_ASSERT(h("(-2)**3") == -8);
            AMREX_ALWAYS_ASSERT(h("(-3)**-1") == 1/(-3));
            AMREX_ALWAYS_ASSERT(h("(-5)^0") == 1);

            amrex::Print() << count++ << ". Testing \"a // b\"\n";
            for (int a = -15; a <= 15; ++a) {
                for (int b = -5; b <= 5; ++b) {
                    if (b != 0) {
                        IParser iparser("a//b");
                        iparser.setConstant("a",a);
                        iparser.registerVariables({"b"});
                        auto exe = iparser.compile<1>();
                        AMREX_ALWAYS_ASSERT(exe(b) ==
                                            static_cast<int>(std::floor(double(a)/double(b))));
                    }
                }
            }
        }
        amrex::Print() << "\nAll IParser tests passed\n\n";
    }

    amrex::Finalize();
}
