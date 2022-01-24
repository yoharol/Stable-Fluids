#ifndef FLUID_SOLVER
#define FLUID_SOLVER
#define min(x, y) (x) < (y) ? (x) : (y)
#define max(x, y) (x) > (y) ? (x) : (y)

#include <Eigen/Core>
#include <cmath>
#include <vector>

template <typename T> struct TexPair
{
  public:
    T *cur;
    T *nxt;
    TexPair(T &a, T &b)
    {
        cur = &a;
        nxt = &b;
    }
    void Swap()
    {
        T *tmp = cur;
        cur = nxt;
        nxt = tmp;
    }
};

inline int IXY(int i, int j, int N)
{
    return N * j + i;
}
template <typename T, typename SCALAR> inline T lerp(T l, T r, SCALAR t)
{
    return l + t * (r - l);
}

template <typename T, typename SCALAR>
T sample(const int &N, const std::vector<T> &qf, const SCALAR &u, const SCALAR &v)
{
    int x = static_cast<int>(u);
    x = max(0, min(N - 1, x));
    int y = static_cast<int>(v);
    y = max(0, min(N - 1, y));
    return qf[IXY(x, y, N)];
}

template <typename T, typename VEC> T bilerp(const int &N, const std::vector<T> &qf, const VEC &p)
{
    using SCALAR = typename VEC::Scalar;
    SCALAR s = p(0) - 0.5f;
    SCALAR t = p(1) - 0.5f;
    SCALAR iu = floor(s);
    SCALAR iv = floor(t);
    SCALAR fu = s - iu;
    SCALAR fv = t - iv;
    T a = sample<T, SCALAR>(N, qf, iu, iv);
    T b = sample<T, SCALAR>(N, qf, iu + 1, iv);
    T c = sample<T, SCALAR>(N, qf, iu, iv + 1);
    T d = sample<T, SCALAR>(N, qf, iu + 1, iv + 1);
    return lerp<T, SCALAR>(lerp<T, SCALAR>(a, b, fu), lerp<T, SCALAR>(c, d, fu), fv);
}

template <typename VEC> VEC backtrace(const int &N, VEC p, typename VEC::Scalar dt, const std::vector<VEC> &vel)
{
    VEC v1(bilerp<VEC, VEC>(N, vel, p));
    VEC p1(p(0) - 0.5 * dt * v1(0), p(1) - 0.5 * dt * v1(1));
    VEC v2(bilerp<VEC, VEC>(N, vel, p1));
    VEC p2(p(0) - 0.75 * dt * v2(0), p(1) - 0.75 * dt * v2(1));
    VEC v3(bilerp<VEC, VEC>(N, vel, p2));
    p = p + (-1.f) * dt * ((2.f / 9.f) * v1 + (1.f / 3.f) * v2 + (4.f / 9.f) * v3);
    return p;
}

template <typename T, typename VEC>
void advection(const int &N, const std::vector<VEC> &vel, const std::vector<T> &qf, std::vector<T> &new_qf,
               typename VEC::Scalar dt)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            VEC p(i, j);
            p += VEC(.5f, .5f);
            p = backtrace<VEC>(N, p, dt, vel);

            new_qf[IXY(i, j, N)] = bilerp<T, VEC>(N, qf, p);
        }
    }
}

template <typename VEC>
void apply_force(const int &N, std::vector<VEC> &vel, const typename VEC::Scalar &damping,
                 const typename VEC::Scalar &g, const typename VEC::Scalar &dt)
{
    for (int i = 0; i < N * N; i++)
    {
        vel.at(i) += VEC(0, g * 1.0 * dt);
        vel.at(i) *= damping;
    }
}

template <typename VEC>
void get_divergence(const int &N, const std::vector<VEC> &vel, std::vector<typename VEC::Scalar> &divergence)
{
    using SCALAR = typename VEC::Scalar;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            SCALAR vl = sample<VEC, SCALAR>(N, vel, i - 1, j)(0);
            SCALAR vr = sample<VEC, SCALAR>(N, vel, i + 1, j)(0);
            SCALAR vb = sample<VEC, SCALAR>(N, vel, i, j - 1)(1);
            SCALAR vt = sample<VEC, SCALAR>(N, vel, i, j + 1)(1);
            SCALAR vc_u = sample<VEC, SCALAR>(N, vel, i, j)(0);
            SCALAR vc_v = sample<VEC, SCALAR>(N, vel, i, j)(1);
            if (i == 0)
                vl = -vc_u;
            else if (i == N - 1)
                vr = -vc_u;
            if (j == 0)
                vb = -vc_v;
            else if (j == N - 1)
                vt = -vc_v;
            divergence[IXY(i, j, N)] = (vr - vl + vt - vb) * .5f;
        }
}

template <typename VEC>
void pressure_gauss_sidel(const int &N, const std::vector<typename VEC::Scalar> &divergence,
                          const std::vector<typename VEC::Scalar> &pressure,
                          std::vector<typename VEC::Scalar> &new_pressure)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            SCALAR pl = sample<SCALAR, SCALAR>(N, pressure, i - 1, j);
            SCALAR pr = sample<SCALAR, SCALAR>(N, pressure, i + 1, j);
            SCALAR pb = sample<SCALAR, SCALAR>(N, pressure, i, j - 1);
            SCALAR pt = sample<SCALAR, SCALAR>(N, pressure, i, j + 1);
            SCALAR diver = divergence[IXY(i, j, N)];
            new_pressure.at(IXY(i, j, N)) = (pl + pr + pb + pt + (-1.f) * diver) * 0.25f;
        }
}

template <typename VEC>
void subtract_gradient(const int &N, std::vector<VEC> &vel, const std::vector<typename VEC::Scalar> &pressure)
{
    using SCALAR = typename VEC::Scalar;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            SCALAR pl = sample<SCALAR, SCALAR>(N, pressure, i - 1, j);
            SCALAR pr = sample<SCALAR, SCALAR>(N, pressure, i + 1, j);
            SCALAR pb = sample<SCALAR, SCALAR>(N, pressure, i, j - 1);
            SCALAR pt = sample<SCALAR, SCALAR>(N, pressure, i, j + 1);
            vel.at(IXY(i, j, N)) -= 0.5f * VEC(pr - pl, pt - pb);
        }
}

template <typename SCALAR> SCALAR smooth_step(SCALAR a, SCALAR b, SCALAR x)
{
    SCALAR y = (x - a) / (b - a);
    if (y < 0.0)
        y = 0.0;
    if (y > 1.0)
        y = 1.0;
    SCALAR rst = y * y * (3.0 - 2.0 * y);
    return rst;
}

template <typename VEC>
void add_source(const int &N, int x, int y, int r, typename VEC::Scalar value, std::vector<typename VEC::Scalar> &dye,
                std::vector<VEC> &vel)
{
    using SCALAR = typename VEC::Scalar;
    for (int index = 0; index < (2 * r + 1) * (2 * r + 1); index++)
    {
        int i = (SCALAR)index / (2 * r + 1) - r;
        int j = index % (2 * r + 1) - r;
        SCALAR smooth = smooth_step<SCALAR>(r * r, 0.0, i * i + j * j);
        dye.at(IXY(x + i, y + j, N)) += value * smooth;
        vel.at(IXY(x + i, y + j, N)) += VEC(0, -value * smooth * 100.0f);
    }
}

template <typename SCALAR>
void fill_color_s(const int &N, const std::vector<SCALAR> &sf, std::vector<SCALAR> &color_buffer)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            SCALAR s = log(sf[IXY(i, j, N)] * 0.25f + 1.0f);
            SCALAR s3 = s * s * s;
            color_buffer[IXY(i, j, N)] = abs(1.5f * s);
        }
    }
}

#endif