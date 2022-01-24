#ifndef STABLE_FLUID
#define STABLE_FLUID
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
T sample(const int N, const int M, const std::vector<T> &qf, const SCALAR u, const SCALAR v)
{
    int x = static_cast<int>(u);
    x = max(0, min(N - 1, x));
    int y = static_cast<int>(v);
    y = max(0, min(M - 1, y));
    return qf[IXY(x, y, N)];
}

template <typename T, typename VEC, typename SCALAR = typename VEC::Scalar>
T bilerp(const int N, const int M, const std::vector<T> &qf, const VEC &p)
{
    SCALAR s = p(0) - 0.5f;
    SCALAR t = p(1) - 0.5f;
    SCALAR iu = floor(s);
    SCALAR iv = floor(t);
    SCALAR fu = s - iu;
    SCALAR fv = t - iv;
    T a = sample<T, SCALAR>(N, M, qf, iu, iv);
    T b = sample<T, SCALAR>(N, M, qf, iu + 1, iv);
    T c = sample<T, SCALAR>(N, M, qf, iu, iv + 1);
    T d = sample<T, SCALAR>(N, M, qf, iu + 1, iv + 1);
    return lerp<T, SCALAR>(lerp<T, SCALAR>(a, b, fu), lerp<T, SCALAR>(c, d, fu), fv);
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
VEC backtrace(const int N, const int M, VEC p, SCALAR dt, const std::vector<VEC> &vel)
{
    VEC v1(bilerp<VEC, VEC>(N, M, vel, p));
    VEC p1(p(0) - 0.5 * dt * v1(0), p(1) - 0.5 * dt * v1(1));
    VEC v2(bilerp<VEC, VEC>(N, M, vel, p1));
    VEC p2(p(0) - 0.75 * dt * v2(0), p(1) - 0.75 * dt * v2(1));
    VEC v3(bilerp<VEC, VEC>(N, M, vel, p2));
    p = p + (-1.f) * dt * ((2.f / 9.f) * v1 + (1.f / 3.f) * v2 + (4.f / 9.f) * v3);
    return p;
}

template <typename T, typename VEC, typename SCALAR = typename VEC::Scalar>
void advection(const int N, const int M, const std::vector<VEC> &vel, const std::vector<T> &qf, std::vector<T> &new_qf,
               SCALAR dt)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            VEC p(i, j);
            p += VEC(.5f, .5f);
            p = backtrace<VEC>(N, M, p, dt, vel);

            new_qf[IXY(i, j, N)] = bilerp<T, VEC>(N, M, qf, p);
        }
    }
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void apply_force(const int N, const int M, std::vector<VEC> &vel, const SCALAR &damping, const SCALAR g,
                 const SCALAR dt)
{
    for (int i = 0; i < N * M; i++)
    {
        vel.at(i) += VEC(0, g * 1.0 * dt);
        vel.at(i) *= damping;
    }
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void get_divergence(const int N, const int M, const std::vector<VEC> &vel, std::vector<SCALAR> &divergence)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            SCALAR vl = sample<VEC, SCALAR>(N, M, vel, i - 1, j)(0);
            SCALAR vr = sample<VEC, SCALAR>(N, M, vel, i + 1, j)(0);
            SCALAR vb = sample<VEC, SCALAR>(N, M, vel, i, j - 1)(1);
            SCALAR vt = sample<VEC, SCALAR>(N, M, vel, i, j + 1)(1);
            SCALAR vc_u = sample<VEC, SCALAR>(N, M, vel, i, j)(0);
            SCALAR vc_v = sample<VEC, SCALAR>(N, M, vel, i, j)(1);
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

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void subtract_gradient(const int N, const int M, std::vector<VEC> &vel, const std::vector<SCALAR> &pressure)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            SCALAR pl = sample<SCALAR, SCALAR>(N, M, pressure, i - 1, j);
            SCALAR pr = sample<SCALAR, SCALAR>(N, M, pressure, i + 1, j);
            SCALAR pb = sample<SCALAR, SCALAR>(N, M, pressure, i, j - 1);
            SCALAR pt = sample<SCALAR, SCALAR>(N, M, pressure, i, j + 1);
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

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void add_source(const int N, int x, int y, int r, SCALAR value, std::vector<SCALAR> &dye, std::vector<VEC> &vel)
{
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
void fill_color_s(const int N, const int M, const std::vector<SCALAR> &sf, std::vector<SCALAR> &color_buffer)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            SCALAR s = log(sf[IXY(i, j, N)] * 0.25f + 1.0f);
            SCALAR s3 = s * s * s;
            color_buffer[IXY(i, j, N)] = abs(1.5f * s);
        }
    }
}

#endif