#ifndef STABLE_FLUID
#define STABLE_FLUID

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

/**
 * @brief Data aggregation of current and next frame
 *
 * @tparam T the type of data stored in aggregation
 */
template <typename T> struct TexPair
{
  public:
    T &cur;
    T &nxt;
    TexPair(T &a, T &b) : cur(a), nxt(b)
    {
    }
    void Swap()
    {
        std::swap(cur, nxt);
    }
};

template <typename T> inline T min(T a, T b)
{
    return a < b ? a : b;
}
template <typename T> inline T max(T a, T b)
{
    return a > b ? a : b;
}
/**
 * @brief Get index of 1d array out of 2d index (i,j)
 *
 * @param i index of x-axis
 * @param j index of y-axis
 * @param N width of 2d array
 */
inline int IXY(int i, int j, int N)
{
    return N * j + i;
}

template <typename T, typename SCALAR> inline T lerp(T l, T r, SCALAR t)
{
    return l + t * (r - l);
}

template <typename SCALAR> inline bool check_in_range(SCALAR a, SCALAR b, SCALAR t)
{
    return (a > t) || (b < t);
}

/**
 * @brief Get qf[u, v]
 */
template <typename T, typename SCALAR>
T sample(const int N, const int M, const std::vector<T> &qf, const SCALAR u, const SCALAR v)
{
    int x = static_cast<int>(u);
    x = max<int>(0, min<int>(N - 1, x));
    int y = static_cast<int>(v);
    y = max<int>(0, min<int>(M - 1, y));
    return qf[IXY(x, y, N)];
}

/**
 * @brief Bilinear interpolation
 *
 * @param N width
 * @param M height
 * @param p position to interpolate at
 */
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

/**
 * @brief Locate which point will move to position p in next dt time
 *
 * @param N width
 * @param M height
 * @param vel velocity field
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
VEC backtrace(const int N, const int M, VEC p, SCALAR dt, const std::vector<VEC> &vel)
{
    VEC v1(bilerp<VEC, VEC>(N, M, vel, p));
    VEC p1(p(0) - 0.5 * dt * v1(0), p(1) - 0.5 * dt * v1(1));
    VEC v2(bilerp<VEC, VEC>(N, M, vel, p1));
    VEC p2(p(0) - 0.75 * dt * v2(0), p(1) - 0.75 * dt * v2(1));
    VEC v3(bilerp<VEC, VEC>(N, M, vel, p2));
    p = p + (-1.f) * dt * ((2.f / 9.f) * v1 + (1.f / 3.f) * v2 + (4.f / 9.f) * v3);
    return p * 0.9999f;
}

/**
 * @brief update advection of fluid property qf
 * @param N width
 * @param M height
 * @param vel velocity
 * @param dt time interval
 */
template <typename T, typename VEC, typename SCALAR = typename VEC::Scalar>
void advection(const int N, const int M, const std::vector<VEC> &vel, const std::vector<T> &qf, std::vector<T> &new_qf,
               SCALAR dt, SCALAR damping = 0.9999f)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            VEC p(i, j);
            p += VEC(.5f, .5f);
            p = backtrace<VEC>(N, M, p, dt, vel) * damping;

            new_qf[IXY(i, j, N)] = bilerp<T, VEC>(N, M, qf, p);
        }
    }
}

/**
 * @brief apply global forces and damping
 * @param N width
 * @param M height
 * @param vel velocity field
 * @param g scale of gravity(-9.8f as default)
 * @param dt time interval
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
void apply_force(const int N, const int M, std::vector<VEC> &vel, const SCALAR g = -9.8f, const SCALAR dt = 0.03f)
{
    for (int i = 0; i < N * M; i++)
    {
        vel.at(i) += VEC(0, g * 1.0 * dt);
    }
}

/**
 * @brief Get the divergence of velocity field
 * @param N width
 * @param M height
 * @param vel velocity field
 * @param divergence divergence of velocity field
 */
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

/**
 * @brief Single iteration of gauss-sidel method
 * @param N width
 * @param M height
 * @param pressure pressure field
 * @param divergence divergence of velocity field
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
void pressure_gauss_sidel(const int N, const int M, const std::vector<SCALAR> &divergence,
                          const std::vector<SCALAR> &pressure, std::vector<SCALAR> &new_pressure)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            SCALAR pl = sample<SCALAR, SCALAR>(N, M, pressure, i - 1, j);
            SCALAR pr = sample<SCALAR, SCALAR>(N, M, pressure, i + 1, j);
            SCALAR pb = sample<SCALAR, SCALAR>(N, M, pressure, i, j - 1);
            SCALAR pt = sample<SCALAR, SCALAR>(N, M, pressure, i, j + 1);
            SCALAR diver = divergence[IXY(i, j, N)];
            new_pressure.at(IXY(i, j, N)) = (pl + pr + pb + pt + (-1.f) * diver) * 0.25f;
        }
}

/**
 * @brief Apply pressure to velocity field
 * @param N width
 * @param M height
 * @param vel velocity field
 * @param pressure pressure field
 */
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

template <typename SCALAR> SCALAR smooth_step(SCALAR a, SCALAR x)
{
    SCALAR y = (a - x) / a;
    if (y < 0.0)
        y = 0.0;
    if (y > 1.0)
        y = 1.0;
    SCALAR rst = y * y;
    return rst;
}

/**
 * @brief add round-shape fluid momentum at point [x,y]
 *
 * @param N width
 * @param x x-position
 * @param y y-position
 * @param r radius
 * @param dir direction of added velocity
 * @param dye color
 * @param vel velocity field
 * @param value scale of added color and velocity
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
void add_source(const int N, int x, int y, int r, SCALAR value, VEC dir, std::vector<SCALAR> &dye,
                std::vector<VEC> &vel)
{
    for (int i = -r; i <= r; i++)
        for (int j = -r; j <= r; j++)
        {
            int index = IXY(x + i, y + j, N);
            SCALAR smooth = smooth_step<SCALAR>(r * r, i * i + j * j);
            smooth *= value;
            if (index < 0 || index >= dye.size())
                printf("Error info: index out of range {%d, %d, %d, %d, %d}\n", x, y, i, j, r);
            dye[index] = min(smooth + dye[index], 3.0f);
            vel[index] += dir * smooth * 100.0f;
        }
}

/**
 * @brief Fill scalar field sf into color buffer
 * @param N width
 * @param M height
 * @param sf scalar field to show
 */
template <typename SCALAR>
void fill_color_buffer(const int N, const int M, const std::vector<SCALAR> &sf, unsigned char color_buffer[][3],
                       const SCALAR coe[2])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            int index = IXY(i, j, N);
            SCALAR s = log(abs(sf[index]) * coe[0] + 1.0f);
            SCALAR s3 = coe[1] * s * s * s;
            color_buffer[index][0] = static_cast<unsigned char>(255.f * min(s3, 1.0f));
            color_buffer[index][1] = SCALAR(0);
            color_buffer[index][2] = SCALAR(0);
        }
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void fill_color_buffer(const int N, const int M, const std::vector<VEC> &sf, unsigned char color_buffer[][3],
                       const SCALAR coe[2])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            int index = IXY(i, j, N);
            SCALAR s = log(abs(sf[index](0)) * coe[0] + 1.0f);
            SCALAR s3 = coe[1] * s * s * s;
            color_buffer[index][0] = static_cast<unsigned char>(255.f * min(s3, 1.0f));
            s = log(abs(sf[index](1)) * coe[0] + 1.0f);
            s3 = coe[1] * s * s * s;
            color_buffer[index][1] = static_cast<unsigned char>(255.f * min(s3, 1.0f));
        }
}

#endif