#ifndef FLUID_SOLVER
#define FLUID_SOLVER

#include <vector>

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

#endif