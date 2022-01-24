#include <Eigen/Core>
#include <gl_api.h>
#include <solver.h>
#include <stable_fluid.h>

using VEC = Eigen::Vector2f;
using SCALAR = typename VEC::Scalar;

int main()
{

    const int N = 256; // width
    const int M = 128; // height
    const int NN = N * M;

    SCALAR dt = 0.03f;
    SCALAR dx = 1.f;
    SCALAR damping = 0.9999f;
    int p_solver_iters = 40;
    SCALAR g = -9.8f;

    std::vector<VEC> vel(NN, VEC::Zero());
    std::vector<VEC> new_vel(NN, VEC::Zero());
    std::vector<SCALAR> divergence(NN, 0);
    std::vector<SCALAR> pressure(NN, 0);
    std::vector<SCALAR> new_pressure(NN, 0);
    std::vector<SCALAR> dye(NN, 0);
    std::vector<SCALAR> new_dye(NN, 0);
    std::vector<SCALAR> color_buffer(NN, 0);

    TexPair<std::vector<VEC>> vel_tex(vel, new_vel);
    TexPair<std::vector<SCALAR>> pressure_tex(pressure, new_pressure);
    TexPair<std::vector<SCALAR>> dye_tex(dye, new_dye);

    GLFWwindow *window = CreateWindow(N, M, "stable fuild");

    FPS_Counter counter;
    counter.ResetCounter();
    printf("FPS:");

    while (!glfwWindowShouldClose(window))
    {
        advection<VEC, VEC>(N, M, *vel_tex.cur, *vel_tex.cur, *vel_tex.nxt, dt);
        advection<SCALAR, VEC>(N, M, *vel_tex.cur, *dye_tex.cur, *dye_tex.nxt, dt);
        vel_tex.Swap();
        dye_tex.Swap();

        add_source<VEC>(N, N * 0.5f, M * 0.8f, 5, 0.8f, *dye_tex.cur, *vel_tex.cur);
        apply_force<VEC>(N, M, *vel_tex.cur, damping, g, dt);

        get_divergence<VEC>(N, M, *vel_tex.cur, divergence);
        for (int i = 0; i < p_solver_iters; i++)
        {
            pressure_gauss_sidel<VEC>(N, M, divergence, *pressure_tex.cur, *pressure_tex.nxt);
            pressure_tex.Swap();
        }
        subtract_gradient<VEC>(N, M, *vel_tex.cur, *pressure_tex.cur);

        fill_color_s<SCALAR>(N, M, *dye_tex.cur, color_buffer);

        counter.AddCount();
        if (counter.GetTime() >= 1.0)
        {
            printf("\rFPS: %d", counter.GetFPS());
            counter.ResetCounter();
        }
        DrawArray<VEC>(N, M, color_buffer);
        UpdateWindow(window);
    }
    EndGL(window);
}