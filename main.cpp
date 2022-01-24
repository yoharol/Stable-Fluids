#include <Eigen/Core>
#include <gl_api.h>
#include <stable_fluid.h>

using VEC = Eigen::Vector2f;
using SCALAR = typename VEC::Scalar;

template <typename VEC> VEC world_to_screen(int x, int y, typename VEC::Scalar halfdx, int N)
{
    VEC p = typename VEC::Scalar(2.0) * VEC(float(x) / N + halfdx, float(y) / N + halfdx) + VEC(-1.0, -1.0);
    return p;
}

template <typename VEC> void DrawArray(const int &N, const std::vector<typename VEC::Scalar> &img)
{
    glClearColor(0.f, 0.f, 0.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);
    glPointSize(1.);
    glBegin(GL_POINTS);
    for (int i = 0; i < N * N; i++)
    {
        int x = i % N;
        int y = i / N;
        VEC p = world_to_screen<VEC>(x, y, 1.0 / (typename VEC::Scalar(N)), N);
        typename VEC::Scalar c = img[i];
        glColor4f(c, 0., 0., (typename VEC::Scalar)x / N);
        glVertex3f(p(0), p(1), 0);
    }
    glEnd();
}

int main()
{

    const int N = 128;
    const int NN = N * N;

    SCALAR dt = 0.03f;
    SCALAR dx = 1.f;
    SCALAR damping = 0.9999f;
    int p_solver_iters = 60;
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

    GLFWwindow *window = CreateWindow(N, N, "stable fuild");

    // FluidSolver2D<VEC, N> fluid_solver(0.03f, 60);

    int count = 0;
    double lastTime = glfwGetTime();
    printf("FPS:");

    while (!glfwWindowShouldClose(window))
    {
        advection<VEC, VEC>(N, *vel_tex.cur, *vel_tex.cur, *vel_tex.nxt, dt);
        advection<SCALAR, VEC>(N, *vel_tex.cur, *dye_tex.cur, *dye_tex.nxt, dt);
        vel_tex.Swap();
        dye_tex.Swap();

        add_source<VEC>(N, static_cast<int>(N * 0.5f), N - 25, 5, 0.8f, *dye_tex.cur, *vel_tex.cur);
        apply_force<VEC>(N, *vel_tex.cur, damping, g, dt);

        get_divergence<VEC>(N, *vel_tex.cur, divergence);
        for (int i = 0; i < p_solver_iters; i++)
        {
            pressure_gauss_sidel<VEC>(N, divergence, *pressure_tex.cur, *pressure_tex.nxt);
            pressure_tex.Swap();
        }
        subtract_gradient<VEC>(N, *vel_tex.cur, *pressure_tex.cur);

        fill_color_s<SCALAR>(N, *dye_tex.cur, color_buffer);

        count++;
        double currentTime = glfwGetTime();
        if (currentTime - lastTime >= 1.0)
        {
            printf("\rFPS: %d", count);
            lastTime = currentTime;
            count = 0;
        }
        DrawArray<VEC>(N, color_buffer);
        UpdateWindow(window);
    }
    EndGL(window);
}