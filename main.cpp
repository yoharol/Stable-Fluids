//#include <Eigen/Core>
#include <Eigen/Core>
#include <random>

#include "gl_api.h"
#include "stable_fluid.h"

using VEC = Eigen::Vector2f;
using SCALAR = typename VEC::Scalar;

int main()
{

    const int N = 256; // width
    const int M = 256; // height
    const int NN = N * M;

    SCALAR dt = 0.03f;
    SCALAR dx = 1.f;
    SCALAR damping = 0.9999f;
    int p_solver_iters = 80;
    SCALAR g = -9.8f;

    std::vector<VEC> vel(NN, VEC::Zero());
    std::vector<VEC> new_vel(NN, VEC::Zero());
    std::vector<SCALAR> divergence(NN, 0);
    std::vector<SCALAR> pressure(NN, 0);
    std::vector<SCALAR> new_pressure(NN, 0);
    std::vector<SCALAR> dye(NN, 0);
    std::vector<SCALAR> new_dye(NN, 0);
    // std::vector<SCALAR> color_buffer(NN, 0);
    unsigned char color_buffer[N * M];

    TexPair<std::vector<VEC>> vel_tex(vel, new_vel);
    TexPair<std::vector<SCALAR>> pressure_tex(pressure, new_pressure);
    TexPair<std::vector<SCALAR>> dye_tex(dye, new_dye);

    GLFWwindow *window = glapi::gl_create_window(N, M, "stable fuild");

    FPS_Counter counter;
    counter.ResetCounter();
    printf("FPS:");

    GLuint texName;
    glapi::gl_allocate_gltex(texName);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, N, M, 0, GL_RED, GL_UNSIGNED_BYTE, color_buffer);

    std::mt19937 rndeng(std::random_device{}());
    std::uniform_real_distribution<SCALAR> rnda(-1.0f, 1.0f);
    std::uniform_real_distribution<SCALAR> rndr(static_cast<SCALAR>(N) * 0.08f, static_cast<SCALAR>(N) * 0.23f);

    while (!glfwWindowShouldClose(window))
    {
        advection<VEC, VEC>(N, M, vel_tex.cur, vel_tex.cur, vel_tex.nxt, dt, damping);
        advection<SCALAR, VEC>(N, M, vel_tex.cur, dye_tex.cur, dye_tex.nxt, dt, damping);
        vel_tex.Swap();
        dye_tex.Swap();

        int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        if (state == GLFW_PRESS)
        {
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            int pos_x = static_cast<int>(xpos);
            int pos_y = M - static_cast<int>(ypos);
            int r;
            r = min(pos_x - 1, N - pos_x - 1);
            r = min(r, min(pos_y - 1, M - pos_y - 1));
            r = min(r, rndr(rndeng));
            VEC dir(rnda(rndeng), rnda(rndeng));
            dir.normalize();
            assert(pos_x - r >= 0);
            assert(pos_x + r < N);
            assert(pos_y - r >= 0);
            assert(pos_y + r < M);
            add_source<VEC>(N, pos_x, pos_y, r, 3.0f, dir, dye_tex.cur, vel_tex.cur);
        }

        apply_force<VEC>(N, M, vel_tex.cur, g, dt);

        get_divergence<VEC>(N, M, vel_tex.cur, divergence);
        for (int i = 0; i < p_solver_iters; i++)
        {
            pressure_gauss_sidel<VEC>(N, M, divergence, pressure_tex.cur, pressure_tex.nxt);
            pressure_tex.Swap();
        }
        subtract_gradient<VEC>(N, M, vel_tex.cur, pressure_tex.cur);

        fill_color_buffer(N, M, dye_tex.cur, color_buffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, N, M, 0, GL_RED, GL_UNSIGNED_BYTE, color_buffer);

        glapi::gl_draw_tex2d(texName);

        counter.AddCount();
        if (counter.GetTime() >= 1.0)
        {
            printf("\rFPS: %d", counter.GetFPS());
            counter.ResetCounter();
        }
        // DrawArray<VEC>(N, M, color_buffer);
        glapi::gl_update_window(window);
    }
    glapi::gl_end(window);
}