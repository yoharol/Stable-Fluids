//#include <Eigen/Core>
#include <Eigen/Core>
#include <algorithm>
#include <random>

#include "gl_api.h"
#include "stable_fluid.h"

using VEC = Eigen::Vector2f;
using SCALAR = typename VEC::Scalar;

int main()
{
    // width
    const int N = 256;

    // height
    const int M = 256;

    const int NN = N * M;

    // update time delta
    SCALAR dt = 0.03f;

    // interval between each grid
    SCALAR dx = 1.f;

    SCALAR damping = 0.9999f;

    // number of iteration in Gauss-Sidel method
    int p_solver_iters = 80;

    SCALAR g = -9.8f;

    std::vector<VEC> vel(NN, VEC::Zero());
    std::vector<VEC> new_vel(NN, VEC::Zero());
    std::vector<SCALAR> divergence(NN, 0);
    std::vector<SCALAR> pressure(NN, 0);
    std::vector<SCALAR> new_pressure(NN, 0);
    std::vector<SCALAR> dye(NN, 0);
    std::vector<SCALAR> new_dye(NN, 0);

    // color to draw in each frame
    unsigned char color_buffer[N * M][3];

    TexPair<std::vector<VEC>> vel_tex(vel, new_vel);
    TexPair<std::vector<SCALAR>> pressure_tex(pressure, new_pressure);
    TexPair<std::vector<SCALAR>> dye_tex(dye, new_dye);

    // create window
    char window_name[13] = "stable fluid";
    GLFWwindow *window = glapi::gl_create_window(N, M, window_name);
    if (window == NULL)
        return 0;
    glClearColor(0.0, 0.0, 0.0, 1.0);

    // count fps
    FPS_Counter counter;
    counter.ResetCounter();
    printf("FPS:");

    // allocate texture to draw
    GLuint texName;
    glapi::gl_allocate_gltex(texName);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, N, M, 0, GL_RGB, GL_UNSIGNED_BYTE, color_buffer);
    // position of texture in screen
    const double texpos[4][2] = {{-1.0, -1.0}, {-1.0, 1.0}, {1.0, 1.0}, {1.0, -1.0}};

    // initialize random number generator
    std::mt19937 rndeng(std::random_device{}());
    std::uniform_real_distribution<SCALAR> rnda(-1.0f, 1.0f);
    std::uniform_real_distribution<SCALAR> rndr(static_cast<SCALAR>(N) * 0.08f, static_cast<SCALAR>(N) * 0.23f);

    // add source to fluid as a demonstration
    add_source<VEC>(N, N / 2, M / 2, N / 10, 3.0f, -VEC::UnitY(), dye_tex.cur, vel_tex.cur);

    // define keyboard input
    unsigned char input_key = 0;
    unsigned char draw_num = 0;
    // coefficient to be used in fill_color_buffer(stable_fluid.h)
    SCALAR coe_draw[4][3] = {{0.25f, 6.0f}, {0.01f, 1.0f}, {0.01f, 4.0f}};
    std::map<unsigned char, unsigned char> key_binding = {
        {GLFW_KEY_R, 'r'}, {GLFW_KEY_D, 'd'}, {GLFW_KEY_P, 'p'}, {GLFW_KEY_V, 'v'}};

    // time of previous input(to avoid continuous invoking same function)
    double previous_click = 0.;
    double previous_press = 0.;

    // update loop
    while (!glfwWindowShouldClose(window))
    {
        // advection of velocity field and dye field
        advection<VEC, VEC>(N, M, vel_tex.cur, vel_tex.cur, vel_tex.nxt, dt, damping);
        advection<SCALAR, VEC>(N, M, vel_tex.cur, dye_tex.cur, dye_tex.nxt, dt, damping);
        vel_tex.Swap();
        dye_tex.Swap();

        // check mouse input
        int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        if (state == GLFW_PRESS && glfwGetTime() - previous_click > 0.1)
        {
            previous_click = glfwGetTime();

            // add source of fliud at moust position with a random momentum
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            int pos_x = static_cast<int>(xpos);
            int pos_y = M * 2 - static_cast<int>(ypos);
            if (pos_y > M)
            {
                pos_y -= M;
                int r;
                r = min<int>(pos_x - 1, N - pos_x - 1);
                r = min<int>(r, min(pos_y - 1, M - pos_y - 1));
                r = min<int>(r, rndr(rndeng));
                VEC dir(rnda(rndeng), rnda(rndeng));
                dir.normalize();
                assert(pos_x - r >= 0);
                assert(pos_x + r < N);
                assert(pos_y - r >= 0);
                assert(pos_y + r < M);
                add_source<VEC>(N, pos_x, pos_y, r, 3.0f, dir, dye_tex.cur, vel_tex.cur);
            }
        }

        // check keyboard input
        if (glfwGetTime() - previous_press > 0.1)
        {
            input_key = glapi::gl_get_key(window, key_binding);
            if (input_key != 0)
                previous_press = glfwGetTime();
            switch (input_key)
            {
            case 'r':
                // reset all
                std::fill(vel_tex.cur.begin(), vel_tex.cur.end(), VEC::Zero());
                std::fill(pressure_tex.cur.begin(), pressure_tex.cur.end(), SCALAR(0));
                std::fill(dye_tex.cur.begin(), dye_tex.cur.end(), SCALAR(0));
                break;
            case 'd':
                // draw dye
                draw_num = 0;
                break;
            case 'p':
                // draw pressure
                draw_num = 1;
                break;
            case 'v':
                // draw velocity
                draw_num = 2;
                break;
            default:
                break;
            }
        }

        // apply global force to velocity field
        apply_force<VEC>(N, M, vel_tex.cur, g, dt);

        // solve divergence and pressure
        get_divergence<VEC>(N, M, vel_tex.cur, divergence);
        for (int i = 0; i < p_solver_iters; i++)
        {
            pressure_gauss_sidel<VEC>(N, M, divergence, pressure_tex.cur, pressure_tex.nxt);
            pressure_tex.Swap();
        }
        // apply pressure to velocity field
        subtract_gradient<VEC>(N, M, vel_tex.cur, pressure_tex.cur);

        // pick velocity, pressure or dye to draw
        switch (draw_num)
        {
        case 1:
            fill_color_buffer<SCALAR>(N, M, pressure_tex.cur, color_buffer, coe_draw[1]);
            break;
        case 2:
            fill_color_buffer<VEC>(N, M, vel_tex.cur, color_buffer, coe_draw[2]);
            break;
        case 0:
        default:
            fill_color_buffer<SCALAR>(N, M, dye_tex.cur, color_buffer, coe_draw[0]);
            break;
        }
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, N, M, 0, GL_RGB, GL_UNSIGNED_BYTE, color_buffer);
        glapi::gl_draw_tex2d(texName, texpos);

        // count fps
        counter.AddCount();
        if (counter.GetTime() >= 1.0)
        {
            printf("\rFPS: %d", counter.GetFPS());
            counter.ResetCounter();
        }
        glapi::gl_update_window(window);
    }
    glapi::gl_end(window);
}