#ifndef GL_API
#define GL_API

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>

class FPS_Counter
{
  private:
    int count;
    double start_time;

  public:
    void ResetCounter()
    {
        count = 0;
        start_time = glfwGetTime();
    }

    double GetTime()
    {
        return glfwGetTime() - start_time;
    }

    void AddCount()
    {
        count++;
    }

    int GetFPS()
    {
        return static_cast<int>(count / (glfwGetTime() - start_time));
    }
};

namespace glapi
{

bool glew_init();

GLFWwindow *gl_create_window(const int width, const int height, char *window_name);

void gl_update_window(GLFWwindow *window);

void gl_end(GLFWwindow *window);

void gl_allocate_gltex(GLuint &texName);

void gl_draw_tex2d(const GLuint texName);

} // namespace glapi

/* Implementation of method with templates*/

template <typename VEC, typename SCALAR = typename VEC::Scalar>
inline VEC world_to_screen(const int N, const int M, int x, int y, SCALAR halfdx, SCALAR halfdy)
{
    VEC p = SCALAR(2.0) * VEC(SCALAR(x) / N + halfdx, SCALAR(y) / M + halfdy) + VEC(-1.0, -1.0);
    return p;
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void DrawArray(const int N, const int M, const std::vector<SCALAR> &img)
{
    glClearColor(0.f, 0.f, 0.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);
    glPointSize(1.);
    glBegin(GL_POINTS);
    for (int i = 0; i < N * M; i++)
    {
        int x = i % N;
        int y = i / N;
        VEC p = world_to_screen<VEC>(N, M, x, y, 1.0 / SCALAR(N), 1.0 / SCALAR(M));
        typename VEC::Scalar c = img[i];
        glColor4f(c, 0., 0., (typename VEC::Scalar)x / N);
        glVertex3f(p(0), p(1), 0);
    }
    glEnd();
}

#endif
