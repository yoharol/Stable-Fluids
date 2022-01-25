#include "gl_api.h"

bool glapi::glad_init()
{
    // Load glad
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return false;
    }
    return true;
}

GLFWwindow *glapi::gl_create_window(const int width, const int height, char *window_name)
{
    glfwInit();
    GLFWwindow *window = glfwCreateWindow(width, height, window_name, NULL, NULL);
    if (!window)
    {
        std::cerr << "Cannot create window\n";
        glfwTerminate();
        return NULL;
    }
    // Introduce window to current context
    glfwMakeContextCurrent(window);

    if (!glad_init())
        return NULL;

    glViewport(0, 0, width, height);

    return window;
}

void glapi::gl_allocate_gltex(GLuint &texName)
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

void glapi::gl_update_window(GLFWwindow *window)
{
    glfwSwapBuffers(window);
    glfwPollEvents();
}

void glapi::gl_end(GLFWwindow *window)
{
    glfwDestroyWindow(window);
    glfwTerminate();
}

void glapi::gl_draw_tex2d(const GLuint texName)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glBindTexture(GL_TEXTURE_2D, texName);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(-1.0, -1.0, 0.0);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(-1.0, 1.0, 0.0);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(1.0, 1.0, 0.0);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(1.0, -1.0, 0.0);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}
