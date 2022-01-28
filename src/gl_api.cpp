#include "gl_api.h"

/*bool glapi::glew_init()
{
    // Load glad
    if (glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return false;
    }
    return true;
}*/

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

    // if (!glew_init())
    //    return NULL;

    glViewport(0, 0, width, height);

    return window;
}

void glapi::gl_allocate_gltex(GLuint &texName)
{
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

void glapi::gl_draw_tex2d(const GLuint texName, const double texpos[4][2])
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glBindTexture(GL_TEXTURE_2D, texName);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(texpos[0][0], texpos[0][1], 0.0);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(texpos[1][0], texpos[1][1], 0.0);
    glTexCoord2f(1.0, 1.0);
    glVertex3f(texpos[2][0], texpos[2][1], 0.0);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(texpos[3][0], texpos[3][1], 0.0);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

unsigned char glapi::gl_get_key(GLFWwindow *window, std::map<unsigned char, unsigned char> &key_binding)
{
    for (auto it : key_binding)
    {
        if (glfwGetKey(window, it.first) == GLFW_PRESS)
            return it.second;
    }
    return 0;
}
