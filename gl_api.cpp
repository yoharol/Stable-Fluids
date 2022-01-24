#include <gl_api.h>

void InitGlew()
{
    // Load glew
    if (glewInit() != GLEW_OK)
    {
        std::cerr << "GLEW init failed. It might be caused if the current context is not set.\n";
    }
    else
        printf("GLFW and GLEW successfully initialized with OpenGL verison %s\n", glGetString(GL_VERSION));
}

GLFWwindow *CreateWindow(const int width, const int height, char *window_name)
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

    InitGlew();

    glViewport(0, 0, width, height);

    return window;
}

void UpdateWindow(GLFWwindow *window)
{
    glfwSwapBuffers(window);
    glfwPollEvents();
}

void EndGL(GLFWwindow *window)
{
    glfwDestroyWindow(window);
    glfwTerminate();
}
