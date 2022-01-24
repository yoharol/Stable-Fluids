This project applies the most basic fluid simulation:

- All physics properties are defined in centered-grid
- Linear equations are solved with Gauss-Seidel method
- Single thread

![thumbnail](thumbnail.png)

## Compile

Build Tool: Cmake

External libraries required:

```cmake
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(Eigen3 REQUIRED)
```

(To build this project on windows, [vcpkg]([microsoft/vcpkg: C++ Library Manager for Windows, Linux, and MacOS (github.com)](https://github.com/Microsoft/vcpkg)) is recommended.)

## Future Plan

- Add GPU-based programming to improve performance
- Include vortex method

