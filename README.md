This project implements the most basic fluid simulation:

- All physics properties are defined in centered-grid
- Linear equations are solved with Gauss-Seidel method
- Single thread but good performance

![thumbnail](thumbnail.png)

## Run this program

Download [windows build](https://github.com/yoharol/Stable-Fluids/releases/tag/v1.0)

Click on screen to add random momentum

Press key to set different modes(which physical quantity to show):

| key | mode     |
| --- | -------- |
| 'd' | dye      |
| 'v' | velocity |
| 'p' | pressure |
| 'r' | reset    |

## Compile and Run

External package required:

- GLFW
- OpenGL
- Eigen3

Build tool: CMake

## Future Plan
- [x] MAC grid
- [x] Vortex Confinement
- [ ] Better large-scale linear equation solver
- [ ] Calculate on GPU

