This project implements the most basic fluid simulation:

- All physics properties are defined in centered-grid
- Linear equations are solved with Gauss-Seidel method
- Single thread but good performance

![thumbnail](thumbnail.png)

## How to play

Download [windows build](https://github.com/yoharol/Stable-Fluids/releases/tag/v1.0)

Click on screen to add random momentum

## Compile and Run

There's no external package required.

```shell
mkdir build
cd ./build
cmake ..
cmake --build .
./demo
```

## Future Plan

- Add GPU-based programming to improve performance
- Include vortex method

