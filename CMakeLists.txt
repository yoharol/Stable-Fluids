cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE "Release")

project(fluid_simulation)

find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(demo)

target_sources(demo
    PRIVATE
    ${CMAKE_SOURCE_DIR}/main.cpp
    ${CMAKE_SOURCE_DIR}/gl_api.cpp
)

target_include_directories(demo
    PRIVATE
    "${CMAKE_SOURCE_DIR}"
)

target_link_libraries(demo
    PRIVATE
    glfw
    GLEW::GLEW
    Eigen3::Eigen
)