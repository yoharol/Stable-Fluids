cmake_minimum_required(VERSION 3.17)

project(fluid_simulation)

set(CMAKE_BUILD_TYPE "Release")

find_package(glfw3 REQUIRED)
if(glfw3_FOUND)
    message("glfw3 loaded: ${glfw3_DIR}")
endif()

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message("Eigen3 loaded: ${Eigen3_DIR}")
endif()

find_package(GLEW REQUIRED)
if(GLEW_FOUND)
    message("GLEW loaded: ${glfw3_DIR}")
endif()

add_executable(demo)

target_sources(demo
    PRIVATE
    ${CMAKE_SOURCE_DIR}/main.cpp
    ${CMAKE_SOURCE_DIR}/src/gl_api.cpp
)

target_include_directories(demo
    PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
)

target_link_libraries(demo
    PRIVATE
    glfw
    Eigen3::Eigen
    GLEW::GLEW
)
