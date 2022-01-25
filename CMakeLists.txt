cmake_minimum_required(VERSION 3.17)

set(CMAKE_BUILD_TYPE "Release")

project(fluid_simulation)

find_package(glfw3 REQUIRED PATHS ${CMAKE_SOURCE_DIR}/lib)
if(glfw3_FOUND)
    message("glfw3 loaded: ${glfw3_DIR}")
endif()
find_package(GLEW REQUIRED PATHS ${CMAKE_SOURCE_DIR}/lib)
if(GLEW_FOUND)
    message("GLEW loaded: ${GLEW_DIR}")
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

target_link_directories(demo
    PRIVATE
)

target_link_libraries(demo
    PRIVATE
    glfw
    GLEW::GLEW
    #Eigen3::Eigen
)

add_custom_command(TARGET demo POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/lib/bin"
        $<TARGET_FILE_DIR:demo>)