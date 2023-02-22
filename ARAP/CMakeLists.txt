cmake_minimum_required(VERSION 3.16)
project(arap)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
option(LIBIGL_OPENGL "Build target igl::opengl" ON)
option(LIBIGL_IMGUI  "Build target igl::imgui" OFF)
include(libigl)

# Add your project files
AUX_SOURCE_DIRECTORY(src SOURCES)
add_executable(${PROJECT_NAME} main.cpp ${SOURCES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl igl::glfw)
