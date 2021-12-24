cmake_minimum_required(VERSION 3.16)
project(arap)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option (LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui" ON)
set(IMGUIZMO_PATH ${LIBIGL_SOURCE_DIR}/igl/opengl/glfw/imgui/imguizmo)
include(libigl)


# Add your project files
file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui)
