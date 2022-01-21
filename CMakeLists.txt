cmake_minimum_required(VERSION 3.16)
project(arap)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option (LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui" ON)
set(IMGUIZMO_PATH ${LIBIGL_SOURCE_DIR}/igl/opengl/glfw/imgui/imguizmo)
include(libigl)



# Add your project files
include_directories("include/")
file(GLOB SRC_FILES source/*.cpp)

#add_library(core ${SRC_FILES})
#target_link_libraries(core igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui)

add_executable(${PROJECT_NAME} "main.cpp"  "include/arapStuff.h" "source/arapStuff.cpp"
        "include/arapDeform.h" "source/arapDeform.cpp" "source/precompute.cpp")
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui)

add_executable(test "arapTest.cpp" "include/arapStuff.h" "source/arapStuff.cpp"
        "include/arapDeform.h" "source/arapDeform.cpp" "source/precompute.cpp")
target_link_libraries(test PUBLIC igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui)
