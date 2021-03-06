cmake_minimum_required(VERSION 3.16)
project(arap)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option (LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui" ON)
set(IMGUIZMO_PATH ${LIBIGL_SOURCE_DIR}/igl/opengl/glfw/imgui/imguizmo)
include(libigl)

# OpenMP
#find_package(OpenMP REQUIRED)

add_subdirectory(external/spdlog)
find_package(spdlog)



# Add your project files
include_directories("include/")
file(GLOB SRC_FILES source/*.cpp)

#add_library(core ${SRC_FILES})
#target_link_libraries(core igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui)

#if(OpenMP_CXX_FOUND)
add_executable(${PROJECT_NAME} "main.cpp"  "include/arapStuff.h" "source/arapStuff.cpp"
        "include/arapDeformer.h" "source/arapDeformer.cpp" "source/precompute.cpp" source/WeightTable.cpp include/WeightTable.h source/MatrixLConstructor.cpp include/MatrixLConstructor.h source/MatrixBConstructor.cpp include/MatrixBConstructor.h)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui spdlog)
#target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui PRIVATE OpenMP::OpenMP_CXX)
#endif()

#add_executable(test "arapTest.cpp" "include/arapStuff.h" "source/arapStuff.cpp"
#        "include/arapDeformer.h" "source/arapDeformer.cpp" "source/precompute.cpp" source/WeightTable.cpp include/WeightTable.h source/MatrixLConstructor.cpp include/MatrixLConstructor.h source/MatrixBConstructor.cpp include/MatrixBConstructor.h )
#target_link_libraries(test PUBLIC igl::core igl::opengl_glfw igl::opengl igl::opengl_glfw_imgui)

