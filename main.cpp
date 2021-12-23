#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXd Vertices;
Eigen::MatrixXi Faces;

int main(int argc, char *argv[])
{
    igl::readOFF("../meshes/bar1.off", Vertices, Faces);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(Vertices, Faces);
  viewer.data().set_face_based(true);
  viewer.launch();
}