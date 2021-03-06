#include <igl/opengl/glfw/Viewer.h>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/stopwatch.h"
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <stack>
#include <Eigen/Core>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/opengl/glfw/imgui/ImGuizmoPlugin.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui_internal.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imguizmo/ImGuizmo.h>
#include "precompute.h"
#include <stdexcept>
#include <Eigen/SparseCholesky>
#include <memory>
#include "arapDeformer.h"
#include <cmath>


enum class Mode
{
    kPlacingAnchors = 0,
    kPlacingHandles = 1,
    kDeform = 2,
    kSetRotationAxis = 3,
};
struct State
{
    Eigen::MatrixXd anchors;
    Eigen::MatrixXd handles;
    Eigen::MatrixXd axis;

    // used for adding constraints in a least-square way. (see O.Sorkine, Differential Representations for Mesh Processing, p.4)
    std::vector<Triplet> anchor_constraints;
    std::vector<Triplet> handle_constraints;
    //Eigen::MatrixXd anchor_constraints;
    //Eigen::MatrixXd handle_constraints;

    // used for indexing e.g. L matrix
    std::vector<int> anchorIndices;
    std::vector<int> handleIndices;

    bool isTest = true;

    Mode mode = Mode::kPlacingAnchors;
    // bool placing_anchors = true;
    // bool placing_handles = false;
    // bool in_deformation_mode = false;
} state;


Eigen::Matrix3d calculate_rodrigue_matrix(double angle, Eigen::Vector3d unit_vector)
{   
        Eigen::Matrix3d m;
        double c = cos(angle);
        double s = sin(angle);
        double x = unit_vector.x();
        double y = unit_vector.y();
        double z = unit_vector.z();

        m <<  c + x*x*(1-c), x*y*(1-c)-z*s, x*z*(1-c)+y*s, 
              y*x*(1-c)+z*s, c+y*y*(1-c), y*z*(1-c)-x*s,
              z*x*(1-c)-y*s, z*y*(1-c)+x*s, c+z*z*(1-c);
        return m;
}

int main(int argc, char *argv[])
{

    /* -------------------------------------------------------------------------- */
    /*                                  Variables                                 */
    /* -------------------------------------------------------------------------- */
    Eigen::MatrixXd Vertices;
    Eigen::MatrixXd Color;
    Eigen::MatrixXi Faces;
    const int kPointSize = 10;
    long closet_pt_id = -1;
    Eigen::RowVector3f last_mouse;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double, Eigen::RowMajor> L_hat;
    Eigen::SparseMatrix<double, Eigen::RowMajor> L_hat_T;
    Eigen::MatrixXd rhs;
    
    Eigen::MatrixXd Vertices_new;

    // arap stuff
    std::shared_ptr<arapDeformer> deformer;
    std::vector<int> fixedPts;
    std::vector<Eigen::Vector3d> fixedPositions;
    std::vector<Eigen::Vector3d> axis_points;

    bool constraintsChange = true;


    /* -------------------------------------------------------------------------- */
    /*                                   Data IO                                  */
    /* -------------------------------------------------------------------------- */
    igl::readOFF("../meshes/bar1.off", Vertices, Faces);
    //igl::readOFF("../meshes/dino.off", Vertices, Faces);
    //igl::readOFF("../meshes/test_mesh.off", Vertices, Faces);
    igl::opengl::glfw::Viewer viewer;

    // register the ImGuizmo with viewer
    // can't use more than one imgui plugins now

    //igl::opengl::glfw::imgui::ImGuizmoPlugin imguizmo;
    //viewer.plugins.push_back(&imguizmo);
    //// Initialize ImGuizmo at mesh centroid
    //imguizmo.T.block(0,3,3,1) = 
        //0.5*(Vertices.colwise().maxCoeff() + Vertices.colwise().minCoeff()).transpose().cast<float>();
    //const Eigen::Matrix4f T0 =imguizmo.T;
      //imguizmo.callback = [&](const Eigen::Matrix4f & T)
  //{
    //const Eigen::Matrix4d TT = (T*T0.inverse()).cast<double>().transpose();
    //viewer.data().set_vertices(
      //(Vertices.rowwise().homogeneous()*TT).rowwise().hnormalized());
    //viewer.data().compute_normals();
  //};

    spdlog::info("load complete.");
    spdlog::info("Vertices: {}x{}", Vertices.rows(), Vertices.cols());
    spdlog::info("Faces: {}x{}", Faces.rows(), Faces.cols());
    

    viewer.data().point_size = 10;


    viewer.data().point_size = kPointSize;
    
    //LaplacianPair lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::kUniformWeight);
    deformer = std::make_shared<arapDeformer>(Faces, Vertices);

    /* -------------------------------------------------------------------------- */
    /*                              Lambda functions                              */
    /* -------------------------------------------------------------------------- */

    // Undo Management
    std::stack<State> undo_stack;
    std::stack<State> redo_stack;

    const auto push_undo = [&](State &_state = state) {
        undo_stack.push(_state);
        // clear
        redo_stack = std::stack<State>();
    };

    const auto undo = [&]() {
        if (!undo_stack.empty())
        {
            redo_stack.push(state);
            state = undo_stack.top();
            undo_stack.pop();
        }
    };

    const auto redo = [&]() {
        if (!redo_stack.empty())
        {
            undo_stack.push(state);
            state = redo_stack.top();
            redo_stack.pop();
        }
    };

    const auto &update = [&]() {
        // predefined colors
        const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
        const Eigen::RowVector3d red(1.0, 0, 0);
        const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
        const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
        const Eigen::RowVector3d green(0.2, 0.6, 0.3);
        const Eigen::RowVector3d grey(0.8, 0.8, 0.8);
        
        if (state.mode == Mode::kDeform)
        {
//            long n_handle = state.handles.rows();
//            rhs.bottomRows(n_handle) = state.handles;
//            //spdlog::info("n_lhs: {}, n_rhs: {}", solver.rows(), rhs.rows());
//            Eigen::MatrixXd rhs_tmp = L_hat_T * rhs;
//
//            Vertices = solver.solve(rhs_tmp);
            //


            // both anchor and handles are considered to be constraints.


            fixedPositions.clear();
            int chi; // constraint - handle - i
            for (chi = 0; chi < state.handles.rows(); chi++)
            {
                //spdlog::info ("processing vertex {} as constraint handle", state.handleIndices[chi]);
                //spdlog::info ("its content: {}" , state.handles.row(chi));

                fixedPositions.push_back(state.handles.row(chi));
            }
            int cai;
            for (cai = 0; cai < state.anchors.rows(); cai++)
            {
                //spdlog::info ("processing vertex {} as constraint anchor", state.anchorIndices[cai]);
                //spdlog::info ("its content: {}" , state.anchors.row(cai));

                fixedPositions.push_back(state.anchors.row(cai));

            }

            //deformer->setVertices(Vertices);
            deformer->updateConstraints(fixedPositions);

            //spdlog::stopwatch stopwatch;
            deformer->compute(1);
            //spdlog::info("deformation costs: {}", stopwatch);
            Vertices = deformer->getVertices();
            //std::cout << Vertices << std::endl;




        }

        viewer.data().set_vertices(Vertices);
        viewer.data().set_colors(grey);
        viewer.data().set_points(state.anchors, red);
        viewer.data().add_points(state.handles, orange);

        viewer.data().compute_normals();
    };

    /* -------------------------------------------------------------------------- */
    /*                            Menu plugin                                     */
    /* -------------------------------------------------------------------------- */
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    // Customize the menu
    double doubleVariable = 0.1f; // Shared between two menus

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");

            // ... or using a custom callback
            static bool boolVariable = true;
            if (ImGui::Checkbox("bool", &boolVariable))
            {
                // do something
                std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
            }

            // Expose an enumeration type
            enum Orientation
            {
                Up = 0,
                Down,
                Left,
                Right
            };
            static Orientation dir = Up;
            ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

            // We can also use a std::vector<std::string> defined dynamically
            static int num_choices = 3;
            static std::vector<std::string> choices;
            static int idx_choice = 0;
            if (ImGui::InputInt("Num letters", &num_choices))
            {
                num_choices = std::max(1, std::min(26, num_choices));
            }
            if (num_choices != (int)choices.size())
            {
                choices.resize(num_choices);
                for (int i = 0; i < num_choices; ++i)
                    choices[i] = std::string(1, 'A' + i);
                if (idx_choice >= num_choices)
                    idx_choice = num_choices - 1;
            }
            ImGui::Combo("Letter", &idx_choice, choices);

            // Usages
            if (ImGui::Button("Print Usage", ImVec2(-1, 0)))
            {
                std::cout << R"(
                    Q,q               Switch to [place anchor points] mode
                    W,w               Switch to [place handle points] mode
                    E,e               Switch to [deformation] mode
                    R,r               Rotate 45 degrees
                    U,u               Update deformation (run another iteration of solver)
                    Ctrl+Z            Undo
                    Ctrl+Shift+Z      Redo
                )";
            }
            // Print useful information for debugging
            if (ImGui::Button("Print Information", ImVec2(-1, 0)))
            {
                if (state.isTest == true)
                {
                     //Eigen::Matrix4f m;
                     //m << 1, 2, 3, 4,
                         //5, 6, 7, 8,
                         //9, 10,11,12,
                         //13,14,15,16;
                     //spdlog::info("m is:\n{}", m);
                     //spdlog::info("m.middleRows(1,2) is:\n{}", m.middleRows(1, 2));

                    //spdlog::info("state.anchor_constraint:\n{}", state.anchor_constraints);
                    Eigen::SparseMatrix<double, Eigen::RowMajor> L_hat;
                    LaplacianPair lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::kUniformWeight);
                    lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::kUniformWeight);
                    
                    L_hat = add_constraints(lp.second, state.anchor_constraints, state.handle_constraints);
                    spdlog::info("handles:\n{}", state.handles);

                    spdlog::info("rotation axis:\n {}", state.axis);
                }
            }
        }
    };

    //// Draw additional windows
    // menu.callback_draw_custom_window = [&]()
    //{
    //// Define next window position + size
    // ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10),
    // ImGuiCond_FirstUseEver); ImGui::SetNextWindowSize(ImVec2(200, 160),
    // ImGuiCond_FirstUseEver); ImGui::Begin( "New Window", nullptr,
    // ImGuiWindowFlags_NoSavedSettings
    //);

    //// Expose the same variable directly ...
    // ImGui::PushItemWidth(-80);
    // ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0,
    // 0, "%.4f"); ImGui::PopItemWidth();

    // static std::string str = "bunny";
    // ImGui::InputText("Name", str);

    // ImGui::End();
    //};

    /* -------------------------------------------------------------------------- */
    /*                              Mouse Interaction                             */
    /* -------------------------------------------------------------------------- */
    viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer &, int, int) -> bool {
        last_mouse = Eigen::RowVector3f(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
        // TODO: lasso selection, idea -> select all vetrices in the range of [x,y] of (current_mouse - last_mouse)
        if (state.mode != Mode::kDeform)
        {
            // add vertex by finding closest point on mesh to mouse position
            int face_id;
            Eigen::Vector3f bary_coor;
            if (igl::unproject_onto_mesh(last_mouse.head(2), viewer.core().view, viewer.core().proj,
                                         viewer.core().viewport, Vertices, Faces, face_id, bary_coor))
            {
                long max_coeff;
                bary_coor.maxCoeff(&max_coeff);
                // NOTE: get vertex idx here
                int idx_v = Faces(face_id, max_coeff);
                Eigen::RowVector3d new_coor = Vertices.row(idx_v);

                if (state.mode == Mode::kPlacingAnchors)
                {
                    if (state.anchors.size() == 0 ||
                        (state.anchors.rowwise() - new_coor).rowwise().norm().minCoeff() > 0)
                    {
                        push_undo();
                        state.anchors.conservativeResize(state.anchors.rows() + 1, 3);
                        // Snap to closest vertex on hit face
                        state.anchors.row(state.anchors.rows() - 1) = new_coor;
                        
                        state.anchor_constraints.emplace_back(state.anchors.rows() - 1, idx_v, 1);
                        state.anchorIndices.push_back(idx_v);
                        constraintsChange = true;

                    }
                }

                if (state.mode == Mode::kPlacingHandles)
                {
                    if (state.handles.size() == 0 ||
                        (state.handles.rowwise() - new_coor).rowwise().norm().minCoeff() > 0)
                    {
                        push_undo();
                        state.handles.conservativeResize(state.handles.rows() + 1, 3);
                        // Snap to closest vertex on hit face
                        state.handles.row(state.handles.rows() - 1) = new_coor;

                        state.handle_constraints.emplace_back(state.handles.rows() - 1, idx_v, 1);
                        state.handleIndices.push_back(idx_v);
                        constraintsChange = true;
                    }
                }

                if (state.mode == Mode::kSetRotationAxis)
                {
                    state.axis.conservativeResize(state.axis.rows() + 1, 3);
                    // Snap to closest vertex on hit face
                    state.axis.row(state.axis.rows() - 1) = new_coor;
                }

                update();
                return true;
            }
        }
        else
        {
            // save closet point's coordinates for later use of computing displacement
            // return save_closest_handle_point(last_mouse);
            // FIXME: only works if control_pts != empty
            Eigen::MatrixXf control_pts;
            igl::project(Eigen::MatrixXf(state.handles.cast<float>()), viewer.core().view, viewer.core().proj,
                         viewer.core().viewport, control_pts);
            Eigen::VectorXf distance = (control_pts.rowwise() - last_mouse).rowwise().norm();
            closet_pt_id = (distance.minCoeff(&closet_pt_id) < 30) ? closet_pt_id : -1;
            if (closet_pt_id != -1)
            {
                last_mouse(2) = control_pts(closet_pt_id, 2);
                push_undo();
                //update();
                // spdlog::info("")
                return true;
            }
        }
        return false;
    };

    // FIXME: implement ImGuizmo
    viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &, int, int) -> bool {
        // apply_displacement(last_mouse, closet_pt_id, state);
        if (closet_pt_id != -1)
        {
            Eigen::RowVector3f drag_mouse(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, last_mouse(2));
            Eigen::RowVector3f drag_scene, last_scene;
            igl::unproject(drag_mouse, viewer.core().view, viewer.core().proj, viewer.core().viewport, drag_scene);
            igl::unproject(last_mouse, viewer.core().view, viewer.core().proj, viewer.core().viewport, last_scene);
            state.handles.rowwise() += (drag_scene - last_scene).cast<double>();
            last_mouse = drag_mouse;
            update();
            return true;
        }
        return false;
    };

    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer &, int, int) -> bool {
        //spdlog::info("mouse_up");
        //update();
        closet_pt_id = -1;
        return false;
    };
    

    /* -------------------------------------------------------------------------- */
    /*                            Keyboard Interaction                            */
    /* -------------------------------------------------------------------------- */

    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod) {
        switch (key)
        {
        // Switch modes
        case 'A':
        case 'a': {
            state.mode = Mode::kSetRotationAxis;
            spdlog::info("Mode: [Setting Rotation Axis]");
            break;
        }
        case 'Q':
        case 'q': {
            state.mode = Mode::kPlacingAnchors;
            spdlog::info("Mode: [Placing Anchors]");
            break;
        }
        case 'W':
        case 'w': {
            state.mode = Mode::kPlacingHandles;
            spdlog::info("Mode: [Placing Handles]");
            break;
        }
        case 'E':
        case 'e': {
            state.mode = Mode::kDeform;
            spdlog::info("Mode: [Deformation]");
            spdlog::info("Precomputation Begins");

            fixedPts.clear();
            int chi; // constraint - handle - i
            for (chi = 0; chi < state.handles.rows(); chi++)
            {
                //spdlog::info ("processing vertex {} as constraint handle", state.handleIndices[chi]);
                fixedPts.push_back(state.handleIndices[chi]);
            }
            int cai;
            for (cai = 0; cai < state.anchors.rows(); cai++)
            {
                //spdlog::info ("processing vertex {} as constraint anchor", state.anchorIndices[cai]);
                fixedPts.push_back(state.anchorIndices[cai]);
            }

            deformer->setConstrainedPoints(fixedPts);
            spdlog::info("Precomputation Ends");

            return true;

            break;
        }
        case 'R':
        case 'r': {
            axis_points.clear();
            if (state.axis.rows() == 2)
            {
                Eigen::Vector3d unit_vector = state.axis.row(0) - state.axis.row(1);
                //spdlog::info("unit_vector: {}", unit_vector.normalized());
                Eigen::Matrix3d rodrigue_matrix = calculate_rodrigue_matrix(45, unit_vector.normalized());
                //spdlog::info("r_matrix_row: {}", rodrigue_matrix);
                //state.handle_constraints *= rodrigue_matrix;
                state.handles.resize(state.handles.rows(),3);
                state.handles *= rodrigue_matrix;
                
                update();
                return true;
            }
        }

        // Will just trigger a update
        case 'U':
        case 'u': {
            fixedPositions.clear();
            int chi; // constraint - handle - i
            for (chi = 0; chi < state.handles.rows(); chi++)
            {
                spdlog::info ("processing vertex {} as constraint handle", state.handleIndices[chi]);
                spdlog::info ("its content: {}" , state.handles.row(chi));

                fixedPositions.push_back(state.handles.row(chi));
            }
            int cai;
            for (cai = 0; cai < state.anchors.rows(); cai++)
            {
                spdlog::info ("processing vertex {} as constraint anchor", state.anchorIndices[cai]);
                spdlog::info ("its content: {}" , state.anchors.row(cai));

                fixedPositions.push_back(state.anchors.row(cai));

            }

            deformer->setVertices(Vertices);
            deformer->updateConstraints(fixedPositions);

            deformer->compute(2);
            Vertices = deformer->getVertices();
            return true;
        }
        default:
            return false;
        }
        update();
        return true;
    };

    // Undo Management
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned char key, int mod) -> bool {
        if (key == 'Z' && (mod & GLFW_MOD_SUPER))
        {
            (mod & GLFW_MOD_SHIFT) ? redo() : undo();
            update();
            return true;
        }
        return false;
    };

    // Plot the mesh
    viewer.data().set_mesh(Vertices, Faces);
    viewer.data().set_face_based(true);
    update();
    viewer.launch();
}