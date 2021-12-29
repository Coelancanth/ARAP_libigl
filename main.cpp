#include <igl/opengl/glfw/Viewer.h>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
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




enum class Mode
{
    kPlacingAnchors = 0,
    kPlacingHandles = 1,
    kDeform = 2,
};
struct State
{
    Eigen::MatrixXd anchors;
    Eigen::MatrixXd handles;
    // used for adding constraints in a least-square way. (see O.Sorkine, Differential Representations for Mesh Processing, p.4)
    std::vector<Triplet> anchor_constraints;
    std::vector<Triplet> handle_constraints;
    //Eigen::MatrixXd anchor_constraints;
    //Eigen::MatrixXd handle_constraints;
    bool isTest = true;
    
    Mode mode = Mode::kPlacingAnchors;
    // bool placing_anchors = true;
    // bool placing_handles = false;
    // bool in_deformation_mode = false;
} state;


// bool add_vertices(Mode mode, State &state, Eigen::RowVector3f &last_mouse,
// igl::opengl::glfw::Viewer &viewer, auto& push_undo)
//{
// }

// bool save_closest_handle_point(Eigen::RowVector3f &last_mouse, long
// &closet_pt_id, State state)
//{
// }

// bool apply_displacement(Eigen::RowVector3f &last_mouse, long closet_pt_id,
// State &state)
//{
// }

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



    /* -------------------------------------------------------------------------- */
    /*                                   Data IO                                  */
    /* -------------------------------------------------------------------------- */
    igl::readOFF("../meshes/bar1.off", Vertices, Faces);
    //igl::readOFF("../meshes/test_mesh.off", Vertices, Faces);
    igl::opengl::glfw::Viewer viewer;

    viewer.data().point_size = kPointSize;
    
    LaplacianPair lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::UNIFORM_WEIGHT);
    

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
            long n_handle = state.handles.rows();
            rhs.bottomRows(n_handle) = state.handles;
            //spdlog::info("n_lhs: {}, n_rhs: {}", solver.rows(), rhs.rows());
            Eigen::MatrixXd rhs_tmp = L_hat_T * rhs;

            Vertices = solver.solve(rhs_tmp);
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
                    LaplacianPair lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::UNIFORM_WEIGHT);
                    lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::UNIFORM_WEIGHT);
                    
                    L_hat = add_constraints(lp.second, state.anchor_constraints, state.handle_constraints);
                    spdlog::info("handles:\n{}", state.handles);
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
        // TODO: lasso selection, idea -> select all vetrices in the range of [x,y]
        // of (current_mouse - last_mouse), then do the following
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

                // TODO: invoke ImGuizmo
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
                    }
                }

                // TODO: invoke ImGuizmo
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
                    }
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
                update();
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
            Eigen::RowVector3f drag_mouse(viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y,
                                          last_mouse(2));
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
            // TODO: Encapsulate them into function
            state.mode = Mode::kDeform;
            spdlog::info("Mode: [Deformation]");
            spdlog::info("Precomputation Begins");
            
            // TODO: Generating matrix is too slow now, need optimization.
            // Precomputation for naive laplacian
            LaplacianPair lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::UNIFORM_WEIGHT);
            lp = calculate_laplacian_matrix(Vertices, Faces, WeightType::UNIFORM_WEIGHT);
            // rhs, laplacian coordinates
            Eigen::MatrixXd delta = lp.first * Vertices;
            long n_vertice = delta.rows();
            long n_anchor = state.anchors.rows();
            long n_handle = state.handles.rows();
            rhs.conservativeResize(n_vertice, 3);
            rhs.topRows(n_vertice) = delta;
            rhs.conservativeResize(n_vertice + n_anchor, 3);
            rhs.bottomRows(n_anchor) = state.anchors;
            //spdlog::info("rhs:\n{}", rhs);
            rhs.conservativeResize(rhs.rows() + n_handle, 3);

            
            
            L_hat = add_constraints(lp.second, state.anchor_constraints, state.handle_constraints);
            L_hat_T = Eigen::SparseMatrix<double, Eigen::RowMajor>(L_hat.transpose());
            // Eigen::SparseMatrix<double> M = Eigen::SparseMatrix<double>(L_hat.transpose()) * L_hat;
            Eigen::SparseMatrix<double> M = L_hat_T * L_hat;

            

            // Cholesky decomposition, only need to do once, before entering into deformation mode.
            spdlog::info("Factorization Begins");
            solver.compute(M);
            spdlog::info("Factorization Ends");
            if (solver.info() != Eigen::Success)
            {
                throw std::runtime_error("Deformation failed! Possibly non semi-positive definite matrix!");
            }
            spdlog::info("Precomputation Ends");
            return true;

            break;
        }
        // Will just trigger a update
        case 'U':
        case 'u': {
            break;
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