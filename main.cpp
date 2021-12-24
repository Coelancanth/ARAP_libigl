#include <igl/opengl/glfw/Viewer.h>
#include "spdlog/spdlog.h"
#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <stack>
#include <Eigen/Core>

Eigen::MatrixXd Vertices, Color;
Eigen::MatrixXi Faces;

struct State
{
    Eigen::MatrixXd anchors, handles;
    //bool placing_anchors = true;
    //bool placing_handles = false;
    //bool in_deformation_mode = false;
}   state;

enum Mode
{
    PLACING_ANCHORS = 0,
    PLACING_HANDLES = 1,
    DEFORM = 2,
}   mode = PLACING_ANCHORS;

//bool add_vertices(Mode mode, State &state, Eigen::RowVector3f &last_mouse, igl::opengl::glfw::Viewer &viewer, auto& push_undo)
//{
//}

//bool save_closest_handle_point(Eigen::RowVector3f &last_mouse, long &closet_pt_id, State state)
//{
//}

//bool apply_displacement(Eigen::RowVector3f &last_mouse, long closet_pt_id, State &state)
//{
//}

int main(int argc, char *argv[])
{

    // ================================================================
    // DataIO/ Variables
    // mouse position, used for picking vertices
    // ================================================================
    Eigen::RowVector3f last_mouse;
    long closet_pt_id = -1;
    igl::readOFF("../meshes/bar1.off", Vertices, Faces);
    igl::opengl::glfw::Viewer viewer;
    
    viewer.data().point_size = 10;

    // =================================================================
    // Lambda functions
    // =================================================================
    
    // Undo Management
    std::stack<State> undo_stack, redo_stack;

    const auto push_undo = [&](State & _state = state)
    {
        undo_stack.push(_state);
        // clear
        redo_stack = std::stack<State>();
    };

    const auto undo = [&]()
    {
        if(!undo_stack.empty())
        {
            redo_stack.push(state);
            state = undo_stack.top();
            undo_stack.pop();
        }
    };
    
    const auto redo = [&]()
    {
        if(!redo_stack.empty())
        {
            undo_stack.push(state);
            state = redo_stack.top();
            redo_stack.pop();
        }
    };
    
    const auto & update = [&]()
    {
        // predefined colors
        const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
        const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
        const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
        const Eigen::RowVector3d green(0.2, 0.6, 0.3);
        if (mode == PLACING_ANCHORS)
        {
            viewer.data().set_vertices(Vertices);
            viewer.data().set_colors(blue);
            viewer.data().set_points(state.anchors, orange);
        }
        
        viewer.data().compute_normals();
        
    };



    // =================================================================
    // mouse interaction
    // =================================================================
    viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer&, int, int) -> bool
    {
        last_mouse = Eigen::RowVector3f(
            viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
        // TODO: lasso selection
        if (mode != DEFORM)
        {
            // add vertex by finding closest point on mesh to mouse position
            int face_id;
            Eigen::Vector3f bary_coor;
            if(igl::unproject_onto_mesh(
                last_mouse.head(2),
                viewer.core().view,
                viewer.core().proj,
                viewer.core().viewport,
                Vertices, Faces, 
                face_id, bary_coor))
                {
                    long max_coeff;
                    bary_coor.maxCoeff(&max_coeff);
                    Eigen::RowVector3d new_coor= Vertices.row(Faces(face_id, max_coeff));
            
                    if (mode == PLACING_ANCHORS)
                    {
                        if(state.anchors.size() == 0 || (state.anchors.rowwise() - new_coor).rowwise().norm().minCoeff() > 0)
                        {
                            push_undo();
                            state.anchors.conservativeResize(state.anchors.rows()+1, 3);
                            // Snap to closest vertex on hit face
                            state.anchors.row(state.anchors.rows()-1) = new_coor;
                            //spdlog::info("add {} to the anchors", new_coor);

                        }
                    }

                    if (mode == PLACING_HANDLES)
                    {
                        if(state.handles.size() == 0 || (state.handles.rowwise() - new_coor).rowwise().norm().minCoeff() > 0)
                        {
                            push_undo();
                            state.handles.conservativeResize(state.handles.rows()+1, 3);
                            // Snap to closest vertex on hit face
                            state.handles.row(state.handles.rows()-1) = new_coor;
                            //spdlog::info("add {} to the handles", new_coor);
                        }
                    }
                    update();
                    return true;
                }
        }else{
            // save closet point's coordinates for later use of computing displacement
            //return save_closest_handle_point(last_mouse);
            Eigen::MatrixXf control_pts;
            igl::project(
                Eigen::MatrixXf(state.handles.cast<float>()),
                viewer.core().view,
                viewer.core().proj, viewer.core().viewport, control_pts);
            Eigen::VectorXf distance = (control_pts.rowwise()-last_mouse).rowwise().norm();
            closet_pt_id = (distance.minCoeff(&closet_pt_id) < 30)? closet_pt_id : -1;
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
    
    viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &, int, int) -> bool
    {
        //apply_displacement(last_mouse, closet_pt_id, state);
        if(closet_pt_id != -1)
        {
            Eigen::RowVector3f drag_mouse(
                viewer.current_mouse_x,
                viewer.core().viewport(3) - viewer.current_mouse_y,
                last_mouse(2));
            Eigen::RowVector3f drag_scene, last_scene;
            igl::unproject(
                drag_mouse, 
                viewer.core().view,
                viewer.core().proj,
                viewer.core().viewport, 
                drag_scene);
            igl::unproject(
                last_mouse, 
                viewer.core().view,
                viewer.core().proj,
                viewer.core().viewport, 
                last_scene);
            state.handles.rowwise() += (drag_scene - last_scene).cast<double>();
            last_mouse = drag_mouse;
            update();
            return true;
        }
        return false;
    };
    
    viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int) -> bool
    {
        closet_pt_id = -1;
        return false;
    };
    // =================================================================
    // keyboard interactions
    // =================================================================
    
    
    

    // Plot the mesh
    viewer.data().set_mesh(Vertices, Faces);
    viewer.data().set_face_based(true);
    update();
    viewer.launch();
}