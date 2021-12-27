#include "../include/precompute.h"


Eigen::SparseMatrix<double> calculate_laplacian_matrix(
    const Eigen::MatrixXd &vertices, 
    const Eigen::MatrixXi &faces, 
    const WeightType &weight_type)
    {
        Eigen::SparseMatrix<double> L_s;
        Eigen::SparseMatrix<double> adjacency_matrix;
        Eigen::SparseVector<double> adjacency_matrix_sum;
        Eigen::SparseMatrix<double> adjacency_matrix_diagonal;

        if (weight_type == WeightType::UNIFORM_WEIGHT)
        {
            igl::adjacency_matrix(faces, adjacency_matrix);
            igl::sum(adjacency_matrix, 1, adjacency_matrix_sum);
            igl::diag(adjacency_matrix_sum, adjacency_matrix_diagonal);
            //L_s = adjacency_matrix - adjacency_matrix_diagonal;
            L_s = adjacency_matrix_diagonal - adjacency_matrix;
            return L_s;
        }
        
        // TODO: implement cotan_weight
        //if (weight_type == WeightType::COTAN_WEIGHT)
        //{

        //}
        return L_s;
            
    }
