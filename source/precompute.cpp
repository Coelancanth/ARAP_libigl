#include "../include/precompute.h"
#include <Eigen/src/SparseCore/SparseMatrixBase.h>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"


LaplacianPair calculate_laplacian_matrix(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces,
                                         const WeightType &weight_type)
{
    Eigen::SparseMatrix<double> L_s;

    if (weight_type == WeightType::UNIFORM_WEIGHT)
    {
        Eigen::SparseMatrix<double> A;
        igl::adjacency_matrix(faces, A);

        Eigen::SparseVector<double> A_sum;
        igl::sum(A, 1, A_sum);

        Eigen::SparseMatrix<double> D;
        igl::diag(A_sum, D);

        // L_s = adjacency_matrix - adjacency_matrix_diagonal;
        // Eigen::SparseMatrix<double> L;
        Eigen::MatrixXd L;
        L_s = D - A;
        L = Eigen::MatrixXd(D).inverse() * L_s;

        return {L, L_s};
    }

    // TODO: implement cotan_weight
    // if (weight_type == WeightType::COTAN_WEIGHT)
    //{

    //}
}

// NOTE: b (upper part is laplacian coordinates, while the lower part is vertice coordinates)
Eigen::MatrixXd add_constraints(const Eigen::SparseMatrix<double> &L_s, const Eigen::MatrixXd &anchor_constraint, const Eigen::MatrixXd &handle_constraint)
{
    long n_vertice = L_s.rows();
    long n_anchor = anchor_constraint.rows();
    long n_handle = handle_constraint.rows();
    //spdlog::info("n_vertices: {}, n_anchor: {}, n_handle: {}", n_vertice, n_anchor, n_handle);
    Eigen::MatrixXd L_hat = Eigen::MatrixXd(L_s);
    L_hat.conservativeResize(n_vertice + n_anchor, L_hat.cols());
    L_hat.bottomRows(n_anchor) = anchor_constraint;
    
    L_hat.conservativeResize(L_hat.rows() + n_handle, L_hat.cols());
    L_hat.bottomRows(n_handle) = handle_constraint;


    spdlog::info("L_hat: {}", L_hat);
    //spdlog::info("anchor_constraint: {}", anchor_constraint);
    //L_hat.middleRows(n_vertice, n_vertice + n_anchor - 1) = anchor_constraint;
    //spdlog::info("L_hat: {}", L_hat.middleRows(n_vertice, n_vertice + n_anchor ));
    
    return L_hat;
};