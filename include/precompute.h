#ifndef PRECOMPUTE_H
#define PRECOMPUTE_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>


enum class WeightType
{
    UNIFORM_WEIGHT = 0,
    COTAN_WEIGHT = 1
};

// NOTE: should use template for handling different matrix types
/**
 * @brief Calculate symmetric laplacian matrix, L_s = DL 
 * D is the diagonal matrix such that D_ii = d_i, d_i is the number of immediate neighbors of i (the degree of valence of i)
 * L is laplacian matrix
 * L's entry:
 * - if i and j are neighbors, the ij_th entry is equal to the negative normalzied weight w_ij (w_ij divided by the sum of all weights from i)
 * - if i and j are not neighbors and i is not equal to j, then the ij_th entry is 0
 * - if i = j, then ij_th entry is 1
 * 
 * 
 * @param vertices 
 * @param faces 
 * @return a symmetric laplacian matrix L_s
 *
 */
Eigen::SparseMatrix<double> calculate_laplacian_matrix(
    const Eigen::MatrixXd &vertices, 
    const Eigen::MatrixXi &faces, 
    const WeightType &weight_type);

#endif //