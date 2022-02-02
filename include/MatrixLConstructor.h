//
// Created by wilhelm on 1/28/22.
//

#ifndef ARAP_MATRIXLCONSTRUCTOR_H
#define ARAP_MATRIXLCONSTRUCTOR_H
#include "WeightTable.h"
#include "Eigen/SparseCore"
#include <Eigen/SparseQR>

class MatrixLConstructor
{
// construct the Laplace - Beltrami operator according to (8)
// usage: call constructor and setFixedPoints function, and
// just get invLapMat when done, you may need to invert it.
// should only be constructed once, btw.
// adjList is the adjacency list that you should calculate through
// igl::adjacency_list(*faces, adjList, false); first
// this is not calculated in the function as it is needed everywhere
public:
    MatrixLConstructor(
            const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
            const std::vector<std::vector<int>> & adjList,
            WeightTable wt, const std::vector<int> & fixedPts)
    ;



    // function removed as this is probably too messy to compute in the class.
//    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > getQR();


    Eigen::SparseMatrix<double> laplacianMat; // laplacian beltrami operator in (9) // TODO make private

private:
    // WARNING: fixed points MUST be added in a consistent way OR slack variables will NOT work
    // fixed points order MUST be consistent across L, b
    void setFixedPoints(const std::vector<int> & fixedPts);
};


#endif //ARAP_MATRIXLCONSTRUCTOR_H
