//
// Created by wilhelm on 1/28/22.
//

#ifndef ARAP_MATRIXBCONSTRUCTOR_H
#define ARAP_MATRIXBCONSTRUCTOR_H

#include "WeightTable.h"

class MatrixbConstructor
{
// constructs the b matrix in (8)
// dimension of b is N*3, if I read correctly.
// rotMatrices is the vector with rotational matrices as per (6)
// then just access the bMat...
public:
    // WARNING: fixed points MUST be added in a consistent way OR slack variables will NOT work
    // TODO vertices and positions are one and the same
    MatrixbConstructor(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                       const std::vector<std::vector<int>> &adjList,
                       std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                       //const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                       WeightTable wt,
                       const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions
    );

    Eigen::MatrixXd bMat; // TODO make private

private:

    void setFixedPoints(const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions);
};
#endif //ARAP_MATRIXBCONSTRUCTOR_H
