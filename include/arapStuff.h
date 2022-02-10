//
// Created by wilhelm on 1/21/22.
//

#ifndef ARAP_ARAPSTUFF_H
#define ARAP_ARAPSTUFF_H

#include <memory>
#include "Eigen/Core"
#include "assert.h"
#include "vector"
#include "igl/adjacency_list.h"
#include "exception"
#include "spdlog/spdlog.h"
#include <iostream>
#include "WeightTable.h"

double getAngleBtwnVector(const Eigen::Vector3d & v1, const  Eigen::Vector3d &v2);

Eigen::Matrix3d rotationUpdateSingleVertex(int vertex,
                                           const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                           //const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions,
                                           const std::vector<std::vector<int>> & adjList,
                                           const Eigen::MatrixXd & newPositions,
                                           WeightTable wt);


std::vector<Eigen::Matrix3d> rotationUpdateStep(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& faces,
        //const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions,
        const std::vector<std::vector<int>> & adjList,
        const Eigen::MatrixXd& newPositions,

        const WeightTable&  wt);

Eigen::MatrixXd positionUpdateStep(const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>   & laplacianMatQR,
                                   const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions,
                                   const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                   const std::vector<std::vector<int>> &adjList,
                                   std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                                   //const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                                   WeightTable wt);



#endif //ARAP_ARAPSTUFF_H
