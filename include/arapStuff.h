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

class WeightTable
{
public:
    WeightTable ();
    WeightTable(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                       const std::vector<std::vector<int>> & adjList);

    double getWeight(int i, int j);


private:
    bool useCotanWeight;
    Eigen::MatrixXd weights;
};


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




    Eigen::ColPivHouseholderQR< Eigen::MatrixXd > getQR();


    Eigen::MatrixXd laplacianMat; // laplacian beltrami operator in (9) // TODO make private

private:
    // WARNING: fixed points MUST be added in a consistent way OR slack variables will NOT work
    // fixed points order MUST be consistent across L, b
    void setFixedPoints(const std::vector<int> & fixedPts);
};


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
                       const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                       WeightTable wt,
                       const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions
                       );

    Eigen::MatrixXd bMat; // TODO make private

private:

    void setFixedPoints(const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions);
};

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

Eigen::MatrixXd positionUpdateStep(const Eigen::ColPivHouseholderQR< Eigen::MatrixXd > & laplacianMatQR,
                                   const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions,
                                   const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                   const std::vector<std::vector<int>> &adjList,
                                   std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                                   const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                                   WeightTable wt);



#endif //ARAP_ARAPSTUFF_H