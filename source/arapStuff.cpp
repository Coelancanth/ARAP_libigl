//
// Created by wilhelm on 1/21/22.
//

#include "arapStuff.h"
#include "iostream"
#include "math.h"
#include "MatrixLConstructor.h"
#include "MatrixBConstructor.h"
double getAngleBtwnVector(const Eigen::Vector3d & v1, const  Eigen::Vector3d & v2)
{
    double angle = v1.dot(v2) / v1.norm() / v2.norm();
    return acos(angle);
}



Eigen::Matrix3d rotationUpdateSingleVertex(int vertex,
                                           const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                           const std::vector<std::vector<int>> & adjList,
                                           const Eigen::MatrixXd & newPositions,
                                           WeightTable wt)

{
    // Note: Vertices positions newPositions are Nx3 matrices
    // Vertices are initial positions (step=0), not sure if duplicating positions though...

    const Eigen::MatrixXd& positions = vertices; // alias, in case the note above is incorrect...

    int vCount = vertices.rows();
    Eigen::Matrix3d Si;
    Eigen::Matrix3d Ri;

    Si.setZero();
    for (int vi: adjList[vertex]) {
        // get the edge vi-vertex
        Eigen::Matrix<double, 1, 3> eij = positions.row(vi) - positions.row(vertex);
        Eigen::Matrix<double, 1, 3> eijPrime = newPositions.row(vi) - newPositions.row(vertex);
        double wij = wt.getWeight(vi, vertex);
        Si += wij * eij.transpose() * eijPrime;
    }
    //std::cout << "vert" << vertex << "covmat\n" << Si << std::endl;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Si, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double det = svd.singularValues()(0) * svd.singularValues()(1) * svd.singularValues()(2);

    if (det > 0)
    {
        return svd.matrixV() * svd.matrixU().transpose();
    }
    if (det < 0)
    {
        Eigen::Matrix3d u = svd.matrixU();
        u.row(2) *= -1;
        return svd.matrixV() * svd.matrixU().transpose();
    }
    else
    {
        Eigen::Matrix3d res = svd.matrixV() * svd.matrixU().transpose();
        if (res.determinant() > 0) return res;
        else if (res.determinant() < 0)
        {
            Eigen::Matrix3d u = svd.matrixU();
            u.row(2) *= -1;
            return svd.matrixV() * svd.matrixU().transpose();
        }
        else
        {
            spdlog::warn("covariance matrix of vert {} is not full rank even after reconstruction! trying to return at least something...", vertex);
            return res;
        }

    }




}

std::vector<Eigen::Matrix3d> rotationUpdateStep(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
        const std::vector<std::vector<int>> & adjList,
        const Eigen::MatrixXd & newPositions,
        const WeightTable & wt
)
{
    int vCount = vertices.rows();
    std::vector<Eigen::Matrix3d> result;
    result.resize(vCount);

    for (int vi = 0; vi < vCount; vi++)
    {
        result[vi] = rotationUpdateSingleVertex(vi, vertices, faces, adjList, newPositions, wt);
    }
    return result;
}

Eigen::MatrixXd positionUpdateStep(const Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  & laplacianMatQR,
                                   const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions,
                                   const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                   const std::vector<std::vector<int>> &adjList,
                                   std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                                   //const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                                   WeightTable wt)
{

    // constructs the b matrix and solves for p' as in (9)
    int vertSize = vertices.rows();
    const Eigen::MatrixXd& positions = vertices;

    MatrixbConstructor bConstr(vertices, faces, adjList, rotMatrices,  wt, fixedPts, fixedPositions);

    // std::cout << "bMat: " <<  bConstr.bMat << std::endl;
    auto res = laplacianMatQR.solve(bConstr.bMat);
    // std::cout << "solved position~\n" << res << std::endl;

    return res.block(0,0,vertSize,3);
}

