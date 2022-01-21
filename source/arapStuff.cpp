//
// Created by wilhelm on 1/21/22.
//

#include "arapStuff.h"
#include "iostream"

MatrixLConstructor::MatrixLConstructor(
                                       const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                       const std::vector<std::vector<int>> & adjList,
                                       WeightTable wt)

    {
        // construct the matrix according to (8)
        // TODO: make matrix sparse
        int vCount = vertices.rows();

        laplacianMat.resize(vCount, vCount);
        laplacianMat.setZero();

        for (int vi = 0; vi < vCount; vi++)
        {
            // construct the sparse matrix, one row at a time.
            const std::vector<int> neighbors = adjList[vi];
            for (int nei: neighbors)
            {


                laplacianMat(vi, nei) = -wt.getWeight(vi, nei);
                laplacianMat(vi, vi) += wt.getWeight(vi, nei);

            }
        }

        std::cout <<"Lap Matrix: before constr" << std::endl << laplacianMat << std::endl;


    };

void MatrixLConstructor::setFixedPoints(const std::vector<int> & fixedPts)
{
    // overwrite the corresponding places in the matrix
    for (int vi: fixedPts)
    {
        laplacianMat.row(vi).setZero();
        laplacianMat(vi, vi) = 1;
    }
    std::cout <<"Lap Matrix: after constr" << laplacianMat << std::endl;

}

Eigen::ColPivHouseholderQR< Eigen::MatrixXd > MatrixLConstructor::getQR()
{
    std::cout <<"Lap Matrix: {}" << laplacianMat << std::endl;
    return laplacianMat.colPivHouseholderQr();
}


MatrixbConstructor::MatrixbConstructor(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                   const std::vector<std::vector<int>> &adjList,
                   std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                   const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                   WeightTable wt)
{
    int vCount = vertices.rows();
    bMat.resize(vCount, 3);
    bMat.setZero();
    for (int vi = 0; vi < vCount; vi++)
    {
        // first compute the vector, then transpose and put into b.
        const std::vector<int> neighbors = adjList[vi];
        for (int nei: neighbors)
        {
            double wij = wt.getWeight(vi, nei);
            Eigen::Matrix3d RiRjT = (rotMatrices[vi] + rotMatrices[nei]).transpose();
            Eigen::Matrix<double, 1,3> pipjT = (positions.row(vi) - positions.row(nei)); // should make a row vector of 1x3...

            //spdlog::info("bMat row nei: {}x{}", bMat.row(nei).rows(),bMat.col(nei).cols() );
            Eigen::Matrix<double, 1, 3> rhs = wij / 2 * (pipjT * RiRjT);
            //spdlog::info("bMat rhs: {}x{}", rhs.rows(),rhs.cols() );


            bMat.row(nei) += rhs;


        }
    }
}

void MatrixbConstructor::setFixedPoints(const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions)
{
    // overwrite the corresponding places in the matrix with the fixed positions
    for (int vi = 0; vi < fixedPts.size() ; vi++)
    {
        // TODO: this following line does not work, read Eigen::Map to rewrite this line
        // bMat.row(vi) = fixedPositions[vi].resize(1,3);

        bMat(vi,0) = fixedPositions[vi](0);
        bMat(vi,1) = fixedPositions[vi](1);
        bMat(vi,2) = fixedPositions[vi](2);

    }
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
    for (int vi: adjList[vertex])
    {
        // get the edge vi-vertex
        Eigen::Matrix<double, 1, 3> eij = positions.row(vi) - positions.row(vertex);
        Eigen::Matrix<double, 1, 3> eijPrime = newPositions.row(vi) - newPositions.row(vertex);
        double wij = wt.getWeight(vi, vertex);
        Si += wij * eij.transpose() * eijPrime;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Si, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double det = svd.singularValues()(0) * svd.singularValues()(1) * svd.singularValues()(2);
    if (det > 0)
    {
        return svd.matrixV() * svd.matrixU().transpose();
    }
    else
    {
        spdlog::info("vertex {} svd negative det = {}", vertex, det);
        Eigen::Matrix3d u = svd.matrixU();
        u.row(2) *= -1;
        return svd.matrixV() * u.transpose();
    }
}

std::vector<Eigen::Matrix3d> rotationUpdateStep(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
        const std::vector<std::vector<int>> & adjList,
        const Eigen::MatrixXd & newPositions,
        WeightTable wt
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

Eigen::MatrixXd positionUpdateStep(const Eigen::ColPivHouseholderQR< Eigen::MatrixXd > & laplacianMatQR,
                                   const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions,
                                   const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                   const std::vector<std::vector<int>> &adjList,
                                   std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                                   const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                                   WeightTable wt)
{

    // constructs the b matrix and solves for p' as in (9)
    MatrixbConstructor bConstr(vertices, faces, adjList, rotMatrices, positions, wt);
    bConstr.setFixedPoints(fixedPts, fixedPositions);
    std::cout << "bMat: " <<  bConstr.bMat << std::endl;
    return laplacianMatQR.solve(bConstr.bMat);

}

