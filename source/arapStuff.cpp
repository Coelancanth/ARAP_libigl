//
// Created by wilhelm on 1/21/22.
//

#include "arapStuff.h"
#include "iostream"
#include "math.h"
#include "MatrixLConstructor.h"
#include "MatrixBConstructor.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <Eigen/src/Core/Matrix.h>
#include <igl/adjacency_list.h>

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

    spdlog::stopwatch sw3;
    const Eigen::MatrixXd& positions = vertices; // alias, in case the note above is incorrect...

    int vCount = vertices.rows();
    Eigen::Matrix3d Si;
    // Eigen::Matrix3d Ri;
    Eigen::MatrixXd Pi;
    Eigen::MatrixXd PiPrime;
    Eigen::MatrixXd Di;
    Di.resize(adjList[vertex].size(), adjList[vertex].size());
    Pi.resize(adjList[vertex].size(),3);
    PiPrime.resize(adjList[vertex].size(),3);
    Di.setZero();
    PiPrime.setZero();
    Pi.setZero();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Pi); // for debugging, CLion gives an annoying error if I dont do this

    spdlog::stopwatch sw0;
    for (int i = 0; i < adjList[vertex].size(); i++) 
    {
        int vi = adjList[vertex][i];
        // get the edge vi-vertex
        Pi.row(i) = positions.row(vi) - positions.row(vertex);
        PiPrime.row(i) = newPositions.row(vi) - newPositions.row(vertex);
        double wij = wt.getWeight(vi, vertex);
        Di(i,i) = wij/2;
    }
    spdlog::info("constructing P_i, D_i, P_i' costs: {}", sw0);

    Si = Pi.transpose()*Di*PiPrime;
    //std::cout << "vert" << vertex << "covmat\n" << Si << std::endl;

    spdlog::stopwatch sw1;
    svd = Eigen::JacobiSVD<Eigen::MatrixXd>(Si, Eigen::ComputeThinU | Eigen::ComputeThinV); // for debugging
    spdlog::info("svd costs: {}", sw1);
    double det = svd.singularValues()(0) * svd.singularValues()(1) * svd.singularValues()(2);
    Eigen::Matrix3d r =  svd.matrixV() * svd.matrixU().transpose();
    spdlog::info("elapsed: {}", sw3);
    return r;

    if (det > 0)
    {
        spdlog::info("det > 0,");
    }
    if (det < 0)
    {
        spdlog::info("det < 0,");
        Eigen::Matrix3d u = svd.matrixU();
        u.row(2) *= -1;
        return svd.matrixV() * svd.matrixU().transpose();
    }
    else
    {
        spdlog::info("det = 0,");
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

//std::vector<Eigen::Matrix3d> rotationUpdateStep(
        //const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
        //const std::vector<std::vector<int>> & adjList,
        //const Eigen::MatrixXd & newPositions,
        //const WeightTable & wt
//)
//{
    //int vCount = vertices.rows();
    //std::vector<Eigen::Matrix3d> result;
    //result.resize(vCount);

    //for (int vi = 0; vi < vCount; vi++)
    //{
        //// NOTE:
        //spdlog::stopwatch sw2;
        //result[vi] = rotationUpdateSingleVertex(vi, vertices, faces, adjList, newPositions, wt);
        ////Eigen::Matrix3d m = rotationUpdateSingleVertex(vi, vertices, faces, adjList, newPositions, wt);
        ////spdlog::info("m for rotation_update_single_vertex costs: {}", sw2);
        ////result[vi] = m;
        //spdlog::info("rotation_update_single_vertex costs: {}", sw2);
        ////        if (vi < 5) std::cout << "V" << vi << " rotmat:\n" << result[vi] << std::endl;
    //}
    //return result;
//}

std::vector<Eigen::Matrix3d> rotationUpdateStep(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
        const std::vector<std::vector<int>> & adjList,
        const Eigen::MatrixXd & newPositions,
        WeightTable  wt
)
{
    int vCount = vertices.rows();
    std::vector<Eigen::Matrix3d> result;
    result.resize(vCount);
    for (int vi = 0; vi < vCount; vi++)
    {
        spdlog::stopwatch sw3;
        const Eigen::MatrixXd& positions = vertices; // alias, in case the note above is incorrect...

        int vCount = vertices.rows();
        Eigen::Matrix3d Si;
        // Eigen::Matrix3d Ri;
        Eigen::MatrixXd Pi;
        Eigen::MatrixXd PiPrime;
        Eigen::MatrixXd Di;
        Di.resize(adjList[vi].size(), adjList[vi].size());
        Pi.resize(adjList[vi].size(),3);
        PiPrime.resize(adjList[vi].size(),3);
        Di.setZero();
        PiPrime.setZero();
        Pi.setZero();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Pi); // for debugging, CLion gives an annoying error if I dont do this

        //spdlog::stopwatch sw0;
        //spdlog::info("vi: {}", vi);
        // FIXME: change vii to something else... naming is too ugly
        for (int i = 0; i < adjList[vi].size(); i++) 
        {
            //int vi = adjList[vi][i];
            //// get the edge vi-vi
            //Pi.row(i) = positions.row(vi) - positions.row(vi);
            //PiPrime.row(i) = newPositions.row(vi) - newPositions.row(vi);
            //double wij = wt.getWeight(vi, vi);
            //Di(i,i) = wij/2;
            int vii = adjList[vi][i];
            // get the edge vii-vertex
            Pi.row(i) = positions.row(vii) - positions.row(vi);
            PiPrime.row(i) = newPositions.row(vii) - newPositions.row(vi);
            double wij = wt.getWeight(vii,vi);
            Di(i,i) = wij/2;
        }
        //spdlog::info("constructing P_i, D_i, P_i' costs: {}", sw0);

        Si = Pi.transpose()*Di*PiPrime;
        //std::cout << "vert" << vertex << "covmat\n" << Si << std::endl;

        //spdlog::stopwatch sw1;
        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(Si, Eigen::ComputeThinU | Eigen::ComputeThinV); // for debugging
        //spdlog::info("svd costs: {}", sw1);
        double det = svd.singularValues()(0) * svd.singularValues()(1) * svd.singularValues()(2);
        Eigen::Matrix3d r =  svd.matrixV() * svd.matrixU().transpose();
        //spdlog::info("elapsed: {}", sw3);
        result[vi] = r;
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

    //spdlog::stopwatch sw;
    MatrixbConstructor bConstr(vertices, faces, adjList, rotMatrices,  wt, fixedPts, fixedPositions);
    //spdlog::info("constructing rhs costs: {}", sw);

//    std::cout << "bMat: " <<  bConstr.bMat.block(0,0,5,3) << std::endl;
//    std::cout << "bMat: " <<  bConstr.bMat << std::endl;
    //spdlog::stopwatch sw1;
    auto res = laplacianMatQR.solve(bConstr.bMat);
    //spdlog::info("solving costs: {}", sw1);
    // std::cout << "solved position~\n" << res << std::endl;
//    std::cout << "Head of vertexPrime\n" << res.block(0,0,5,3) << "\n";
    return res.block(0,0,vertSize,3);
}

