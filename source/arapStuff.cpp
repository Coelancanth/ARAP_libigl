//
// Created by wilhelm on 1/21/22.
//

#include "arapStuff.h"
#include "iostream"
#include "math.h"
double getAngleBtwnVector(const Eigen::Vector3d & v1, const  Eigen::Vector3d & v2)
{
    double angle = v1.dot(v2) / v1.norm() / v2.norm();
    return acos(angle);
}

WeightTable::WeightTable()
{
    useCotanWeight = false;
}
double WeightTable::getWeight(int i, int j)
{
    if (!useCotanWeight) return 1.0;
    return weights(i,j);
}
WeightTable::WeightTable(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
            const std::vector<std::vector<int>> & adjList)
{
    int vertCount = vertices.rows();
    weights.resize(vertCount, vertCount);
    weights.setZero();
    // calculate the cotan weight matrix
    // the cotan weight is calculated as follows ...
    // for each edge, find the 2 angles that are facing it. then substitute into the w_ij equation

    // to save some time, I use this approach:
    // for each triangle
    //  for each edge in triganle:
    //   let i,j be its 2 vertices
    //     weight[i,j] and weight[j,i] += 0.5cot(angle opposite of edge)
    useCotanWeight = true;
    for (int ti = 0; ti < faces.rows(); ti++)
    {
        int v1 = faces(ti,0);
        int v2 = faces(ti,1);
        int v3 = faces(ti,2);
        const auto& v1p = vertices.row(v1);
        const auto& v2p = vertices.row(v2);
        const auto& v3p = vertices.row(v3);
        double ang1 = getAngleBtwnVector(v3p-v1p, v2p-v1p);
        double ang2 = getAngleBtwnVector(v3p-v2p, v1p-v2p);
        double ang3 = getAngleBtwnVector(v1p-v3p, v2p-v3p);
        double cot1 = 1 / tan(ang1);
        double cot2 = 1 / tan(ang2);
        double cot3 = 1 / tan(ang3);

        weights(v1,v2) += cot3/2;
        weights(v2,v1) += cot3/2;
        weights(v1,v3) += cot2/2;
        weights(v3,v1) += cot2/2;
        weights(v3,v2) += cot1/2;
        weights(v2,v3) += cot1/2;
    }
    std::cout << "weight init complete.\n" << weights << std::endl;

}







MatrixLConstructor::MatrixLConstructor(
                                       const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                       const std::vector<std::vector<int>> & adjList,
                                       WeightTable wt , const std::vector<int> & fixedPts)

    {
        // construct the matrix according to (8)
        // TODO: make matrix sparse
        int vCount = vertices.rows();
        int cCount = fixedPts.size();

        laplacianMat.resize(vCount+cCount, vCount+cCount);
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

        setFixedPoints(fixedPts);

    };

void MatrixLConstructor::setFixedPoints(const std::vector<int> & fixedPts)
{
    // overwrite the corresponding places in the matrix, according to paper..
//    for (int vi: fixedPts)
//    {
//
//        laplacianMat.row(vi).setZero();
//        //laplacianMat.col(vi).setZero();
//        laplacianMat(vi, vi) = 1;
//    }
//    std::cout <<"Lap Matrix: after fixed ptr constr\n" << laplacianMat << std::endl;

    // add "slack" variables..
    int constrLen = fixedPts.size();
    int vertexLen = this->laplacianMat.rows() - constrLen;
    for (int i = 0; i < constrLen; i++)
    {
        int point = fixedPts[i];
        laplacianMat.col(i  + vertexLen).setZero();
        laplacianMat.row(i  + vertexLen).setZero();

        laplacianMat(i+vertexLen, point) = 1;
        laplacianMat(point, i+vertexLen) = 1;


    }

}

Eigen::ColPivHouseholderQR< Eigen::MatrixXd > MatrixLConstructor::getQR()
{
    std::cout <<"Lap Matrix: \n" << laplacianMat << std::endl;
    return laplacianMat.colPivHouseholderQr();
}


MatrixbConstructor::MatrixbConstructor(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                   const std::vector<std::vector<int>> &adjList,
                   std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                   const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                   WeightTable wt,
                   const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions
                   )
{
    std::cout << "ROT\n";
    for (int i = 0; i < rotMatrices.size(); i++)
    {

        std::cout << "MAT" << i << "\n" <<  rotMatrices[i] << std::endl;
    }

    int vCount = vertices.rows();
    int cCount = fixedPts.size();
    bMat.resize(vCount+cCount, 3);
    bMat.setZero();
    for (int vi = 0; vi < vCount; vi++)
    {
        // first compute the vector, then transpose and put into b.
        const std::vector<int> neighbors = adjList[vi];
        for (int nei: neighbors)
        {
            double wij = wt.getWeight(vi, nei);
            Eigen::Matrix3d RiRj = (rotMatrices[vi] + rotMatrices[nei]);
            Eigen::Matrix<double, 3,1> pipj = (positions.row(vi) - positions.row(nei)).transpose(); // should make a row vector of 1x3...

            //spdlog::info("bMat row nei: {}x{}", bMat.row(nei).rows(),bMat.col(nei).cols() );
            Eigen::Matrix<double, 3,1> rhsT = wij / 2 * (RiRj)*pipj;
            //spdlog::info("bMat rhs: {}x{}", rhs.rows(),rhs.cols() );


            bMat.row(vi) += rhsT.transpose();


        }
    }
    this->setFixedPoints(fixedPts, fixedPositions);
}

void MatrixbConstructor::setFixedPoints(const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions)
{

    // overwrite the corresponding places in the matrix with the fixed positions
    std::cout << "before adding fixed points, b mat:\n" << bMat << std::endl;

    int constrCount = fixedPts.size();
    int idxCount = this->bMat.rows() - constrCount;


//  old version ...
//    for (int vi = 0; vi < fixedPts.size() ; vi++)
//    {
//        int idx = fixedPts[vi];
//        bMat.row(idx) = fixedPositions[vi];
//    }
//   new version ...
    for (int i = 0; i < constrCount; i++)
    {
        bMat.row(i + idxCount) = fixedPositions[i];
    }
    std::cout << "after adding fixed points, b mat:\n" << bMat << std::endl;

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
    std::cout << "vert" << vertex << "covmat\n" << Si << std::endl;
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

Eigen::MatrixXd positionUpdateStep(const Eigen::ColPivHouseholderQR< Eigen::MatrixXd > & laplacianMatQR,
                                   const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions,
                                   const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                   const std::vector<std::vector<int>> &adjList,
                                   std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                                   const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                                   WeightTable wt)
{

    // constructs the b matrix and solves for p' as in (9)
    int vertSize = vertices.rows();

    MatrixbConstructor bConstr(vertices, faces, adjList, rotMatrices, positions, wt, fixedPts, fixedPositions);

    std::cout << "bMat: " <<  bConstr.bMat << std::endl;
    auto res = laplacianMatQR.solve(bConstr.bMat);
    std::cout << "solved position~\n" << res << std::endl;

    return res.block(0,0,vertSize,3);
}

