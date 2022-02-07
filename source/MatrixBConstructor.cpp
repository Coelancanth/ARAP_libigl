//
// Created by wilhelm on 1/28/22.
//

#include "MatrixBConstructor.h"
#include "arapStuff.h"





MatrixbConstructor::MatrixbConstructor(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                                       const std::vector<std::vector<int>> &adjList,
                                       std::vector<Eigen::Matrix3d> rotMatrices, // contains Ri
                                       // const Eigen::MatrixXd& positions, // containes pi, shape is nx3
                                       WeightTable wt,
                                       const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions
)
{
    const Eigen::MatrixXd& positions = vertices;
    // std::cout << "ROT\n";
    for (int i = 0; i < rotMatrices.size(); i++)
    {

        // std::cout << "MAT" << i << "\n" <<  rotMatrices[i] << std::endl;
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
    // std::cout << "before adding fixed points, b mat:\n" << bMat << std::endl;

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
    // std::cout << "after adding fixed points, b mat:\n" << bMat << std::endl;

}