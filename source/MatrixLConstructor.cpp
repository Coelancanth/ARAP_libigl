//
// Created by wilhelm on 1/28/22.
//

#include "MatrixLConstructor.h"



#include "arapStuff.h"

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

    // std::cout <<"Lap Matrix: before constr" << std::endl << laplacianMat << std::endl;

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
    // std::cout <<"Lap Matrix: \n" << laplacianMat << std::endl;
    return laplacianMat.colPivHouseholderQr();
}
