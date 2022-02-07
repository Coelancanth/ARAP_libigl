//
// Created by wilhelm on 1/28/22.
//

#include "MatrixLConstructor.h"



#include "arapStuff.h"
#include "Eigen/SparseCore"
#include "Eigen/SparseQR"

MatrixLConstructor::MatrixLConstructor(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
        const std::vector<std::vector<int>> & adjList,
        WeightTable wt , const std::vector<int> & fixedPts)

{
    // construct the matrix according to (8)
    // TODO: make matrix sparse
    int vCount = vertices.rows();
    int cCount = fixedPts.size();




    // it is expected that each row has 5 nonzero entries -> 3 neighbors and 1 on diagonal and 1 optional slack variable.
    // the last few rows should have 1 entry each
    laplacianMat.resize(vCount+cCount, vCount+cCount);        // default is column major
    laplacianMat.reserve(Eigen::VectorXi::Constant(vCount+cCount,5));
    for (int vi = 0; vi < vCount; vi++)
    {
        // construct the sparse matrix, one row at a time.
        const std::vector<int> neighbors = adjList[vi];
        for (int nei: neighbors)
        {
            laplacianMat.coeffRef(vi, nei) -= wt.getWeight(vi, nei);
            laplacianMat.coeffRef(vi, vi) += wt.getWeight(vi, nei);
        }
    }
    // std::cout <<"Lap Matrix: before constr" << std::endl << laplacianMat << std::endl;

    setFixedPoints(fixedPts);
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> detcalc;
    detcalc.compute(laplacianMat);

//    std::cout <<" determinant of L matrix " << detcalc.determinant() << " or " << detcalc.logAbsDeterminant() << std::endl;

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
    // seems setting elements to zero in Eigen::SparseMatrix
    // is sufficient and the makecompressed will remove the zeroes entries.
    // however we dont know anymore which columns should be removed so...




    for (int i = 0; i < constrLen; i++)
    {
        int point = fixedPts[i];

        laplacianMat.coeffRef(i+vertexLen, point) += 1;
        laplacianMat.coeffRef(point, i+vertexLen) += 1;


    }
    laplacianMat.makeCompressed();

//    std::cout << "Head of Laplacian\n";
//    for ( int i = 0; i < 5; i++)
//    {
//        std::cout << laplacianMat.row(i) << std::endl;
//    }


}

//Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > MatrixLConstructor::getQR()
//{
//    // std::cout <<"Lap Matrix: \n" << laplacianMat << std::endl;
//    //return laplacianMat.colPivHouseholderQr();
//    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > qrsolver;
//    qrsolver.compute(laplacianMat);
//    return qrsolver;
//}
