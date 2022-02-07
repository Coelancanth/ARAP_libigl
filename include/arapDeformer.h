//
// Created by wilhelm on 2/7/22.
//

#ifndef ARAP_ARAPDEFORMER_H
#define ARAP_ARAPDEFORMER_H

#include "arapStuff.h"
#include "spdlog/spdlog.h"
class arapDeformer {
public:
    arapDeformer(const Eigen::MatrixXi & sourcefaces, const Eigen::MatrixXd sourcevert) ;

    // this initializes the L matrix and the QR structure only
    // you have to provide the coor of the constrained points when (re)using the L matrix
    void setConstrainedPoints(const std::vector<int>& fixedPts);

    // only put in vertices of the same mesh!
    void setVertices(const Eigen::MatrixXd& sourcevert);

    void updateConstraints(const std::vector<Eigen::Vector3d>&  fixedPositions);
    // run for a few iterations
    void compute(int iter);

    Eigen::MatrixXd getVertices() {return result;};
private:

    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   laplacianMatQR;
    std::vector<int> fixedPts;
    std::vector<Eigen::Vector3d>  fixedPositions;
    std::vector<std::vector<int>> adjList;

    Eigen::MatrixXi faces; // topology of mesh
    Eigen::MatrixXd vertices;
    Eigen::MatrixXd result;

    WeightTable wt;
};


#endif //ARAP_ARAPDEFORMER_H
