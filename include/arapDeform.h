//
// Created by wilhelm on 1/21/22.
//


#ifndef ARAP_ARAPDEFORM_H
#define ARAP_ARAPDEFORM_H

#include "arapStuff.h"
#include "spdlog/spdlog.h"
class arapDeform {
public:
    arapDeform(const Eigen::MatrixXd& sourcevert,const Eigen::MatrixXi & sourcefaces);
    void setConstraints(const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions);
    void updateConstraints(const std::vector<Eigen::Vector3d>&  fixedPositions);
    void compute(int iter);
    Eigen::MatrixXd getVertices() {return vertices;};
private:
    //void makeGuess(Eigen::MatrixXd & storage);
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   laplacianMatQR;
    std::vector<int> fixedPts;
    std::vector<Eigen::Vector3d>  fixedPositions;
    std::vector<std::vector<int>> adjList;
    Eigen::MatrixXd vertices; // vertices may change between passes...
    const Eigen::MatrixXi & faces; // but faces do not
    WeightTable wt;
};


#endif //ARAP_ARAPDEFORM_H
