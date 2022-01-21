//
// Created by wilhelm on 1/21/22.
//


#ifndef ARAP_ARAPDEFORM_H
#define ARAP_ARAPDEFORM_H

#include "arapStuff.h"
#include "spdlog/spdlog.h"
class arapDeform {
public:
    arapDeform(const Eigen::MatrixXd& sourcevert,const Eigen::MatrixXi & sourcefaces;
    void setConstraints(const std::vector<int>& fixedPts, const std::vector<Eigen::Vector3d>&  fixedPositions);
    void compute( Eigen::MatrixXd& targetmesh, Eigen::MatrixXd initialGuess);
private:
    Eigen::ColPivHouseholderQR< Eigen::MatrixXd >  laplacianMatQR;
    std::vector<int> fixedPts;
    std::vector<Eigen::Vector3d>  fixedPositions;
    std::vector<std::vector<int>> adjList;
    const Eigen::MatrixXd& vertices;
    const Eigen::MatrixXi& faces;
    WeightTable wt;
};


#endif //ARAP_ARAPDEFORM_H
