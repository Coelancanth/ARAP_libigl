//
// Created by wilhelm on 1/21/22.
//

#include "arapDeform.h"
#include "igl/adjacency_list.h"
arapDeform::arapDeform(const Eigen::MatrixXd &sourcevert, const Eigen::MatrixXi &sourcefaces):
    adjList(), faces(sourcefaces), vertices(sourcevert)


{
    spdlog::info("creating adj matrix, using libigl functions");
    igl::adjacency_list(faces, adjList, false);
    spdlog::info("creating weight table");
    wt = WeightTable(vertices, faces, adjList);

}

void arapDeform::setConstraints(const std::vector<int> &fixedPt, const std::vector<Eigen::Vector3d> &fixedPos) {
    fixedPositions = fixedPos;
    fixedPts = fixedPt;

    spdlog::info("creating laplacian (QR)");

    MatrixLConstructor lconstr (vertices, faces, adjList, wt, fixedPts);
    // lconstr.setFixedPoints(fixedPts);
    laplacianMatQR = lconstr.getQR();

    // replace the corresponding vertices with the fixed value
//    for (int i=0; i < fixedPos.size(); i++)
//    {
//        int pt = fixedPt[i];
//        this->vertices.row(pt) = fixedPos[i];
//    }

}

Eigen::MatrixXd arapDeform::compute(const Eigen::MatrixXd & initialGuess) {
    // how to guess though?
    std::vector<Eigen::Matrix3d> rot = rotationUpdateStep(
            vertices, faces, adjList, initialGuess, wt
            );

    // use all identity for first iteration?
//    int idxcount = initialGuess.rows();
//    std::vector<Eigen::Matrix3d> rot (idxcount, Eigen::Matrix3d::Identity());
    Eigen::MatrixXd pos = vertices;
    Eigen::MatrixXd posprime = positionUpdateStep (laplacianMatQR, fixedPts, fixedPositions,
                                                   vertices, faces,
                                                   adjList, rot, initialGuess, wt);;
    this->vertices = posprime;
    return posprime;
}