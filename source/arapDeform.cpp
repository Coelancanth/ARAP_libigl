//
// Created by wilhelm on 1/21/22.
//

#include "arapDeform.h"
#include "igl/adjacency_list.h"
#include "MatrixBConstructor.h"
#include "MatrixLConstructor.h"
#include "arapStuff.h"
arapDeform::arapDeform(const Eigen::MatrixXd &sourcevert, const Eigen::MatrixXi &sourcefaces):
    adjList(), faces(sourcefaces), vertices(sourcevert)


{
    spdlog::info("creating adj matrix, using libigl functions");
    igl::adjacency_list(faces, adjList, false);
    spdlog::info("creating weight table");
    wt = WeightTable(vertices, faces, adjList);

}
void arapDeform::updateConstraints(const std::vector<Eigen::Vector3d> &fixedPos)
{
    fixedPositions = fixedPos;
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

void arapDeform::compute(int iter) {
    // how to guess though?
//    std::vector<Eigen::Matrix3d> rot = rotationUpdateStep(
//            vertices, faces, adjList, initialGuess, wt
//            );

    // use all identity for first iteration?
    int idxcount = this->vertices.rows();
    std::vector<Eigen::Matrix3d> rot (idxcount, Eigen::Matrix3d::Identity());

    Eigen::MatrixXd pos = vertices;
    Eigen::MatrixXd posprime = positionUpdateStep (laplacianMatQR, fixedPts, fixedPositions,
                                                   vertices, faces,
                                                   adjList, rot,  wt);

    for (int i = 1; i < iter; i++)
    {
        rot = rotationUpdateStep(pos, faces, adjList, posprime, wt);
        pos = posprime;
        posprime = positionUpdateStep(laplacianMatQR, fixedPts, fixedPositions,
                                      pos, faces, adjList, rot, wt);
    }
    this->vertices = posprime;

}