//
// Created by wilhelm on 2/7/22.
//

#include "arapDeformer.h"
#include "WeightTable.h"
#include "MatrixLConstructor.h"
#include "spdlog/stopwatch.h"
arapDeformer::arapDeformer(const Eigen::MatrixXi &sourcefaces, const Eigen::MatrixXd sourcevert)
    : faces(sourcefaces), vertices(sourcevert), wt(vertices, faces, adjList) {
    spdlog::info("creating adj matrix");
    igl::adjacency_list(faces, adjList, false);
    spdlog::info("creating weight table");
}


void arapDeformer::setVertices(const Eigen::MatrixXd &sourcevert) {
    vertices = sourcevert;
    //spdlog::info("creating weight table");
    //wt = WeightTable(vertices, faces, adjList);
}

void arapDeformer::setConstrainedPoints(const std::vector<int> &fixedPt) {
    fixedPts = fixedPt;
    spdlog::info("creating laplacian (sparse QR)");
    MatrixLConstructor lconstr (vertices, faces, adjList, wt, fixedPts);
    laplacianMatQR.compute(lconstr.laplacianMat);

}

void arapDeformer::updateConstraints(const std::vector<Eigen::Vector3d> &fixedPos) {
    fixedPositions = fixedPos;
}

void arapDeformer::compute(int iter) {
    int idxcount = vertices.rows();
    std::vector<Eigen::Matrix3d> rot;
    Eigen::MatrixXd pos;
    Eigen::MatrixXd& posprime = this->result; // this->result = posprime;
    spdlog::stopwatch sw0;
    rot = rotationUpdateStep(vertices, faces, adjList, vertices, wt);
    spdlog::info("rotation_update costs: {}", sw0);
    spdlog::stopwatch sw1;
    posprime = positionUpdateStep(laplacianMatQR, fixedPts, fixedPositions,
                                  vertices, faces, adjList, rot, wt);
    spdlog::info("position_update costs: {}", sw1);

    for (int i = 1; i < iter; i++)
    {
//        std::cout << "iter" << i << ":\n";
        rot = rotationUpdateStep(vertices, faces, adjList, posprime, wt);

        pos = posprime;
        posprime = positionUpdateStep(laplacianMatQR, fixedPts, fixedPositions,
                                      vertices, faces, adjList, rot, wt);
    }

}

