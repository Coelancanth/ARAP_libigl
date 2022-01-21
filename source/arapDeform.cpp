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
    wt = WeightTable();

}

void arapDeform::setConstraints(const std::vector<int> &fixedPt, const std::vector<Eigen::Vector3d> &fixedPos) {
    fixedPositions = fixedPos;
    fixedPts = fixedPt;

    spdlog::info("creating laplacian (QR)");

    MatrixLConstructor lconstr (vertices, faces, adjList, wt);
    laplacianMatQR = lconstr.getQR();

}

Eigen::MatrixXd arapDeform::compute(Eigen::MatrixXd &targetmesh, Eigen::MatrixXd initialGuess) {
    // how to guess though?
    std::vector<Eigen::Matrix3d> rot = rotationUpdateStep(
            vertices, faces, adjList, initialGuess, wt
            );
    Eigen::MatrixXd pos = positionUpdateStep (laplacianMatQR, fixedPts, fixedPositions,
                              vertices, faces,
                              adjList, rot, initialGuess, wt);
    bool converged = false;
    int stepcount = 0;
    while (!converged)
    {
        rot = rotationUpdateStep(vertices, faces, adjList, pos, wt);
        pos = positionUpdateStep (laplacianMatQR, fixedPts, fixedPositions,
                                 vertices, faces,
                                 adjList, rot, pos, wt);
        stepcount += 1;

        if (stepcount >= 50)
        {
            converged = true;
            // TODO: stopping condition?
        }
    }
    return pos;
}