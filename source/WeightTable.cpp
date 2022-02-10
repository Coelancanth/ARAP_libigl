//
// Created by wilhelm on 1/28/22.
//
#include "arapStuff.h"
#include "WeightTable.h"

WeightTable::WeightTable()
{
    useCotanWeight = false;
    spdlog::warn("using uniform weight. this may be problematic");
}
double WeightTable::getWeight(int i, int j)
{
    if (!useCotanWeight) return 1.0;
    return weights(i,j);
}
WeightTable::WeightTable(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                         const std::vector<std::vector<int>> & adjList)
{
    int vertCount = vertices.rows();
    weights.resize(vertCount, vertCount);
    weights.setZero();
    // calculate the cotan weight matrix
    // the cotan weight is calculated as follows ...
    // for each edge, find the 2 angles that are facing it. then substitute into the w_ij equation

    // to save some time, I use this approach:
    // for each triangle
    //  for each edge in triangle:
    //   let i,j be its 2 vertices
    //     weight[i,j] and weight[j,i] += 0.5cot(angle opposite of edge)
    useCotanWeight = true;
    for (int ti = 0; ti < faces.rows(); ti++)
    {
        int v1 = faces(ti,0);
        int v2 = faces(ti,1);
        int v3 = faces(ti,2);
        const auto& v1p = vertices.row(v1);
        const auto& v2p = vertices.row(v2);
        const auto& v3p = vertices.row(v3);
        double ang1 = getAngleBtwnVector(v3p-v1p, v2p-v1p);
        double ang2 = getAngleBtwnVector(v3p-v2p, v1p-v2p);
        double ang3 = getAngleBtwnVector(v1p-v3p, v2p-v3p);
        double cot1 = 1 / tan(ang1);
        double cot2 = 1 / tan(ang2);
        double cot3 = 1 / tan(ang3);

        weights(v1,v2) += cot3/2;
        weights(v2,v1) += cot3/2;
        weights(v1,v3) += cot2/2;
        weights(v3,v1) += cot2/2;
        weights(v3,v2) += cot1/2;
        weights(v2,v3) += cot1/2;
    }
    // std::cout << "weight init complete.\n" << weights << std::endl;

}

