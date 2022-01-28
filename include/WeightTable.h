//
// Created by wilhelm on 1/28/22.
//

#ifndef ARAP_WEIGHTTABLE_H
#define ARAP_WEIGHTTABLE_H
#include <memory>
#include "Eigen/Core"
#include "assert.h"
#include "vector"
#include "igl/adjacency_list.h"
#include "exception"
#include "spdlog/spdlog.h"
class WeightTable
{
public:
    WeightTable ();
    WeightTable(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi&faces,
                const std::vector<std::vector<int>> & adjList);

    double getWeight(int i, int j);


private:
    bool useCotanWeight;
    Eigen::MatrixXd weights;
};


#endif //ARAP_WEIGHTTABLE_H
