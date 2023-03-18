#ifndef NODE_H
#define NODE_H

#include <QtGlobal>
#include <eigen3/Eigen/Geometry>

/**
* \brief A point in 3d space.
*/
class Node
{
public:
    Node(){
//        index = 0;
        coordinates = {0,0,0};
    }
    Node(Eigen::Vector3d coordinates){
        this->coordinates = coordinates;
    }

//    int index;
    /**
    * \brief The 3D coordinates of the node.
    */
    Eigen::Vector3d coordinates;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // NODE_H
