#ifndef CUBOID_H
#define CUBOID_H

#include <eigen3/Eigen/Geometry>

static constexpr float inf = std::numeric_limits<float>::infinity();

/**
* \class Cuboid
* \brief Simple 3d axis aligned cuboid class.
*/
class Cuboid
{
public:

    Cuboid(){
        minPoint = {inf,inf,inf};
        maxPoint = {-inf,-inf,-inf};
    }

    Cuboid(Eigen::Vector3d minPoint, Eigen::Vector3d maxPoint){
        this->minPoint = minPoint;
        this->maxPoint = maxPoint;
    }

    /**
    * \brief The cuboid corner with the largest x,y and z values.
    */
    Eigen::Vector3d minPoint;

    /**
    * \brief The cuboid corner with the smallest x,y and z values.
    */
    Eigen::Vector3d maxPoint;

    inline Cuboid & operator |= (const Eigen::Vector3d & p)
    {
        minPoint(0) = std::min(minPoint(0), p(0));
        minPoint(1) = std::min(minPoint(1), p(1));
        minPoint(2) = std::min(minPoint(2), p(2));
        maxPoint(0) = std::max(maxPoint(0), p(0));
        maxPoint(1) = std::max(maxPoint(1), p(1));
        maxPoint(2) = std::max(maxPoint(2), p(2));
        return *this;
    }
};

#endif // CUBOID_H
