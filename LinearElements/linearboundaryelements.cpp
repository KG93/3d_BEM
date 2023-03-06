#include "linearboundaryelements.h"

LinearBoundaryElements::LinearBoundaryElements()
{

}

void LinearBoundaryElements::calculateAssociatedTriangles()
{
    long numberOfNodes = nodes.size();
    long numberOfTriangles = triangles.size();

    for( int i=0; i<numberOfTriangles; i++)
    {
        LinearTriangle tmpTriangle = triangles.at(i);
        if(tmpTriangle.node1Index > 0 && tmpTriangle.node1Index < numberOfNodes)
        {
            nodes[tmpTriangle.node1Index].associatedTriangles.push_back(i);
        }
        if(tmpTriangle.node2Index > 0 && tmpTriangle.node2Index < numberOfNodes)
        {
            nodes[tmpTriangle.node2Index].associatedTriangles.push_back(i);
        }
        if(tmpTriangle.node3Index > 0 && tmpTriangle.node3Index < numberOfNodes)
        {
            nodes[tmpTriangle.node3Index].associatedTriangles.push_back(i);
        }
    }
}

void LinearBoundaryElements::calculateCollocationPoints()
{
    collocationPoints.resize(triangles.size());
    for(int i=0;i<collocationPoints.size();i++)
    {
        LinearTriangle tmpLinTriangle = triangles.at(i);
        Eigen::Vector3d NodeCoord1 = nodes.at(tmpLinTriangle.node1Index).coordinates;
        Eigen::Vector3d NodeCoord2 = nodes.at(tmpLinTriangle.node2Index).coordinates;
        Eigen::Vector3d NodeCoord3 = nodes.at(tmpLinTriangle.node3Index).coordinates;
        VectorTriangle tmpTriangle = VectorTriangle(NodeCoord1, NodeCoord2, NodeCoord3);
        collocationPoints[i]=global::midpointOfTriangle(tmpTriangle);
    }
}

void LinearBoundaryElements::calculateTrianglesArea()
{
    long numberOfTriangles = triangles.size();
    trianglesArea.resize(numberOfTriangles);

    for(int i=0; i < trianglesArea.size(); i++)
    {
        LinearTriangle tmpLinTriangle = triangles.at(i);
        Eigen::Vector3d NodeCoord1 = nodes.at(tmpLinTriangle.node1Index).coordinates;
        Eigen::Vector3d NodeCoord2 = nodes.at(tmpLinTriangle.node2Index).coordinates;
        Eigen::Vector3d NodeCoord3 = nodes.at(tmpLinTriangle.node3Index).coordinates;
        VectorTriangle tmpTriangle = VectorTriangle(NodeCoord1, NodeCoord2, NodeCoord3);
        trianglesArea(i) = global::areaOfTriangle(tmpTriangle);
    }
}

void LinearBoundaryElements::calculateTrianglesAreaAndAverageDimension()
{
    calculateTrianglesArea();
    if(trianglesArea.size() > 0)
    {
        averageTriangleDim = trianglesArea.array().sqrt().mean() / 2;
    }
    else
    {
        averageTriangleDim = 0;
    }
}
