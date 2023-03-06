#include "boundaryelements.h"
#include "qmatrix4x4.h"

BoundaryElements::BoundaryElements()
{

}

void BoundaryElements::calculateCollocationPoints()
{
    collocationPoints.resize(triangles.size());
    for(int i=0;i<collocationPoints.size();i++)
    {
        collocationPoints[i]=global::midpointOfTriangle(triangles.at(i));
    }
}

void BoundaryElements::calculateTriangleMidPoints()
{
    for(int i=0; i<triangles.size(); i++)
    {
        triangles[i].triangleMidpoint = global::midpointOfTriangle(triangles.at(i));
    }
}

void BoundaryElements::calculateTrianglesArea()
{
    trianglesArea.resize(triangles.size());
    for(int i=0; i < trianglesArea.size(); i++)
    {
        trianglesArea(i) = global::areaOfTriangle(triangles.at(i));
    }
}

void BoundaryElements::calculateTrianglesAreaAndAverageDimension()
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

//void BoundaryElements::setUpBoundaryConditions()
//{
//    robinBoundaryConditions.resize(triangles.size());
//    for(int i=0;i<robinBoundaryConditions.size();i++)
//    {
//        robinBoundaryConditions[i]=triangles.at(i).robinBoundaryCondition;
//    }
//}

void BoundaryElements::flipNormals()
{
    for(int i=0;i<triangles.size();i++)
    {
        triangles[i].normal = -triangles.at(i).normal;
    }
}

void BoundaryElements::swapElementIndexes(int index1, int index2)
{
    int numberOfTriangles = triangles.length();
    int triangleAreaVecLength = trianglesArea.size();
    int collocationPointsVecLength = collocationPoints.length();
    int triangleColorsVecLength = triangleColors.length();
    int phiSolutionVecLength = phiSolution.size();
    int dPhiSolutionVecLength = dPhiSolution.size();
    int soundPressureVecLength = soundPressure.size();

    bool swapTriangleColors = true;
    bool swapSolutions = true;
    if(index1 >= numberOfTriangles || index2 >= numberOfTriangles || index1 < 0 || index2 < 0)     ///check what to swap
    {
        std::cout<<"Indexes out of range in swapElementIndexes call!"<<std::endl;
        return;
    }

    if(index1 >= triangleAreaVecLength || index2 >= triangleAreaVecLength)
    {
        calculateTrianglesArea();
    }

    if(index1 >= collocationPointsVecLength || index2 >= collocationPointsVecLength)
    {
        calculateCollocationPoints();
    }

    if(index1 >= triangleColorsVecLength || index2 >= triangleColorsVecLength)
    {
        swapTriangleColors = false;
    }


    if(index1 >= phiSolutionVecLength || index2 >= phiSolutionVecLength || index1 >= dPhiSolutionVecLength ||
            index2 >= dPhiSolutionVecLength || index1 >= soundPressureVecLength || index2 >= soundPressureVecLength)
    {
        swapSolutions = false;
    }



    //do the swapping
    VectorTriangle tmpTriangle;
    double tmpTriangleArea;
    Eigen::Vector3d tmpCollocationPoint;
    Eigen::Vector4d tmpTriangleColor;
    std::complex<double> tmpPhiSolution;
    std::complex<double> tmpDPhiSolution;
    std::complex<double> tmpSoundPressure;

    tmpTriangle = triangles.at(index1);
    triangles[index1] = triangles.at(index2);
    triangles[index2] = tmpTriangle;

    tmpTriangleArea = trianglesArea(index1);
    trianglesArea(index1) = trianglesArea(index2);
    trianglesArea(index2) = tmpTriangleArea;

    tmpCollocationPoint = collocationPoints.at(index1);
    collocationPoints[index1] = collocationPoints.at(index2);
    collocationPoints[index2] = tmpCollocationPoint;

    if(swapTriangleColors)
    {
        tmpTriangleColor = triangleColors.at(index1);
        triangleColors[index1] = triangleColors.at(index2);
        triangleColors[index2] = tmpTriangleColor;
    }

    if(swapSolutions)
    {
        tmpPhiSolution = phiSolution(index1);
        phiSolution(index1) = phiSolution(index2);
        phiSolution(index2) = tmpPhiSolution;

        tmpDPhiSolution = dPhiSolution(index1);
        dPhiSolution(index1) = dPhiSolution(index2);
        dPhiSolution(index2) = tmpDPhiSolution;

        tmpSoundPressure = soundPressure(index1);
        soundPressure(index1) = soundPressure(index2);
        soundPressure(index2) = tmpSoundPressure;
    }
}

void BoundaryElements::reflectGeometry(ImpedancePlane reflectionPlane)
{
    geometryIsReflected = true;
    MeshFunctions::reflectTriangles(triangles, reflectionPlane.halfSpacePlanePoint, reflectionPlane.halfSpacePlaneNormal); // reflect triangles
}



