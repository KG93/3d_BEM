#include "meshfunctions.h"

MeshFunctions::MeshFunctions()
{

}
void MeshFunctions::swapTriangleNormals(QVector<VectorTriangle>& elementsVectorTriangles)
{
    for(int i = 0; i < elementsVectorTriangles.length(); i++)
    {
        VectorTriangle tmp = elementsVectorTriangles.at(i);
        elementsVectorTriangles[i] = VectorTriangle(tmp.elementIndex, tmp.node3, tmp.node2, tmp.node1, -tmp.normal);
    }
}

void MeshFunctions::swapQuadrilateralNormals(QVector<VectorQuadrilateral>& elementsVectorQuadrilaterals)
{
    for(int i=0; i<elementsVectorQuadrilaterals.length(); i++)
    {
        VectorQuadrilateral tmp = elementsVectorQuadrilaterals.at(i);
        elementsVectorQuadrilaterals[i] = VectorQuadrilateral(tmp.elementIndex, tmp.node4, tmp.node3, tmp.node2, tmp.node1);
    }
}

QVector<VectorTriangle> MeshFunctions::refineTriangles(const QVector<VectorTriangle>& elementsVectorTriangles, const double EdgeLength)
{
    QVector<VectorTriangle> returnVector;
    int vectorLength = elementsVectorTriangles.length();
    for(int i = 0; i < vectorLength; i++)
    {
        VectorTriangle tmp = elementsVectorTriangles.at(i);
        double edgeLength1 = global::lineLength(tmp.node1,tmp.node2);
        double edgeLength2 = global::lineLength(tmp.node2,tmp.node3);
        double edgeLength3 = global::lineLength(tmp.node3,tmp.node1);
        double maxEdge = 0;
//        Eigen::Vector3d maxEdgeMidpoint;
        int max = global::indexOfMaxOfThree(edgeLength1, edgeLength2, edgeLength3);
        switch(max) {
            case 1 : maxEdge = edgeLength1;
                     break;
            case 2 : maxEdge = edgeLength2;
                     break;
            case 3 : maxEdge = edgeLength3;
                     break;
        }
        if(maxEdge <= EdgeLength)
        {
            returnVector.push_back(tmp);
            continue;
        }
        else
        {
            recursiveRefineTriangle(tmp, max, returnVector, EdgeLength);
        }
    }
    return returnVector;
}

void MeshFunctions::recursiveRefineTriangle(VectorTriangle parentTriangle, int longestEdgeIndex, QVector<VectorTriangle>& vectorToStoreTriangles, const double EdgeLength)
{
    double maxEdge = 0;
    switch(longestEdgeIndex) {
        case 1 : maxEdge=global::lineLength(parentTriangle.node1,parentTriangle.node2);
                 break;
        case 2 : maxEdge=global::lineLength(parentTriangle.node2,parentTriangle.node3);
                 break;
        case 3 : maxEdge=global::lineLength(parentTriangle.node3,parentTriangle.node1);
                 break;
    }
    if(maxEdge <= EdgeLength)
    {
        vectorToStoreTriangles.push_back(parentTriangle);
        return;
    }
    else
    {
        VectorTriangle vectorTriangle1;
        VectorTriangle vectorTriangle2;
        Eigen::Vector3d maxEdgeMidpoint;

        switch(longestEdgeIndex) {
            case 1 :  maxEdgeMidpoint=middleOfLine(parentTriangle.node1,parentTriangle.node2);
                      vectorTriangle1=VectorTriangle(parentTriangle.elementIndex,parentTriangle.node1,maxEdgeMidpoint,parentTriangle.node3,parentTriangle.normal);
                      vectorTriangle2=VectorTriangle(parentTriangle.elementIndex,maxEdgeMidpoint,parentTriangle.node2,parentTriangle.node3,parentTriangle.normal);
                     break;
            case 2 : maxEdgeMidpoint=middleOfLine(parentTriangle.node2,parentTriangle.node3);
                     vectorTriangle1=VectorTriangle(parentTriangle.elementIndex,parentTriangle.node1,parentTriangle.node2,maxEdgeMidpoint,parentTriangle.normal);
                     vectorTriangle2=VectorTriangle(parentTriangle.elementIndex,maxEdgeMidpoint,parentTriangle.node3,parentTriangle.node1,parentTriangle.normal);
                     break;
            case 3 : maxEdgeMidpoint=middleOfLine(parentTriangle.node3,parentTriangle.node1);
                     vectorTriangle1=VectorTriangle(parentTriangle.elementIndex,parentTriangle.node1,parentTriangle.node2,maxEdgeMidpoint,parentTriangle.normal);
                     vectorTriangle2=VectorTriangle(parentTriangle.elementIndex,maxEdgeMidpoint,parentTriangle.node2,parentTriangle.node3,parentTriangle.normal);
                     break;
        }
//        vectorTriangle1.normal=parentTriangle.normal;
//        vectorTriangle2.normal=parentTriangle.normal;
//        vectorTriangle2.robinBoundaryCondition=parentTriangle.robinBoundaryCondition;
        double edgeLength1=global::lineLength(vectorTriangle1.node1,vectorTriangle1.node2);
        double edgeLength2=global::lineLength(vectorTriangle1.node2,vectorTriangle1.node3);
        double edgeLength3=global::lineLength(vectorTriangle1.node3,vectorTriangle1.node1);
        int max=global::indexOfMaxOfThree(edgeLength1, edgeLength2, edgeLength3);
        recursiveRefineTriangle(vectorTriangle1,max,vectorToStoreTriangles,EdgeLength);

         edgeLength1=global::lineLength(vectorTriangle2.node1,vectorTriangle2.node2);
         edgeLength2=global::lineLength(vectorTriangle2.node2,vectorTriangle2.node3);
         edgeLength3=global::lineLength(vectorTriangle2.node3,vectorTriangle2.node1);
         max=global::indexOfMaxOfThree(edgeLength1, edgeLength2, edgeLength3);
         recursiveRefineTriangle(vectorTriangle2,max,vectorToStoreTriangles,EdgeLength);
    }
}

QVector<VectorTriangle> MeshFunctions::quadrilateralsToTriangles(const QVector<VectorQuadrilateral>& quadrilaterals)
{
    QVector<VectorTriangle> triangles;
    for(int i = 0; i < quadrilaterals.length(); i++) //split quadrilaterals with too high aspect ratio recursively in two quadrilaterals
    {
        VectorQuadrilateral parentQuadrilateral=quadrilaterals.at(i);

        double diagonal1=global::lineLength(parentQuadrilateral.node1,parentQuadrilateral.node3);
        double diagonal2=global::lineLength(parentQuadrilateral.node2,parentQuadrilateral.node4);
        VectorTriangle vectorTriangle1;
        VectorTriangle vectorTriangle2;
        if(diagonal1 > diagonal2) //split quadrilateral by shorter diagonal
        {
            vectorTriangle1=VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,parentQuadrilateral.node2,parentQuadrilateral.node4);
            vectorTriangle2=VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node2,parentQuadrilateral.node3,parentQuadrilateral.node4);
        }
        else
        {
            vectorTriangle1=VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node4,parentQuadrilateral.node1,parentQuadrilateral.node3);
            vectorTriangle2=VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,parentQuadrilateral.node2,parentQuadrilateral.node3);
        }
//        vectorTriangle1.robinBoundaryCondition=parentQuadrilateral.robinBoundaryCondition;
//        vectorTriangle2.robinBoundaryCondition=parentQuadrilateral.robinBoundaryCondition;
        triangles.push_back(vectorTriangle1);
        triangles.push_back(vectorTriangle2);
    }
    return triangles;
}

QVector<VectorTriangle> MeshFunctions::refineQuadrilaterals(const QVector<VectorQuadrilateral>& elementsVectorQuadrilaterals, const double EdgeLength)
{
    QVector<VectorQuadrilateral> refinedQuadrilaterals;
    for(int i = 0; i < elementsVectorQuadrilaterals.length(); i++) //split quadrilaterals with too high aspect ratio recursively in two quadrilaterals
    {
        VectorQuadrilateral parentQuadrilateral=elementsVectorQuadrilaterals.at(i);

        double edgeLength1=global::lineLength(parentQuadrilateral.node1,parentQuadrilateral.node2);
        double edgeLength2=global::lineLength(parentQuadrilateral.node2,parentQuadrilateral.node3);
        double edgeLength3=global::lineLength(parentQuadrilateral.node3,parentQuadrilateral.node4);
        double edgeLength4=global::lineLength(parentQuadrilateral.node4,parentQuadrilateral.node1);
        QPair<double,int> longestEdgeWithIndex=maxOfFour(edgeLength1, edgeLength2, edgeLength3, edgeLength4);
        double aspectRatio=global::quadrilateralAspectRatio(parentQuadrilateral);
        if(aspectRatio>3 && longestEdgeWithIndex.first>EdgeLength) //refine Quadrilateral if aspect ratio is high end the longer edge is longer than the specified max edge length
        {
            recursiveRefineQuadrilateral(parentQuadrilateral, longestEdgeWithIndex.second, refinedQuadrilaterals,EdgeLength);
        }
        else
        {
            refinedQuadrilaterals.push_back(parentQuadrilateral);
        }
    }
    QVector<VectorTriangle> triangles;
    triangles = global::quadrilateralsToTriangles(refinedQuadrilaterals);
    refinedQuadrilaterals.clear();
//    triangles=refineTriangles(triangles,EdgeLength);
    return triangles;
}

void MeshFunctions::recursiveRefineQuadrilateral(VectorQuadrilateral parentQuadrilateral, int longestEdgeIndex, QVector<VectorQuadrilateral>& refinedQuadrilateralsStorage, const double EdgeLength)
{
    VectorQuadrilateral vectorQuadrilateral1;
    VectorQuadrilateral vectorQuadrilateral2;
    Eigen::Vector3d maxEdgeMidpoint;
    Eigen::Vector3d oppositeEdgeMidpoint;

    switch(longestEdgeIndex) {
        case 1 :  maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node1,parentQuadrilateral.node2);
                  oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node3,parentQuadrilateral.node4);
                  vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node4);
                  vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node2,parentQuadrilateral.node3,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
        case 2 : maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node2,parentQuadrilateral.node3);
                 oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node4,parentQuadrilateral.node1);
                 vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node2,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node1);
                 vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node3,parentQuadrilateral.node4,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
        case 3 : maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node3,parentQuadrilateral.node4);
                 oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node1,parentQuadrilateral.node2);
                 vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node3,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node2);
                 vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node4,parentQuadrilateral.node1,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
        case 4 : maxEdgeMidpoint=middleOfLine(parentQuadrilateral.node4,parentQuadrilateral.node1);
                 oppositeEdgeMidpoint=middleOfLine(parentQuadrilateral.node2,parentQuadrilateral.node3);
                 vectorQuadrilateral1=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node4,maxEdgeMidpoint,oppositeEdgeMidpoint,parentQuadrilateral.node3);
                 vectorQuadrilateral2=VectorQuadrilateral(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,parentQuadrilateral.node2,oppositeEdgeMidpoint,maxEdgeMidpoint);
                 break;
    }
//    vectorQuadrilateral1.robinBoundaryCondition=parentQuadrilateral.robinBoundaryCondition;
//    vectorQuadrilateral2.robinBoundaryCondition=parentQuadrilateral.robinBoundaryCondition;

    double aspectRatio=global::quadrilateralAspectRatio(vectorQuadrilateral1);
    if(aspectRatio>3 ) //refine Quadrilateral if aspect ratio is high end the longer edge is longer than the specified max edge length
    {
        double edgeLength1=global::lineLength(vectorQuadrilateral1.node1,vectorQuadrilateral1.node2);
        double edgeLength2=global::lineLength(vectorQuadrilateral1.node2,vectorQuadrilateral1.node3);
        double edgeLength3=global::lineLength(vectorQuadrilateral1.node3,vectorQuadrilateral1.node4);
        double edgeLength4=global::lineLength(vectorQuadrilateral1.node4,vectorQuadrilateral1.node1);
        QPair<double,int> longestEdgeWithIndex=maxOfFour(edgeLength1, edgeLength2, edgeLength3, edgeLength4);
        if(longestEdgeWithIndex.first>EdgeLength)
        {
            recursiveRefineQuadrilateral(vectorQuadrilateral1, longestEdgeWithIndex.second, refinedQuadrilateralsStorage,EdgeLength);
        }
        else
        {
            refinedQuadrilateralsStorage.push_back(vectorQuadrilateral1);
        }
    }
    else
    {
        refinedQuadrilateralsStorage.push_back(vectorQuadrilateral1);
    }
    double aspectRatio2 = global::quadrilateralAspectRatio(vectorQuadrilateral2);
    if(aspectRatio2 > 3)
    {
        double edgeLength1=global::lineLength(vectorQuadrilateral2.node1,vectorQuadrilateral2.node2);
        double edgeLength2=global::lineLength(vectorQuadrilateral2.node2,vectorQuadrilateral2.node3);
        double edgeLength3=global::lineLength(vectorQuadrilateral2.node3,vectorQuadrilateral2.node4);
        double edgeLength4=global::lineLength(vectorQuadrilateral2.node4,vectorQuadrilateral2.node1);
        QPair<double,int> longestEdgeWithIndex=maxOfFour(edgeLength1, edgeLength2, edgeLength3, edgeLength4);
        if(longestEdgeWithIndex.first>EdgeLength)
        {
        recursiveRefineQuadrilateral(vectorQuadrilateral2, longestEdgeWithIndex.second, refinedQuadrilateralsStorage,EdgeLength);
        }
        else
        {
            refinedQuadrilateralsStorage.push_back(vectorQuadrilateral2);
        }

    }
    else
    {
        refinedQuadrilateralsStorage.push_back(vectorQuadrilateral2);
    }
}

Eigen::Vector3d MeshFunctions::middleOfLine(Eigen::Vector3d point1, Eigen::Vector3d point2)
{    
    return (point1 + point2) / 2;
}

QPair<double,quint64> MeshFunctions::maxOfFour(double a, double b, double c,double d)
{
    if(a>b)
    {
        if(a>c)
        {
            if(a>d)
            {
                return qMakePair(a,1);
            }
            else
            {
                return qMakePair(d,4);
            }
        }
        else
        {
            if(c>d)
            {
                return qMakePair(c,3);
            }
            else
            {
                return qMakePair(d,4);
            }
        }
    }
    else
    {
        if(b>c)
        {
            if(b>d)
            {
                return qMakePair(b,2);
            }
            else
            {
                return qMakePair(d,4);
            }
        }
        else
        {
            if(c>d)
            {
                return qMakePair(c,3);
            }
            else
            {
                return qMakePair(d,4);
            }
        }

    }
}

void MeshFunctions::transformTriangles(QVector<VectorTriangle> &elementsVectorTriangles, const Eigen::Matrix4d transformationMatrix)
{
    Eigen::Vector4d homogNode;
    for(int i=0;i<elementsVectorTriangles.length();i++)
    {
        homogNode=homogVec(elementsVectorTriangles[i].node1);
        homogNode=transformationMatrix*homogNode;
        elementsVectorTriangles[i].node1=eigenVec(homogNode);

        homogNode=homogVec(elementsVectorTriangles[i].node2);
        homogNode=transformationMatrix*homogNode;
        elementsVectorTriangles[i].node2=eigenVec(homogNode);

        homogNode=homogVec(elementsVectorTriangles[i].node3);
        homogNode=transformationMatrix*homogNode;
        elementsVectorTriangles[i].node3=eigenVec(homogNode);
    }
}

void MeshFunctions::transformQuadrilaterals(QVector<VectorQuadrilateral> &elementsVectorQuadrilaterals,const Eigen::Matrix4d transformationMatrix)
{
    Eigen::Vector4d homogNode;
    for(int i=0;i<elementsVectorQuadrilaterals.length();i++)
    {
        homogNode=homogVec( elementsVectorQuadrilaterals[i].node1);
        homogNode=transformationMatrix*homogNode;
        elementsVectorQuadrilaterals[i].node1=eigenVec(homogNode);

        homogNode=homogVec( elementsVectorQuadrilaterals[i].node2);
        homogNode=transformationMatrix*homogNode;
        elementsVectorQuadrilaterals[i].node2=eigenVec(homogNode);

        homogNode=homogVec( elementsVectorQuadrilaterals[i].node3);
        homogNode=transformationMatrix*homogNode;
        elementsVectorQuadrilaterals[i].node3=eigenVec(homogNode);

        homogNode=homogVec( elementsVectorQuadrilaterals[i].node4);
        homogNode=transformationMatrix*homogNode;
        elementsVectorQuadrilaterals[i].node4=eigenVec(homogNode);
    }
}

void MeshFunctions::reflectTriangles(QVector<VectorTriangle> &elementsVectorTriangles, Eigen::Vector3d planePoint, Eigen::Vector3d planeNormal)
{
    planeNormal = planeNormal.normalized();
    double d = -planeNormal.dot(planePoint);
    Eigen::Matrix4d reflectionMatrix;
    reflectionMatrix << 1 - 2.0 * std::pow(planeNormal(0),2), - 2.0 * planeNormal(0) * planeNormal(1), - 2.0 * planeNormal(0) * planeNormal(2), - 2.0 * planeNormal(0) * d,
                        - 2.0 * planeNormal(0) * planeNormal(1), 1 - 2.0 * std::pow(planeNormal(1),2), - 2.0 * planeNormal(1) * planeNormal(2), - 2.0 * planeNormal(1) * d,
                        - 2.0 * planeNormal(0) * planeNormal(2), - 2.0 * planeNormal(1) * planeNormal(2), 1 - 2.0 * std::pow(planeNormal(2),2), - 2.0 * planeNormal(2) * d,
                        0, 0, 0, 1;

    Eigen::Vector4d homogNode;
    Eigen::Vector4d tmpHomogNode;
    for(int i=0;i<elementsVectorTriangles.length();i++)
    {
        homogNode=homogVec(elementsVectorTriangles[i].node1);
        homogNode=reflectionMatrix*homogNode;

        tmpHomogNode=homogVec(elementsVectorTriangles[i].node3);

        elementsVectorTriangles[i].node3=eigenVec(homogNode);

        homogNode=homogVec(elementsVectorTriangles[i].node2);
        homogNode=reflectionMatrix*homogNode;
        elementsVectorTriangles[i].node2=eigenVec(homogNode);

        homogNode=reflectionMatrix*tmpHomogNode;
        elementsVectorTriangles[i].node1=eigenVec(homogNode);

        elementsVectorTriangles[i].normal = global::normalVectorOfTriangle(elementsVectorTriangles.at(i));
        elementsVectorTriangles[i].triangleMidpoint = global::midpointOfTriangle(elementsVectorTriangles.at(i));
    }
}

Eigen::Matrix4d MeshFunctions::qtToEigenMatrix(const QMatrix4x4& qtMatrix)
{
    Eigen::Matrix4d eigenMatrix;
    eigenMatrix << qtMatrix(0,0),qtMatrix(0,1),qtMatrix(0,2),qtMatrix(0,3),qtMatrix(1,0),qtMatrix(1,1),qtMatrix(1,2),qtMatrix(1,3),qtMatrix(2,0),qtMatrix(2,1),qtMatrix(2,2),qtMatrix(2,3),qtMatrix(3,0),qtMatrix(3,1),qtMatrix(3,2),qtMatrix(3,3);
    return eigenMatrix;
}

