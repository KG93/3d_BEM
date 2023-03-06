#include "observationfield.h"


void ObservationField::calculateTrianglesMidpoints()
{
    triangleMidPoints.resize(triangles.size());
    for(int i=0; i<triangleMidPoints.size(); i++)
    {
        triangleMidPoints[i] = global::midpointOfTriangle(triangles.at(i));
    }
}

