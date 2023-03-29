#include "observationfield.h"


void ObservationField::calculateTrianglesMidpoints()
{
    triangleMidPoints.resize(triangles.size());
    for(int i=0; i<triangleMidPoints.size(); i++)
    {
        triangleMidPoints[i] = global::midpointOfTriangle(triangles.at(i));
    }
}

FieldSolutions ObservationField::getFieldSolutions(const QVector<ObservationField> &fields)
{
    FieldSolutions fieldSolutions;
    fieldSolutions.phiSolutions.resize(fields.size());
    fieldSolutions.soundPressures.resize(fields.size());
    for(int i=0; i<fields.size(); i++)
    {
        fieldSolutions.phiSolutions[i] = fields.at(i).phiSolution;
        fieldSolutions.soundPressures[i] = fields.at(i).soundPressure;
    }
    return fieldSolutions;
}
