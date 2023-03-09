#ifndef WALLIMPEDANCEELEMENT_H
#define WALLIMPEDANCEELEMENT_H

#include "robinboundarycondition.h"

#include <QString>

/**
 * \class WallImpecanceElement
 * \brief This class is used to set up a boundary condition for a boundary element in an ElementSection
 *
 */
class WallImpecanceElement
{
public:
    WallImpecanceElement(){}
//    ~WallImpecanceElement();
    WallImpecanceElement(quint64 index, qint64 elementIndex, QString RefElements, QString DataFileAlias, std::complex<double> Value){
        this->index = index;
        this->elementIndex = elementIndex;
        this-> RefElements = RefElements;
        this-> Value = Value;
        this-> DataFileAlias = DataFileAlias;
    }

    quint64 index;
    qint64 elementIndex; /*!< \brief The index of the element in the referenced element section, whose boundary candition is to be changed to robinBoundaryCondition. */
    QString RefElements; /*!< \brief The referenced element section, which contains the element, whose boundary candition is to be changed to robinBoundaryCondition. */
    std::complex<double> Value;  /*!< \brief The value from the solving script gets converted into the robinBoundaryCondition in WallImpedanceSection::setUpValues(). */
    QString DataFileAlias; /*!< \brief A reference to a datafile to set up frequency dependent boundary conditions. Not implemented yet. */
    RobinBoundaryCondition robinBoundaryCondition; /*!< \brief A boundary condition. */
};
#endif // WALLIMPEDANCEELEMENT_H
