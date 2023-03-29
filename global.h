#ifndef GLOBAL_H
#define GLOBAL_H

#include "vectortriangle.h"
#include "vectorquadrilateral.h"

#include <iostream>
#include <fstream>
#include <QString>
#include <QRegularExpression>
#include <QVector>
#include <QList>
#include <QPair>
#include <QVector3D>
#include <QTextStream>
#include <QFile>
#include <QRandomGenerator>

#ifdef __linux__
    #include <malloc.h>
#endif

/**
* \class global
* \brief This class contains several global constants and utility functions.
*/
class global
{
public:
/**
* \brief Set to true for debugging.
*/
static constexpr bool debugging = false;

/**
* \brief Set to true for active development.
*/
static constexpr bool activeProgramming = true;

/**
* \brief A small value used for triangle partition in the quadrature method for almost singular integrals.
*/
static constexpr double epsilon = 0.05;

/**
* \brief Return \f$ \pi \f$.
*/
static constexpr double PI = std::numbers::pi;

/**
* \brief Tiny numerical value.
*/
static constexpr double tiny = 10.0 * Eigen::NumTraits< std::complex<double> >::epsilon();

/**
* \brief The maximum integer.
*/
static constexpr double maxInt = Eigen::NumTraits<int>::highest();

/**
* \brief Call malloc_trim on Linux only.
*/
static void trimMemory();

/**
* \brief Convert a value in string form with units to a double.
*/
static double stringWithUnitsToDouble(/*const*/ QString input, const QString units, bool &validLine);

/**
* \brief Convert a value in string form with units to an integer.
*/
static int stringWithUnitsToInt(QString input, const QString units, bool &validLine);

/**
* \brief Convert a value in string form with units to an unsigned integer.
*/
static quint64 stringWithUnitsToUInt(QString input, const QString units, bool &validLine);

/**
* \brief Regex to extract a value from a line of the form identifier = value.
*/
static QRegularExpression regExWithOneValue(const QString &identifier);

/**
* \brief Regex to find an identifier in a line. Otionally extract a value from a line of the form identifier = value.
*/
static QRegularExpression regExWithOptionalValue(const QString &identifier);

/**
* \brief Remove quotes and whitespaces from string.
*/
static void removeQuotesAndWhitespace(QString &line);

/**
* \brief Remove quotes from string.
*/
static void removeQuotes(QString &line);

static QStringList stringToStringListQuotesIntact(const QString &currentLine);

/**
* \brief Compares smaller values of two intervals.
*/
static bool smallerFirst(QPair<quint64,quint64> a, QPair<quint64,quint64> b);

static void mergeIntervals(QVector<QPair<quint64,quint64>> &intervalVector);

static QVector<quint64> pairListToSinglesVec(const QVector<QPair<quint64,quint64>> &pairVec);

static QVector<QPair<quint64,quint64>> singlesVectorToPairList(const QVector<quint64> &vector);

static QVector<QPair<quint64,quint64>> diffIntervals(QVector<QPair<quint64,quint64>> &includeVector, QVector<QPair<quint64,quint64>> &excludeVector);

static void printIntervalList(const QVector<QPair<quint64,quint64>> &intervalList);

static bool upperBoundCompare(QPair<quint64,quint64> a, QPair<quint64,quint64> b);

static bool valueIsInInterval(quint64 value, const QVector<QPair<quint64,quint64>> &intervalVector);

static void printQStringList(const QStringList &stringList);

/**
* \brief Create a continuous ascending vector.
*/
static QVector<long> createContiguousIndexVector(const long startIndex, const long endIndex);

/**
* \brief Create a random permutation vector in O(n).
*/
static QVector<long> createRandomPermutationVector(const long startIndex, const long endIndex);

/**
* \brief Write phi-solution from boundary elements into vector.
*/
static Eigen::VectorXcd getTrianglesPhi(const QVector<VectorTriangle> &triangleVector);

/**
* \brief Write dphi-solution from boundary elements into vector.
*/
static Eigen::VectorXcd getTrianglesDPhi(const QVector<VectorTriangle> &triangleVector);

/**
* \brief Write phi-solution from vector into boundary elements.
*/
static void setTrianglesPhi(QVector<VectorTriangle> &triangleVector, const Eigen::VectorXcd &phiVector);

/**
* \brief Write dphi-solution from vector into boundary elements.
*/
static void setTrianglesDPhi(QVector<VectorTriangle> &triangleVector, const Eigen::VectorXcd &phiVector);

static void calculateTriangleMidPoints(QVector<VectorTriangle> &triangleVector);

static Eigen::Matrix3d triangleCoordinateMatrix(const VectorTriangle &triangle);
static Eigen::Matrix3d quadrilateralCoordinateMatrix(const VectorQuadrilateral &quadrilateral);

static Eigen::Vector3d calculateAveragePoint(const QVector<VectorTriangle> &triangleVector, const QVector<VectorQuadrilateral> &quadrilateralVector);
static double calculateMaxRelativeDistance(const QVector<VectorTriangle> &triangleVector, const Eigen::Vector3d &refNode); // calculates the  maximum distance between all the triangle nodes and the reference point
static Eigen::Vector3d calculateMaxNorm(const QVector<VectorTriangle> &triangleVector/*, const QVector<VectorQuadrilateral> &quadrilateralVector*/);
static Eigen::Vector3d triangleCompWiseMaxNorm(const VectorTriangle &triangle);
static Eigen::Vector3d calculateMinNorm(const QVector<VectorTriangle> &triangleVector/*, const QVector<VectorQuadrilateral> &quadrilateralVector*/);
static Eigen::Vector3d triangleCompWiseMinNorm(const VectorTriangle &triangle);

static quint64 indexOfMaxOfThree(double a, double b, double c);
static double maxOfThree(double a, double b, double c);
static double maxClDouble3(const Eigen::Vector3d &node);
static double minOfThree(double a, double b, double c);
static double minClDouble3(const Eigen::Vector3d &node);
static double maxOfFour(double a, double b, double c, double d);
static int allIntsUnique(int a, int b, int c, int d);

/**
* \brief Calculate the distance between two points.
*/
static double lineLength(const Eigen::Vector3d &point1, const Eigen::Vector3d &point2);

static double quadrilateralAspectRatio(const VectorQuadrilateral &parentQuadrilateral);

/**
* \brief Convert a quadrilateral to two triangles.
*/
static QVector<VectorTriangle> quadrilateralsToTriangles(const QVector<VectorQuadrilateral> &quadrilaterals);

/**
* \brief Calculate the volume and center of mass of a triangle mesh object.
*/
static QPair<Eigen::Vector3d, double> centerAndVolumeOfMassOfTriangleObject(const QVector<VectorTriangle> &triangleVector);

/**
* \brief Calculate the normal vector of a triangle.
*/
static Eigen::Vector3d normalVectorOfTriangle(const VectorTriangle &triangle);

/**
* \brief Calculate the midpoint of a triangle.
*/
static Eigen::Vector3d midpointOfTriangle(const VectorTriangle &triangle);

//static Eigen::Vector3d projectPointOnTriangle(const VectorTriangle &triangle, const Eigen::Vector3d &point, bool &inInterior, bool &onNode1, bool &onNode2, bool &onNode3, bool &onSide21, bool &onSide32, bool &onSide13); /*!< Project a point in 3d space onto a triangle. The bools show whethet the pojected point lies on a corner or on a side. (Slow version; Do not use.) */


/**
* \brief The enum represents the different regions of a triangle. Used for projectPointOnTriangleFaster().
*/
enum triangleProjectRegion {inInterior, onNode1, onNode2, onNode3, onSide21, onSide32, onSide13};
/**
* \brief Project a point in 3d space onto a triangle. Adapted from "Distance Between Point and Triangle in 3D" by David Eberly, Geometric Tools
* \param triangle The triangle onto which to project the point
* \param point The point to be projected onto the triangle
* \return the projected point on the triangle and the region of the triangle where the projected point lies.
*/
static std::pair<Eigen::Vector3d, triangleProjectRegion> projectPointOnTriangleFaster(const VectorTriangle &triangle, const Eigen::Vector3d &point);

//static void testProjectOnTriangle(const Eigen::Vector3d &point = {0,0,0}, const VectorTriangle &testTriangle = VectorTriangle({0,0,0}, {1,0,0}, {0,1,0}));

/**
* \brief Calculate the area of a triangle.
*/
static double areaOfTriangle(const VectorTriangle &triangle);

/**
* \brief Calculate the LU-decomposition of a complex matrix. L*U = A
*/
static Eigen::MatrixXcd LUDecompNoPivoting(const Eigen::MatrixXcd &A);

/**
* \brief Solve the system L*x = b, L is lower triangular for x by forward substitution.
*/
static Eigen::VectorXcd forwardSubstitution(const Eigen::MatrixXcd &lowerTriangular, Eigen::VectorXcd rightHandSide);

/**
* \brief Solve the system U*x = b, U is upper triangular for x by backward substitution.
*/
static Eigen::VectorXcd backwardSubstitution(const Eigen::MatrixXcd &upperTriangular, Eigen::VectorXcd rightHandSide);

/**
* \brief Write an Eigen matrix into a file.
*/
static void printMatrixToFile(const QString fileName, const Eigen::MatrixXcd &A);

/**
* \brief Print vector via std::cout.
*/
// templated functions in the header file
template<typename container>
static void printVector(const container &vector);
};

template<typename container>
void global::printVector(const container &vector)
{
    for(typename container::const_iterator it = vector.begin(); it != vector.end(); ++it)
    {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
}

#endif // GLOBAL_H
