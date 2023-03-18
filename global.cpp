#include "global.h"

void global::trimMemory() /*!< Call malloc_trim on Linux. */
{
    #ifdef __linux__
        malloc_trim(0);
    #endif
}

double global::stringWithUnitsToDouble(QString input, const QString units, bool& validLine)
{
    input.remove("\\s");
    QRegularExpression removeChars(units+"$");
    input.remove(removeChars);
    QRegularExpression unitPrefixRegEx = QRegularExpression("([nµumkKMG])$"); /*nµumkKMG*/
    double unitPrefixFactor=1;
//    unitPrefixRegEx.indexIn(input);
    QRegularExpressionMatch match = unitPrefixRegEx.match(input);

    if(unitPrefixRegEx.captureCount()==1)
    {
//        QString actualUnitPrefix = unitPrefixRegEx.cap(1);
        QString actualUnitPrefix = match.captured(1);
        if ( actualUnitPrefix == "n" ) {
            (unitPrefixFactor=1e-9);
        } else
        if ( actualUnitPrefix == "u" || actualUnitPrefix == "µ" ) {
            (unitPrefixFactor=1e-6);
        } else
        if ( actualUnitPrefix == "m" ) {
            (unitPrefixFactor=1e-3);
        } else
        if ( actualUnitPrefix == "K" || actualUnitPrefix == "k") {
            (unitPrefixFactor=1e3);
        } else
        if ( actualUnitPrefix == "M" ) {
            (unitPrefixFactor=1e6);
        } else
        if ( actualUnitPrefix == "G" ) {
            (unitPrefixFactor=1e9);
        }
        input.remove(unitPrefixRegEx);
    }
    bool ok = false;
    double value=input.toDouble(&ok);
    if(input.isEmpty())
    {
        value=1;
        ok=true;
    }
    validLine=ok;
    value=value*unitPrefixFactor;

    return value;
}

int global::stringWithUnitsToInt(QString input, const QString units, bool& validLine)
{
    input.remove("\\s");
//    QRegularExpression removeChars("\\s*"+units+"\\s*$");
    QRegularExpression removeChars(units+"$");
    input.remove(removeChars);
    QRegularExpression unitPrefixRegEx = QRegularExpression("([nµumkKMG])$"); /*nµumkKMG*/
    double unitPrefixFactor=1;
    QRegularExpressionMatch match = unitPrefixRegEx.match(input);
    if(unitPrefixRegEx.captureCount()==1)
    {
        QString actualUnitPrefix = match.captured(1);
        if ( actualUnitPrefix == "n" ) {
            (unitPrefixFactor=1e-9);
        } else
        if ( actualUnitPrefix == "u" || actualUnitPrefix == "µ" ) {
            (unitPrefixFactor=1e-6);
        } else
        if ( actualUnitPrefix == "m" ) {
            (unitPrefixFactor=1e-3);
        } else
        if ( actualUnitPrefix == "K" || actualUnitPrefix == "k") {
            (unitPrefixFactor=1e3);
        } else
        if ( actualUnitPrefix == "M" ) {
            (unitPrefixFactor=1e6);
        } else
        if ( actualUnitPrefix == "G" ) {
            (unitPrefixFactor=1e9);
        }
        input.remove(unitPrefixRegEx);
    }
    bool ok = false;
    int value=input.toInt(&ok);
    if(input.isEmpty())
    {
        value=1;
        ok=true;
    }
    validLine=ok;
    value=value*unitPrefixFactor;
    return value;
}

quint64 global::stringWithUnitsToUInt(QString input, const QString units, bool& validLine)
{
    input.remove("\\s");
//    QRegularExpression removeChars("\\s*"+units+"\\s*$");
    QRegularExpression removeChars(units+"$");
    input.remove(removeChars);
    QRegularExpression unitPrefixRegEx = QRegularExpression("([nµumkKMG])$"); /*nµumkKMG*/
    double unitPrefixFactor=1;
//    unitPrefixRegEx.indexIn(input);
    QRegularExpressionMatch match = unitPrefixRegEx.match(input);
    if(unitPrefixRegEx.captureCount()==1)
    {
//        QString actualUnitPrefix = unitPrefixRegEx.cap(1);
        QString actualUnitPrefix = match.captured(1);
        if ( actualUnitPrefix == "n" ) {
            (unitPrefixFactor=1e-9);
        } else
        if ( actualUnitPrefix == "u" || actualUnitPrefix == "µ" ) {
            (unitPrefixFactor=1e-6);
        } else
        if ( actualUnitPrefix == "m" ) {
            (unitPrefixFactor=1e-3);
        } else
        if ( actualUnitPrefix == "K" || actualUnitPrefix == "k") {
            (unitPrefixFactor=1e3);
        } else
        if ( actualUnitPrefix == "M" ) {
            (unitPrefixFactor=1e6);
        } else
        if ( actualUnitPrefix == "G" ) {
            (unitPrefixFactor=1e9);
        }
        input.remove(unitPrefixRegEx);
    }
    bool ok = false;
    quint64 value=input.toUInt(&ok);
    if(input.isEmpty())
    {
        value=1;
        ok=true;
    }
    validLine = ok;
    value = value * unitPrefixFactor;
    return value;
}

QRegularExpression global::regExWithOneValue(const QString& identifier)
{
//    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][\\w\\s]+[\"\']|\\w+)", QRegularExpression::CaseInsensitiveOption);
    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][^\"\']+[\"\']|.+$)", QRegularExpression::CaseInsensitiveOption);
}

QRegularExpression global::regExWithOptionalValue(const QString& identifier)
{
    return QRegularExpression("^\\s*"+identifier+"\\s*={1}\\s*([\"\'][^\"\']+[\"\']|.+$)?", QRegularExpression::CaseInsensitiveOption);
}

void global::removeQuotesAndWhitespace(QString& line)
{
    QRegularExpression removeChars("[\"\'\\s]");
    line.remove(removeChars);
}

void global::removeQuotes(QString& line)
{
    QRegularExpression removeChars("[\"\']");
    line.remove(removeChars);
}

QStringList global::stringToStringListQuotesIntact(const QString& currentLine)
{
    QStringList list;
    QRegularExpression strings( "[\"\']([^\"\']*)[\"\']|(\\S+)"); //identify any string between "or'-quotes or any string which isn't whitespace
    QRegularExpressionMatchIterator iterator = strings.globalMatch(currentLine);
    while (iterator.hasNext()) {
        QRegularExpressionMatch match = iterator.next();
        list.push_back(match.captured(0));
    }
    return list;
}

bool global::smallerFirst(QPair<quint64,quint64> a, QPair<quint64,quint64> b)
{
    return a.first > b.first;
}

void global::mergeIntervals(QVector<QPair<quint64,quint64>>& intervalList)
{
    for(int i=0;i<intervalList.length();i++)
    {
        if(intervalList.at(i).first > intervalList.at(i).second)
        {
            quint64 first=intervalList.at(i).first;
            quint64 second=intervalList.at(i).second;
            intervalList.replace(i,qMakePair(first,second));
        }
    }
    std::sort(intervalList.begin(),intervalList.end(),smallerFirst);

    int index = 0; // Stores index of last element
    // in output array (modified arr[])

    // Traverse all input Intervals
    for (int i = 0; i < intervalList.length(); i++)
    {
        // If this is not first Interval and overlaps
        // with the previous one
        if (index != 0 && intervalList.at(index-1).first <= intervalList.at(i).second)
        {
            while (index != 0 && intervalList.at(index-1).first <= intervalList.at(i).second)
            {
                // Merge previous and current Intervals
                quint64 second = std::max(intervalList.at(index-1).second, intervalList.at(i).second);
                quint64 first = std::min(intervalList.at(index-1).first, intervalList.at(i).first);
                intervalList.replace(index-1,qMakePair(first,second));
                index--;
            }
        }
        else // Doesn't overlap with previous, add to
            // solution
            intervalList.replace(index,intervalList.at(i));

        index++;
    }
    quint64 vectorLength=index;
//    intervalList=intervalList.mid(0,vectorLength);
    intervalList.resize(vectorLength);
    std::reverse(intervalList.begin(),intervalList.end());

    for(int i=0;i<intervalList.length()-1;i++)
    {
        if(intervalList.at(i).second==(intervalList.at(i+1).first)-1)
        {
            quint64 first=intervalList.at(i).first;
            quint64 second=intervalList.at(i+1).second;
            intervalList.removeAt(i+1);
            intervalList.replace(i,qMakePair(first,second));
            i--;
        }
    }

    std::cout <<"interval list: ";
    printIntervalList(intervalList);
}

QVector<quint64> global::pairListToSinglesVec(const QVector<QPair<quint64,quint64>>& pairVec)
{
    QVector<quint64>singlesVector;
    singlesVector.reserve(pairVec.length() * 2);
    for(int i=0;i<pairVec.length();i++)
    {
        singlesVector.push_back(pairVec.at(i).first);
        singlesVector.push_back(pairVec.at(i).second);
    }
    return singlesVector;
}

QVector<QPair<quint64,quint64>> global::singlesVectorToPairList(const QVector<quint64> &vector)
{
    QVector<QPair<quint64,quint64>> pairList;
    pairList.reserve(vector.length()/2);
    for(int i=0;i<vector.length()/2;i++)
    {
        pairList.push_back(qMakePair(vector.at(i),vector.at(i+1)));
        i++;
    }
    return pairList;
}


QVector<QPair<quint64,quint64>> global::diffIntervals(QVector<QPair<quint64,quint64>> &includeVector, QVector<QPair<quint64,quint64>> &excludeVector)
{
    QVector<QPair<quint64,quint64>> diffVector;
    mergeIntervals(includeVector); //returns list of ordered disjunct intervals with max individual length
    mergeIntervals(excludeVector);
    if(excludeVector.isEmpty())
    {
        return includeVector;
    }
    int includeIndex = 0;
    int excludeIndex = 0;

    quint64 includeStart;
    quint64 includeEnd;
    quint64 excludeStart;
    quint64 excludeEnd;
    std::cout<<"includeList length: "<<includeVector.length()<<std::endl;
    std::cout<<"excludeList length: "<<excludeVector.length()<<std::endl;
    while(includeIndex<includeVector.length())
    {
        if(excludeIndex==excludeVector.length())
        {
            excludeIndex=excludeVector.length()-1;
        }
//        std::cout<<"includeIndex: "<<includeIndex<<std::endl;
//        std::cout<<"excludeIndex: "<<excludeIndex<<std::endl;
//        std::cout<<"includeStart: "<<includeStart<<std::endl;
//        std::cout<<"includeEnd: "<<includeEnd<<std::endl;
//        std::cout<<"excludeStart: "<<excludeStart<<std::endl;
//        std::cout<<"excludeEnd: "<<excludeEnd<<std::endl;

        includeStart=includeVector.at(includeIndex).first;
        includeEnd=includeVector.at(includeIndex).second;
        excludeStart=excludeVector.at(excludeIndex).first;
        excludeEnd=excludeVector.at(excludeIndex).second;

        if(excludeEnd<includeStart) //exclude interval does not overlap with include interval
        {
            if(excludeIndex<excludeVector.length()-1) //haven't reached the last exclude interval yet, therefore the loop goes to the next higher one
            {
                excludeIndex++;
                continue;
            }
            else //reached last excludeInterval
            {
                diffVector.push_back(qMakePair(includeStart,includeEnd)); //current iclude interval doesn't overlap with the last exclude interval
                includeIndex++; //therefore it is added to the difference include list and the loop goes to the next include interval
                continue;
            }
        }
        if(excludeStart>includeEnd) //exclude interval does not overlap with include interval
        {
            diffVector.push_back(qMakePair(includeStart,includeEnd));
            includeIndex++;
            continue;
        }
        if(excludeEnd>=includeStart) //exclude interval overlaps with include
        {
            if(excludeStart<=includeStart)//exclude interval overlaps with include interval from the left
            {
                if(excludeEnd<includeEnd) //exclude interval overlaps the left part of the include interval
                {
                    includeVector[includeIndex].first=excludeEnd+1; //difference of current include and exclude interval
                    excludeIndex++; //loop goes to next exclude interval
                    continue;
                }
                else //exclude interval overlaps the whole include interval
                {
                    includeIndex++; //loop goes to next include interval
                    continue;
                }

            }
            else // excludeStart > includeStart //exclude interval overlaps at least middle part of the include interval, but not the leftmost
            {
                if(excludeEnd<includeEnd)
                {
                    diffVector.push_back(qMakePair(includeStart,excludeStart-1)); //leftmost part of include interval is added to diffVector
                    includeVector[includeIndex].first=excludeEnd+1; //current include interval size is adapted
                    excludeIndex++; //loop goes to next exclude interval
                    continue;
                }
                else //the middle and right parts of the include interval are overlapped
                {
                    diffVector.push_back(qMakePair(includeStart,std::min(includeEnd,excludeStart-1)));
                    includeIndex++; //loop goes to next include interval
                    continue;
                }
            }
        }        
    }
    std::cout<<"diffinterval: ";
    printIntervalList(diffVector);
    std::cout<<std::endl;

    return diffVector;
}

void global::printIntervalList(const QVector<QPair<quint64,quint64>>& intervalList)
{
    for(int i=0; i<intervalList.length(); i++)
    {
        std::cout<<intervalList.value(i).first<<" "<<intervalList.value(i).second<<" ";
    }
    std::cout<<std::endl;
}

bool global::upperBoundCompare(QPair<quint64,quint64> a, QPair<quint64,quint64> b)
{
    return a.second < b.second;
}

bool global::valueIsInInterval(quint64 value,const QVector<QPair<quint64,quint64>>& intervalVector)
{
//    auto interval=std::lower_bound(intervalVector.begin(), intervalVector.end(), qMakePair(value,value), upperBoundCompare);
    for(int i=0;i<intervalVector.size();i++)
    {
        if(value>=intervalVector.at(i).first && value<=intervalVector.at(i).second)
        {
            return true;
        }
    }
    return false;
//    return (value >= interval->first && value <= interval->second);
}

void global::printQStringList(const QStringList &stringList)
{
   for(int i=0; i<stringList.length(); i++)
   {
        std::cout << stringList.value(i).toUtf8().constData() << " ";
   }
   std::cout << std::endl;
}

QVector<long> global::createContiguousIndexVector(const long startIndex, const long endIndex)
{
    if(startIndex<0 || startIndex>endIndex)
    {
        std::cerr<<"Invalid indexes in createContiguousIndexVector call."<<std::endl;
        return QVector<long>();
    }
    long numberOfIndexes = endIndex - startIndex + 1;
    QVector<long> indexes(numberOfIndexes);
    for(int i=0; i<numberOfIndexes; i++)
    {
        indexes[i] = startIndex + i;
    }
    return indexes;
}

QVector<long> global::createRandomPermutationVector(const long startIndex, const long endIndex)
{
    QVector<long> indexes = createContiguousIndexVector(startIndex, endIndex);
    long numberOfIndexes = endIndex - startIndex + 1;

    for(long i = 0; i < numberOfIndexes - 1; i++)
    {
        long randomRestIndex =  QRandomGenerator::system()->bounded((qint32) (numberOfIndexes - i)) + i;
        indexes.swapItemsAt(i, randomRestIndex);
    }
    return indexes;
}

Eigen::VectorXcd global::getTrianglesPhi(const QVector<VectorTriangle>& triangleVector)
{
    long numberOfTriangles = triangleVector.size();
    Eigen::VectorXcd phiVector(numberOfTriangles);
    for(long i = 0; i <  numberOfTriangles; i++)
    {
        phiVector(i) = triangleVector.at(i).phi;
    }
    return phiVector;
}

Eigen::VectorXcd global::getTrianglesDPhi(const QVector<VectorTriangle>& triangleVector)
{
    long numberOfTriangles = triangleVector.size();
    Eigen::VectorXcd dPhiVector(numberOfTriangles);
    for(long i = 0; i <  numberOfTriangles; i++)
    {
        dPhiVector(i) = triangleVector.at(i).dPhi;
    }
    return dPhiVector;
}

void global::setTrianglesPhi(QVector<VectorTriangle>& triangleVector, const Eigen::VectorXcd& phiVector)
{
    if(triangleVector.size() != phiVector.size())
    {
        return;
    }
    else
    {
        for(long i = 0; i < triangleVector.size(); i++)
        {
            triangleVector[i].phi = phiVector(i);
        }
    }
}

void global::setTrianglesDPhi(QVector<VectorTriangle>& triangleVector, const Eigen::VectorXcd& dPhiVector)
{
    if(triangleVector.size() != dPhiVector.size())
    {
        return;
    }
    else
    {
        for(long i = 0; i < triangleVector.size(); i++)
        {
            triangleVector[i].dPhi = dPhiVector(i);
        }
    }
}

void global::calculateTriangleMidPoints(QVector<VectorTriangle> &triangleVector)
{
    for(int i=0; i<triangleVector.size(); i++)
    {
        triangleVector[i].triangleMidpoint = midpointOfTriangle(triangleVector.at(i));
    }
}

Eigen::Matrix3d global::triangleCoordinateMatrix(const VectorTriangle &triangle)
{
    Eigen::Matrix3d coordinateMatrix;
    coordinateMatrix.col(0) = triangle.node1;
    coordinateMatrix.col(1) = triangle.node2;
    coordinateMatrix.col(2) = triangle.node3;
    return coordinateMatrix;
}

Eigen::Matrix3d global::quadrilateralCoordinateMatrix(const VectorQuadrilateral &quadrilateral)
{
    Eigen::MatrixXd coordinateMatrix(3, 4);
    coordinateMatrix.col(0) = quadrilateral.node1;
    coordinateMatrix.col(1) = quadrilateral.node2;
    coordinateMatrix.col(2) = quadrilateral.node3;
    coordinateMatrix.col(3) = quadrilateral.node4;
    return coordinateMatrix;
}

Eigen::Vector3d global::calculateAveragePoint(const QVector<VectorTriangle> &triangleVector, const QVector<VectorQuadrilateral> &quadrilateralVector)
{
    Eigen::MatrixXd coordinates(3, 3 * triangleVector.length() + 4 * quadrilateralVector.length());
    for(int i = 0; i < triangleVector.length(); i++)
    {
        coordinates.block(0, 3 * i, 3 , 3) = triangleCoordinateMatrix(triangleVector.at(i));
    }
    for(int i = 0; i < quadrilateralVector.length(); i++)
    {
        coordinates.block(0, 3 * triangleVector.length() + 4 * i, 3 , 4) = quadrilateralCoordinateMatrix(quadrilateralVector.at(i));
    }
    return coordinates.rowwise().mean();
}

double global::calculateMaxRelativeDistance(const QVector<VectorTriangle> &triangleVector, const Eigen::Vector3d &refNode)
{
    if(triangleVector.length() == 0)
    {
        return 0;
    }
    Eigen::MatrixXd coordinates(3, 3 * triangleVector.length());
    for(int i = 0; i < triangleVector.length(); i++)
    {
        coordinates.block(0, 3 * i, 3 , 3) = triangleCoordinateMatrix(triangleVector.at(i));
    }
    coordinates.colwise() -= refNode;
    return coordinates.colwise().norm().maxCoeff();
}

Eigen::Vector3d global::calculateMaxNorm(const QVector<VectorTriangle> &triangleVector/*, const QVector<VectorQuadrilateral>& quadrilateralVector*/)
{
    Eigen::MatrixXd coordinates(3, 3 * triangleVector.length());
    for(int i = 0; i < triangleVector.length(); i++)
    {
        coordinates.block(0, 3 * i, 3 , 3) = triangleCoordinateMatrix(triangleVector.at(i));
    }
    return coordinates.rowwise().maxCoeff();
}

Eigen::Vector3d global::triangleCompWiseMaxNorm(const VectorTriangle &triangle)
{
    double x,y,z;

    x = maxOfThree(triangle.node1(0), triangle.node2(0), triangle.node3(0));

    y = maxOfThree(triangle.node1(1), triangle.node2(1), triangle.node3(1));

    z = maxOfThree(triangle.node1(2), triangle.node2(2), triangle.node3(2));

    Eigen::Vector3d returnPoint={x,y,z};
    return returnPoint;
}

Eigen::Vector3d global::triangleCompWiseMinNorm(const VectorTriangle &triangle)
{
    double x,y,z;

    x = minOfThree(triangle.node1(0), triangle.node2(0), triangle.node3(0));

    y = minOfThree(triangle.node1(1), triangle.node2(1), triangle.node3(1));

    z = minOfThree(triangle.node1(2), triangle.node2(2), triangle.node3(2));

    Eigen::Vector3d returnPoint={x,y,z};
    return returnPoint;
}



Eigen::Vector3d global::calculateMinNorm(const QVector<VectorTriangle> &triangleVector/*, const QVector<VectorQuadrilateral>& quadrilateralVector*/)
{
    Eigen::MatrixXd coordinates(3, 3 * triangleVector.length());
    for(int i = 0; i < triangleVector.length(); i++)
    {
        coordinates.block(0, 3 * i, 3 , 3) = triangleCoordinateMatrix(triangleVector.at(i));
    }
    return coordinates.rowwise().minCoeff();
}

quint64 global::indexOfMaxOfThree(double a, double b, double c)
{
    if(a>b)
    {
        if(a>c)
        {
            return 1;
        }
        else
        {
            return 3;
        }
    }
    else
    {
        if(b>c)
        {
            return 2;
        }
        else
        {
            return 3;
        }

    }
}

double global::maxOfThree(double a, double b, double c)
{
    if(a>b)
    {
        if(a>c)
        {
            return a;
        }
        else
        {
            return c;
        }
    }
    else
    {
        if(b>c)
        {
            return b;
        }
        else
        {
            return c;
        }
    }
}

double global::maxClDouble3(const  Eigen::Vector3d &node)
{
    double a=node(0);
    double b=node(1);
    double c=node(2);
    if(a>b)
    {
        if(a>c)
        {
            return a;
        }
        else
        {
            return c;
        }
    }
    else
    {
        if(b>c)
        {
            return b;
        }
        else
        {
            return c;
        }
    }
}

double global::minClDouble3(const Eigen::Vector3d &node)
{
    double a=node(0);
    double b=node(1);
    double c=node(2);
    if(a>b)
    {
        if(b>c)
        {
            return c;
        }
        else
        {
            return b;
        }
    }
    else
    {
        if(a>c)
        {
            return c;
        }
        else
        {
            return a;
        }
    }
}

double global::minOfThree(double a, double b, double c)
{
    if(a>b)
    {
        if(b>c)
        {
            return c;
        }
        else
        {
            return b;
        }
    }
    else
    {
        if(a>c)
        {
            return c;
        }
        else
        {
            return a;
        }
    }
}

double global::maxOfFour(double a, double b, double c,double d)
{
    if(a>b)
    {
        if(a>c)
        {
            if(a>d)
            {
                return 1;
            }
            else
            {
                return 4;
            }
        }
        else
        {
            if(c>d)
            {
                return 3;
            }
            else
            {
                return 4;
            }
        }
    }
    else
    {
        if(b>c)
        {
            if(b>d)
            {
                return 2;
            }
            else
            {
                return 4;
            }
        }
        else
        {
            if(c>d)
            {
                return 3;
            }
            else
            {
                return 4;
            }
        }

    }
}

int global::allIntsUnique(int a, int b, int c, int d)
{
    if(a==b || a==c || a==d || b==c || b==d || c==d)
    {
        return false;
    }
    else
    {
        return true;
    }
}


double global::lineLength(const Eigen::Vector3d &point1, const Eigen::Vector3d &point2)
{
    Eigen::Vector3d vector = point1 - point2;
    return vector.norm();
}

double global::quadrilateralAspectRatio(const VectorQuadrilateral &parentQuadrilateral)
{
    QVector3D node1(parentQuadrilateral.node1(0),parentQuadrilateral.node1(1),parentQuadrilateral.node1(2));
    QVector3D node2(parentQuadrilateral.node2(0),parentQuadrilateral.node2(1),parentQuadrilateral.node2(2));
    QVector3D node3(parentQuadrilateral.node3(0),parentQuadrilateral.node3(1),parentQuadrilateral.node3(2));
    QVector3D node4(parentQuadrilateral.node4(0),parentQuadrilateral.node4(1),parentQuadrilateral.node4(2));
    QVector3D midpoint=(1.0/4.0)*(node1+node2+node3+node4);
    QVector3D node1NormalVector=QVector3D::normal(node1,node2,node4);
    QVector3D node2NormalVector=QVector3D::normal(node2,node3,node1);
    QVector3D node3NormalVector=QVector3D::normal(node3,node4,node2);
    QVector3D node4NormalVector=QVector3D::normal(node4,node1,node3);
    QVector3D averageNormal=(1.0/4.0)*(node1NormalVector+node2NormalVector+node3NormalVector+node4NormalVector); //normalized Vector
    double distanceNode1=node1.distanceToPlane(midpoint,averageNormal); //distance of node 1 to plane through the average of the four nodes and perpendicular to the average of the normals at the nodes
    double distanceNode2=node1.distanceToPlane(midpoint,averageNormal);
    double distanceNode3=node1.distanceToPlane(midpoint,averageNormal);
    double distanceNode4=node1.distanceToPlane(midpoint,averageNormal);
    QVector3D projectedNode1=node1-distanceNode1*averageNormal;
    QVector3D projectedNode2=node2-distanceNode2*averageNormal;
    QVector3D projectedNode3=node3-distanceNode3*averageNormal;
    QVector3D projectedNode4=node4-distanceNode4*averageNormal;
    QVector3D projectedMidpoint1=(1.0/2.0)*(projectedNode1+projectedNode2);
    QVector3D projectedMidpoint2=(1.0/2.0)*(projectedNode2+projectedNode3);
    QVector3D projectedMidpoint3=(1.0/2.0)*(projectedNode3+projectedNode4);
    QVector3D projectedMidpoint4=(1.0/2.0)*(projectedNode4+projectedNode1);
    double rectangle1Side1=(projectedMidpoint1-projectedMidpoint3).length();
    double rectangle1Side2=2.0*(projectedMidpoint2.distanceToLine(projectedMidpoint1,projectedMidpoint3-projectedMidpoint1));
    double maxSide1 =std::max(rectangle1Side1,rectangle1Side2);
    double minSide1 =std::min(rectangle1Side1,rectangle1Side2);
    double rectangle2Side1=(projectedMidpoint2-projectedMidpoint4).length();
    double rectangle2Side2=2.0*(projectedMidpoint3.distanceToLine(projectedMidpoint2,projectedMidpoint4-projectedMidpoint2));
    double maxSide2 =std::max(rectangle2Side1,rectangle2Side2);
    double minSide2 =std::min(rectangle2Side1,rectangle2Side2);
    double aspectRatio=std::max(maxSide1/minSide1, maxSide2/minSide2);

    return aspectRatio;
}

QVector<VectorTriangle> global::quadrilateralsToTriangles(const QVector<VectorQuadrilateral> &quadrilaterals)
{
    QVector<VectorTriangle> triangles;
    for(int i = 0; i < quadrilaterals.length(); i++)
    {
        VectorQuadrilateral parentQuadrilateral=quadrilaterals.at(i);

        double diagonal1 = lineLength(parentQuadrilateral.node1, parentQuadrilateral.node3);
        double diagonal2 = lineLength(parentQuadrilateral.node2, parentQuadrilateral.node4);
        VectorTriangle vectorTriangle1;
        VectorTriangle vectorTriangle2;
        if(diagonal1 > diagonal2) //split quadrilateral by shorter diagonal
        {
            vectorTriangle1 = VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,parentQuadrilateral.node2,parentQuadrilateral.node4);
            vectorTriangle2 = VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node2,parentQuadrilateral.node3,parentQuadrilateral.node4);
        }
        else
        {
            vectorTriangle1 = VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node4,parentQuadrilateral.node1,parentQuadrilateral.node3);
            vectorTriangle2 = VectorTriangle(parentQuadrilateral.elementIndex,parentQuadrilateral.node1,parentQuadrilateral.node2,parentQuadrilateral.node3);
        }
        vectorTriangle1.robinBoundaryCondition = parentQuadrilateral.robinBoundaryCondition;
        vectorTriangle2.robinBoundaryCondition = parentQuadrilateral.robinBoundaryCondition;
        triangles.push_back(vectorTriangle1);
        triangles.push_back(vectorTriangle2);
    }
    return triangles;

}

QPair<Eigen::Vector3d, double> global::centerAndVolumeOfMassOfTriangleObject(const QVector<VectorTriangle> &triangleVector)
{
       double totalVolume = 0, currentVolume;
       double xCenter = 0, yCenter = 0, zCenter = 0;
       Eigen::Vector3d center;
       center={0,0,0};

        for (int i = 0; i < triangleVector.length(); i++)
        {
            VectorTriangle currentTriangle=triangleVector.at(i);
            currentVolume = (currentTriangle.node1(0)*currentTriangle.node2(1)*currentTriangle.node3(2) - currentTriangle.node1(0)*currentTriangle.node3(1)*currentTriangle.node2(2) - currentTriangle.node2(0)*currentTriangle.node1(1)*currentTriangle.node3(2) + currentTriangle.node2(0)*currentTriangle.node3(1)*currentTriangle.node1(2) + currentTriangle.node3(0)*currentTriangle.node1(1)*currentTriangle.node2(2) - currentTriangle.node3(0)*currentTriangle.node2(1)*currentTriangle.node1(2)) / 6;
            totalVolume += currentVolume;
            xCenter += ((currentTriangle.node1(0) + currentTriangle.node2(0) + currentTriangle.node3(0)) / 4) * currentVolume;
            yCenter += ((currentTriangle.node1(1) + currentTriangle.node2(1) + currentTriangle.node3(1)) / 4) * currentVolume;
            zCenter += ((currentTriangle.node1(2) + currentTriangle.node2(2) + currentTriangle.node3(2)) / 4) * currentVolume;
        }

        std::cout << std::endl << "Total Volume = " << totalVolume << std::endl;
        std::cout << std::endl << "X center = " << xCenter/totalVolume << std::endl;
        std::cout << std::endl << "Y center = " << yCenter/totalVolume << std::endl;
        std::cout << std::endl << "Z center = " << zCenter/totalVolume << std::endl;
        if(totalVolume != 0)
        {
            center = {xCenter/totalVolume, yCenter/totalVolume, zCenter/totalVolume};
        }
        else
        {
           center = {0,0,0};
        }
       return qMakePair(center, totalVolume);
}

Eigen::Vector3d global::normalVectorOfTriangle(const VectorTriangle &triangle)
{   //node ordering is counter clockwise for triangle front side
    Eigen::Vector3d planeVector1 = triangle.node2 - triangle.node1;
    Eigen::Vector3d planeVector2 = triangle.node3 - triangle.node1;
    return planeVector1.cross(planeVector2).normalized();
}

Eigen::Vector3d global::midpointOfTriangle(const VectorTriangle &triangle)
{
    return Eigen::Vector3d(triangle.node1 + triangle.node2 + triangle.node3) / 3.0;
}

//Eigen::Vector3d global::projectPointOnTriangle(const VectorTriangle &triangle, const Eigen::Vector3d &point, bool &inInterior, bool &onNode1, bool &onNode2, bool &onNode3, bool &onSide21, bool &onSide32, bool &onSide13)
//{
//    Eigen::Vector3d side1 = (triangle.node2-triangle.node1)/*.normalized()*/;
//    Eigen::Vector3d side2 = (triangle.node3-triangle.node1)/*.normalized()*/;
//    Eigen::Vector3d normal = side1.cross(side2).normalized();

//    Eigen::Matrix3d system;
//    system << side1, side2, normal; // find the projection of the point onto the plane of the triangle in the affine coordinates of the triangle sides

//    Eigen::Vector3d rightHandSide = point - triangle.node1;

//    Eigen::Vector3d projectionOnPlaneParameters = system.partialPivLu().solve(rightHandSide);
//    double affineCoordinate0 = projectionOnPlaneParameters(0);
//    double affineCoordinate1 = projectionOnPlaneParameters(1);
//    Eigen::Vector3d projectionOnPlane = triangle.node1 + affineCoordinate0 * side1 + affineCoordinate1 * side2;

//    if(affineCoordinate0 > 0 && affineCoordinate1 > 0 && affineCoordinate0 + affineCoordinate1 < 1) // projection point is in the interior of the triangle;
//    {
//        inInterior = true;
//        onNode1 = false;
//        onNode2 = false;
//        onNode3 = false;
//        onSide21 = false;
//        onSide32 = false;
//        onSide13 = false;
//        return projectionOnPlane;
//    }
//    else // the projected point is in the exterior of the triangle
//    {
//        Eigen::Vector3d lineToNode1 = triangle.node1 - projectionOnPlane;
//        Eigen::Vector3d lineToNode2 = triangle.node2 - projectionOnPlane;
//        Eigen::Vector3d lineToNode3 = triangle.node3 - projectionOnPlane;
//        double distToNode1 = lineToNode1.norm();
//        double distToNode2 = lineToNode2.norm();
//        double distToNode3 = lineToNode3.norm();
//        int nearestNode = indexOfMaxOfThree(distToNode1, distToNode2, distToNode3); // find the nearest node to the projected point
//    }
//    return Eigen::Vector3d(triangle.node1 + triangle.node2 + triangle.node3) / 3.0;
//}

std::pair<Eigen::Vector3d, global::triangleProjectRegion> global::projectPointOnTriangleFaster(const VectorTriangle &triangle, const Eigen::Vector3d &point)
{
    global::triangleProjectRegion projectionRegion;
    Eigen::Vector3d diff = triangle.node1 - point;
    Eigen::Vector3d edge0 = triangle.node2 - triangle.node1;
    Eigen::Vector3d edge1 = triangle.node3 - triangle.node1;
    double a00 = edge0.dot(edge0);
    double a01 = edge0.dot(edge1);
    double a11 = edge1.dot(edge1);
    double b0 = diff.dot(edge0);
    double b1 = diff.dot(edge1);
    double det = std::max(a00 * a11 - a01 * a01, 0.0);
    double s = a01 * b1 - a11 * b0;
    double t = a01 * b0 - a00 * b1;

    if (s + t <= det)
    {
       if (s < 0.0)
       {
           if (t < 0.0)  // region 4
           {
               if (b0 < 0.0)
               {
                   t = 0.0;
                   if (-b0 >= a00)
                   {
                       s = 1.0;
//                       onNode2 = true;
                   }
                   else
                   {
                       s = -b0 / a00;
//                       onSide21 = true;
                   }
               }
               else
               {
                   s = 0.0;
                   if (b1 >= 0.0)
                   {
                       t = 0.0;
//                       onNode1 = true;
                   }
                   else if (-b1 >= a11)
                   {
                       t = 1.0;
//                       onNode3 = true;
                   }
                   else
                   {
                       t = -b1 / a11;
//                       onSide13 = true;
                   }
               }
           }
           else  // region 3
           {
               s = 0.0;
               if (b1 >= 0.0)
               {
                   t = 0.0;
//                   onNode1 = true;
               }
               else if (-b1 >= a11)
               {
                   t = 1.0;
//                   onNode3 = true;
               }
               else
               {
                   t = -b1 / a11;
//                   onSide13 = true;
               }
           }
       }
       else if (t < 0.0)  // region 5
       {
           t = 0.0;
           if (b0 >= 0.0)
           {
               s = 0.0;
//               onNode1 = true;
           }
           else if (-b0 >= a00)
           {
               s = 1.0;
//               onNode2 = true;
           }
           else
           {
               s = -b0 / a00;
//               onSide21 = true;
           }
       }
       else  // region 0
       {
           // minimum at interior point
           s /= det;
           t /= det;
//           inInterior = true;
       }
    }
    else
    {
       double tmp0, tmp1, numer, denom;
       if (s < 0.0)  // region 2
       {
           tmp0 = a01 + b0;
           tmp1 = a11 + b1;
           if (tmp1 > tmp0)
           {
               numer = tmp1 - tmp0;
               denom = a00 - 2.0 * a01 + a11;
               if (numer >= denom)
               {
                   s = 1.0;
                   t = 0.0;
//                   onNode2 = true;
               }
               else
               {
                   s = numer / denom;
                   t = 1.0 - s;
//                   onSide32 = true;
               }
           }
           else
           {
               s = 0.0;
               if (tmp1 <= 0.0)
               {
                   t = 1.0;
//                   onNode3 = true;
               }
               else if (b1 >= 0.0)
               {
                   t = 0.0;
//                   onNode1 = true;
               }
               else
               {
                   t = -b1 / a11;
//                   onSide13 = true;
               }
           }
       }
       else if (t < 0.0)  // region 6
       {
           tmp0 = a01 + b1;
           tmp1 = a00 + b0;
           if (tmp1 > tmp0)
           {
               numer = tmp1 - tmp0;
               denom = a00 - 2.0 * a01 + a11;
               if (numer >= denom)
               {
                   t = 1.0;
                   s = 0.0;
//                   onNode3 = true;
               }
               else
               {
                   t = numer / denom;
                   s = 1.0 - t;
//                   onSide32 = true;
               }
           }
           else
           {
               t = 0.0;
               if (tmp1 <= 0.0)
               {
                   s = 1.0;
//                   onNode2 = true;
               }
               else if (b0 >= 0.0)
               {
                   s = 0.0;
//                   onNode1 = true;
               }
               else
               {
                   s = -b0 / a00;
//                   onSide21 = true;
               }
           }
       }
       else  // region 1
       {
           numer = a11 + b1 - a01 - b0;
           if (numer <= 0.0)
           {
               s = 0.0;
               t = 1.0;
//               onNode3 = true;
           }
           else
           {
               denom = a00 - 2.0 * a01 + a11;
               if (numer >= denom)
               {
                   s = 1.0;
                   t = 0.0;
//                   onNode2 = true;
               }
               else
               {
                   s = numer / denom;
                   t = 1.0 - s;
//                   onSide32 = true;
               }
           }
       }
    }

    double barycentricCoordinateNode1 = 1.0 - s - t;
    double barycentricCoordinateNode2 = s;
    double barycentricCoordinateNode3 = t;

    int  numberOfSignificantPoints = (barycentricCoordinateNode1 > epsilon) + (barycentricCoordinateNode2 > epsilon) + (barycentricCoordinateNode3 > epsilon);

    if( numberOfSignificantPoints == 3)
    {
        projectionRegion = inInterior;
    }
    else if(numberOfSignificantPoints == 2)
    {
        if(barycentricCoordinateNode1 > epsilon && barycentricCoordinateNode2 > epsilon)
        {
            projectionRegion = onSide21;
            barycentricCoordinateNode1 += barycentricCoordinateNode3/2.0;
            barycentricCoordinateNode2 += barycentricCoordinateNode3/2.0;
            barycentricCoordinateNode3 = 0.0;
        }
        else if(barycentricCoordinateNode1 > epsilon && barycentricCoordinateNode3 > epsilon)
        {
            projectionRegion = onSide13;
            barycentricCoordinateNode1 += barycentricCoordinateNode2/2.0;
            barycentricCoordinateNode3 += barycentricCoordinateNode2/2.0;
            barycentricCoordinateNode2 = 0.0;
        }
        else
        {
            projectionRegion = onSide32;
            barycentricCoordinateNode2 += barycentricCoordinateNode1/2.0;
            barycentricCoordinateNode3 += barycentricCoordinateNode1/2.0;
            barycentricCoordinateNode1 = 0.0;
        }
    }
    else
    {
        if(barycentricCoordinateNode1 > epsilon)
        {
            projectionRegion = onNode1;
            barycentricCoordinateNode1 = 1.0;
            barycentricCoordinateNode2 = 0.0;
            barycentricCoordinateNode3 = 0.0;
        }
        else if(barycentricCoordinateNode2 > epsilon)
        {
            projectionRegion = onNode2;
            barycentricCoordinateNode1 = 0.0;
            barycentricCoordinateNode2 = 1.0;
            barycentricCoordinateNode3 = 0.0;
        }
        else
        {
            projectionRegion = onNode3;
            barycentricCoordinateNode1 = 0.0;
            barycentricCoordinateNode2 = 0.0;
            barycentricCoordinateNode3 = 1.0;
        }
    }
    return std::make_pair(barycentricCoordinateNode1 * triangle.node1 + barycentricCoordinateNode2 * triangle.node2 + barycentricCoordinateNode3 * triangle.node3, projectionRegion); // return the projection point and the region where the point lies on the triangle
}

//void global::testProjectOnTriangle(const Eigen::Vector3d &point, const VectorTriangle &testTriangle)
//{
////    Eigen::Vector3d node1(0,0,0);
////    Eigen::Vector3d node2(1,0,0);
////    Eigen::Vector3d node3(0,1,0);

////    VectorTriangle testTriangle(node1, node2, node3);
////    VectorTriangle testTriangle({0,0,0}, {1,0,0}, {0,1,0});

//    bool inInterior = false;
//    bool onNode1 = false;
//    bool onNode2 = false;
//    bool onNode3 = false;
//    bool onSide21 = false;
//    bool onSide32 = false;
//    bool onSide13 = false;

//    Eigen::Vector3d projectionOFPoint = projectPointOnTriangleFaster(testTriangle, point, inInterior, onNode1, onNode2, onNode3, onSide21, onSide32, onSide13);
//    std::cerr << "point: " << point.transpose() << std::endl;
//    std::cerr << "projectionOFPoint: " << projectionOFPoint.transpose() << std::endl;
//    if(inInterior)
//    {
//        std::cerr << "inInterior" << std::endl;
//    }
//    else if(onNode1)
//    {
//        std::cerr << "onNode1" << std::endl;
//    }
//    else if(onNode2)
//    {
//        std::cerr << "onNode2" << std::endl;
//    }
//    else if(onNode3)
//    {
//        std::cerr << "onNode3" << std::endl;
//    }
//    else if(onSide21)
//    {
//        std::cerr << "onSide21" << std::endl;
//    }
//    else if(onSide32)
//    {
//        std::cerr << "onSide32" << std::endl;
//    }
//    else if(onSide13)
//    {
//        std::cerr << "onSide13" << std::endl;
//    }
//    else
//    {
//        std::cerr << "No projection classification!" << std::endl;
//    }
//}

double global::areaOfTriangle(const VectorTriangle &triangle)
{
    Eigen::Vector3d vec1 = (triangle.node1 - triangle.node2);
    Eigen::Vector3d vec2 = (triangle.node1 - triangle.node3);
    return (vec1.cross(vec2)).norm() / 2.0;
}

Eigen::MatrixXcd global::LUDecompNoPivoting(const Eigen::MatrixXcd &A) //L*U = A
{
    Eigen::MatrixXcd LU = Eigen::MatrixXcd::Zero(A.rows(), A.cols());
    std::complex<double> sum = 0;
    long minIndex = std::min(A.rows(), A.cols());
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = i; j < A.cols(); j++)
        {
            sum = 0;
            for (int k = 0; k < i && k < minIndex; k++)
            {
                sum += LU(i, k) * LU(k, j);
            }
            LU(i, j) = A(i, j) - sum;
        }
        for (int j = i + 1; j < A.rows() && i < A.cols(); j++)
        {
            sum = 0;
            for (int k = 0; k < i && k < minIndex; k++)
            {
                sum += LU(j, k) * LU(k, i);
            }
            LU(j, i) = (1.0 / LU(i, i)) * (A(j, i) - sum);
        }
    }
    return LU;
}

Eigen::VectorXcd global::forwardSubstitution(const Eigen::MatrixXcd &lowerTriangular, Eigen::VectorXcd rightHandSide)
{
    long rows = lowerTriangular.rows();
    long cols = lowerTriangular.cols();
    if(rows != cols)
    {
        std::cerr << "Nonsquare matrix in forwardSubstitution call." << std::endl;
        return Eigen::VectorXcd::Zero(rows);
    }
    Eigen::VectorXcd solution(rows);

    for(long rowIndex = 0; rowIndex < rows; ++rowIndex)
    {
        solution(rowIndex) = rightHandSide(rowIndex) / lowerTriangular(rowIndex, rowIndex);
        const long segmentStart = rowIndex + 1;
        const long segmentLength = rows - segmentStart;
        rightHandSide.segment(segmentStart, segmentLength) -= solution(rowIndex) * lowerTriangular.block(segmentStart, rowIndex, segmentLength, 1);
    }
    return solution;
}

Eigen::VectorXcd global::backwardSubstitution(const Eigen::MatrixXcd &upperTriangular, Eigen::VectorXcd rightHandSide)
{
    long rows = upperTriangular.rows();
    long cols = upperTriangular.cols();
    if(rows != cols)
    {
        std::cerr << "Nonsquare matrix in backwardSubstitution call." << std::endl;
        return Eigen::VectorXcd::Zero(rows);
    }
    Eigen::VectorXcd solution(rows);

    for(long rowIndex = rows - 1; rowIndex > 0; --rowIndex)
    {
        solution(rowIndex) = rightHandSide(rowIndex) / upperTriangular(rowIndex, rowIndex);
        const long segmentLength = rowIndex;
        rightHandSide.head(segmentLength) -= solution(rowIndex) * upperTriangular.block(0, rowIndex, segmentLength, 1);
    }
    solution(0) = rightHandSide(0) / upperTriangular(0, 0);
    return solution;
}

void global::printMatrixToFile(const QString fileName, const Eigen::MatrixXcd &A)
{
    long rows = A.rows();
    long cols = A.cols();
    Eigen::IOFormat CommaInitFmt(Eigen::FullPrecision, Eigen::DontAlignCols, ", std::complex<double>", ", ", "std::complex<double>", "\n", "Eigen::MatrixXcd testMatrix(" + std::to_string(rows) + "," + std::to_string(cols) + "); \n testMatrix << ", ";");
    std::ofstream file(fileName.toStdString(), std::ofstream::out | std::ofstream::trunc);
    if (file.is_open())
    {
        file << A.format(CommaInitFmt) << std::endl;
    }
}
