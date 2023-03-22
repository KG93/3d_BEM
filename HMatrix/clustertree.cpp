#include "clustertree.h"

Cluster::Cluster()
{

}

void Cluster::getRandomElementFromEachLeafNode(QVector<long> &rowIndices, const long maxNumberOfElements)
{
    if(maxNumberOfElements > 0 && rowIndices.size() > maxNumberOfElements)
    {
        return;
    }
    if(isLeaf)
    {
        long randomIndex =  QRandomGenerator::system()->bounded( (qint32)(indices.last() - indices.first() + 1)) + indices.first();
        rowIndices.append(randomIndex);
    }
    else
    {
        son1->getRandomElementFromEachLeafNode(rowIndices, maxNumberOfElements);
        son2->getRandomElementFromEachLeafNode(rowIndices, maxNumberOfElements);
    }
}

std::unique_ptr<Cluster> Cluster::returnCopy(Cluster* father) // returns a pointer to a copy of the entire subtree of this
{
    std::unique_ptr<Cluster> returnCluster = std::make_unique<Cluster>();

    returnCluster->indices = this->indices;
    returnCluster->clusterTreeDepth = this->clusterTreeDepth;
    returnCluster->minCuboid = this->minCuboid;
    returnCluster->isLeaf = this->isLeaf;


    if(father != nullptr)
    {
        returnCluster->father = father;
    }
    else
    {
        returnCluster->isRoot = true;
    }
    if(!this->isLeaf)
    {
        if(this->son1 != nullptr)
        {
            returnCluster->son1 = this->son1->returnCopy(returnCluster.get());
        }
        else
        {
            std::cerr << "!this->isLeaf && this->son1 == nullptr in returnCopy() call!" << std::endl;
        }
        if(this->son2 != nullptr)
        {
            returnCluster->son2 = this->son2->returnCopy(returnCluster.get());
        }
        else
        {
            std::cerr << "!this->isLeaf && this->son2 == nullptr in returnCopy() call!" << std::endl;
        }
    }
    else
    {
        returnCluster->isLeaf = true;
    }
    return returnCluster;
}

///////Clustertree
///
ClusterTree::ClusterTree()
{

}

ClusterTree::ClusterTree(ClusterTree &originalTree)
{
    this->rootCluster = originalTree.getRootCluster()->returnCopy(nullptr);
    this->assembleDepthFirstClusterVector();
}

long ClusterTree::size() const
{
    if(rootCluster != nullptr)
    {
        return rootCluster->size();
    }
    return 0;
}

long ClusterTree::startIndex() const
{
    if(rootCluster != nullptr)
    {
        return rootCluster->startIndex();
    }
    return 0;
}

Cuboid ClusterTree::minimalCuboidForTriangles(QVector<VectorTriangle> triangles) const
{
    Cuboid cuboid;
    for(int i=0; i<triangles.size(); i++)
    {
        cuboid |= triangles.at(i).node1;
        cuboid |= triangles.at(i).node2;
        cuboid |= triangles.at(i).node3;
    }
    return cuboid;
//    return Cuboid(global::calculateMinNorm(triangles), global::calculateMaxNorm(triangles));
}

Cuboid ClusterTree::minimalCuboidForIndices(const QVector<long> &triangleIndexes, const BoundaryElements* boundaryElements) const
{
    Cuboid cuboid;
    for(int i=0; i<triangleIndexes.size(); i++)
    {
        long tmpIndex = triangleIndexes.at(i);
        cuboid |= boundaryElements->triangles.at(tmpIndex).triangleMidpoint;
    }
    return cuboid;
}

Cuboid ClusterTree::minimalCuboidForIndices(const QVector<long> &triangleIndexes, const LinearBoundaryElements* linearBoundaryElements) const
{
    Cuboid cuboid;
    for(int i=0; i<triangleIndexes.size(); i++)
    {
        long tmpIndex = triangleIndexes.at(i);
        cuboid |= linearBoundaryElements->nodes.at(tmpIndex).coordinates;
    }
    return cuboid;
}

Cuboid ClusterTree::minimalCuboidForIndices(const QVector<long> &triangleIndices, const QVector<VectorTriangle>* triangles) const
{
    Cuboid cuboid;
    for(int i=0; i<triangleIndices.size(); i++)
    {
        long tmpIndex = triangleIndices.at(i);
        cuboid |= triangles->at(tmpIndex).triangleMidpoint;
    }
    return cuboid;
}

Cuboid ClusterTree::minimalCuboidForIndices(const QVector<long> &pointIndices, const QVector<Eigen::Vector3d>* points) const
{
    Cuboid cuboid;
    for(int i=0; i<pointIndices.size(); i++)
    {
        long tmpIndex = pointIndices.at(i);
        cuboid |= points->at(tmpIndex);
    }
    return cuboid;
}

std::tuple<Cuboid, Cuboid> ClusterTree::splitQuboid(Cuboid quboid) const
{
    Eigen::Vector3d cuboidSideLengths=quboid.maxPoint-quboid.minPoint;
    uint maxIndex;
    cuboidSideLengths.cwiseAbs().maxCoeff(&maxIndex);

    Cuboid returnQuboid1 = quboid;
    returnQuboid1.minPoint(maxIndex) = quboid.minPoint(maxIndex) + cuboidSideLengths(maxIndex)/2;
    Cuboid returnQuboid2 = quboid;
    returnQuboid2.maxPoint(maxIndex) = quboid.minPoint(maxIndex) + cuboidSideLengths(maxIndex)/2;
    return std::make_tuple(returnQuboid1, returnQuboid2);
}

bool ClusterTree::isInCuboid(const long index, const BoundaryElements* boundaryElements, const Cuboid cuboid)
{
    if( ((cuboid.maxPoint - boundaryElements->triangles.at(index).triangleMidpoint).array() >= 0 ).all()  &&  ((cuboid.minPoint - boundaryElements->triangles.at(index).triangleMidpoint).array() <= 0 ).all() )
    {
       return true;
    }
    else
    {
        return false;
    }
}

bool ClusterTree::isInCuboid(const long index, const LinearBoundaryElements* linearBoundaryElements, const Cuboid cuboid)
{
    if( ((cuboid.maxPoint - linearBoundaryElements->nodes.at(index).coordinates).array() >= 0 ).all()  &&  ((cuboid.minPoint - linearBoundaryElements->nodes.at(index).coordinates).array() <= 0 ).all() )
    {
       return true;
    }
    else
    {
        return false;
    }
}

bool ClusterTree::isInCuboid(const long index, const QVector<VectorTriangle>* triangles, const Cuboid cuboid)
{
    if( ((cuboid.maxPoint - triangles->at(index).triangleMidpoint).array() >= 0 ).all()  &&  ((cuboid.minPoint - triangles->at(index).triangleMidpoint).array() <= 0 ).all() )
    {
       return true;
    }
    else
    {
        return false;
    }
}

bool ClusterTree::isInCuboid(const long index, const QVector<Eigen::Vector3d>* points, const Cuboid cuboid)
{
    if( ((cuboid.maxPoint - points->at(index)).array() >= 0 ).all()  &&  ((cuboid.minPoint - points->at(index)).array() <= 0 ).all() )
    {
       return true;
    }
    else
    {
        return false;
    }
}

QVector<Cluster*> ClusterTree::getClusters()
{
    if(rootCluster == nullptr)
    {
        clusters.clear();
        std::cout << "cluster is nullpopinter" << std::endl;
    }
    return clusters;
}

void ClusterTree::assembleDepthFirstClusterVector()
{
    if(rootCluster == nullptr)
    {
        clusters.clear();
        std::cerr << "cluster is nullpopinter in assembleDepthFirstClusterVector() call" << std::endl;
        return;
    }
    clusters.append(rootCluster.get());
    assembleClusterVectorRecursion(*rootCluster);
}

void  ClusterTree::assembleClusterVectorRecursion(Cluster& clusterToAdd)
{
    if(clusterToAdd.isLeaf == true || clusterToAdd.son1 == nullptr || clusterToAdd.son2 == nullptr)
    {
        return;
    }
    else
    {
        clusters.append(clusterToAdd.son1.get());
        assembleClusterVectorRecursion(*clusterToAdd.son1);
        clusters.append(clusterToAdd.son2.get());
        assembleClusterVectorRecursion(*clusterToAdd.son2);
    }
}

Cluster* ClusterTree::getRootCluster()
{
    if(rootCluster == nullptr)
    {
        std::cerr<<"getRootCluster called on empty Clustertree."<<std::endl;
    }
    return rootCluster.get();
}

void ClusterTree::recursiveReorder(Cluster &cluster, long &globalTriangleIndex, BoundaryElements* boundaryElements, const BoundaryElements &boundaryElementsCopy)
{
    if(cluster.isLeaf)
    {
        long numberOfTriangles = cluster.indices.length();
        for(int i = 0; i < numberOfTriangles; i++)
        {
            boundaryElements->triangles[globalTriangleIndex] = boundaryElementsCopy.triangles.at(cluster.indices.at(i));
            cluster.indices[i] = globalTriangleIndex;
            globalTriangleIndex++;
        }
    }
    else
    {
        recursiveReorder(*cluster.son1, globalTriangleIndex, boundaryElements, boundaryElementsCopy);
        recursiveReorder(*cluster.son2, globalTriangleIndex, boundaryElements, boundaryElementsCopy);
        cluster.indices = global::createContiguousIndexVector(cluster.son1->indices.first(), cluster.son2->indices.last());
    }
}

void ClusterTree::recursiveReorder(Cluster &cluster, long &globalTriangleIndex, QVector<VectorTriangle>* triangles, const QVector<VectorTriangle> &trianglesCopy)
{
    if(cluster.isLeaf)
    {
        long numberOfTriangles = cluster.indices.length();
        for(int i = 0; i < numberOfTriangles; i++)
        {
            (*triangles)[globalTriangleIndex] = trianglesCopy.at(cluster.indices.at(i));
            cluster.indices[i] = globalTriangleIndex;
            globalTriangleIndex++;
        }
    }
    else
    {
        recursiveReorder(*cluster.son1, globalTriangleIndex, triangles, trianglesCopy);
        recursiveReorder(*cluster.son2, globalTriangleIndex, triangles, trianglesCopy);
        cluster.indices = global::createContiguousIndexVector(cluster.son1->indices.first(), cluster.son2->indices.last());
    }
}

void ClusterTree::recursiveReorder(Cluster &cluster, long &globalNodeIndex, LinearBoundaryElements* linearBoundaryElements, const LinearBoundaryElements &boundaryElementsCopy)
{
    if(cluster.isLeaf)
    {
        long numberOfNodes = cluster.indices.length();
        for(int i = 0; i < numberOfNodes; i++)
        {
            long oldNodeIndex = cluster.indices.at(i);
            LinearTriangleNode tmpNode = boundaryElementsCopy.nodes.at(oldNodeIndex);
            linearBoundaryElements->nodes[globalNodeIndex] = tmpNode;
            for(int triangleIndex = 0; triangleIndex<tmpNode.associatedTriangles.size(); triangleIndex++)
            {
                reindexLinearTriangleNodes(linearBoundaryElements->triangles[tmpNode.associatedTriangles.at(triangleIndex)], oldNodeIndex, globalNodeIndex);
            }
            cluster.indices[i] = globalNodeIndex;
            globalNodeIndex++;
        }
    }
    else
    {
        recursiveReorder(*cluster.son1, globalNodeIndex, linearBoundaryElements, boundaryElementsCopy);
        recursiveReorder(*cluster.son2, globalNodeIndex, linearBoundaryElements, boundaryElementsCopy);

        cluster.indices = global::createContiguousIndexVector(cluster.son1->indices.first(), cluster.son2->indices.last());
    }
}

void ClusterTree::recursiveReorder(Cluster &cluster, long &globalPointIndex, QVector<Eigen::Vector3d>* points, const QVector<Eigen::Vector3d> &pointsCopy)
{
    if(cluster.isLeaf)
    {
        long numberOfPoints = cluster.indices.length();
        for(int i = 0; i < numberOfPoints; i++)
        {
            (*points)[globalPointIndex] = pointsCopy.at(cluster.indices.at(i));
            cluster.indices[i] = globalPointIndex;
            globalPointIndex++;
        }
    }
    else
    {
        recursiveReorder(*cluster.son1, globalPointIndex, points, pointsCopy);
        recursiveReorder(*cluster.son2, globalPointIndex, points, pointsCopy);
        cluster.indices = global::createContiguousIndexVector(cluster.son1->indices.first(), cluster.son2->indices.last());
    }
}

void ClusterTree::reindexLinearTriangleNodes(LinearTriangle &triangle, const long oldIndex, const long  newIndex)
{
    if(triangle.node1Index == oldIndex)
    {
        triangle.node1Index = newIndex;
    }
    else if(triangle.node2Index == oldIndex)
    {
        triangle.node2Index = newIndex;
    }
    else if(triangle.node3Index == oldIndex)
    {
        triangle.node3Index = newIndex;
    }
    else
    {
        std::cerr << "reindexLinearTriangleNodes called on triangle, where oldIndex was not found." << std::endl;
    }
}

void ClusterTree::clear()
{
    clusters.clear();
    if(rootCluster != nullptr)
    {
        recursiveClear(*rootCluster);
        rootCluster = nullptr;
    }
}

void ClusterTree::recursiveClear(Cluster &cluster)
{
    if(cluster.son1 != nullptr)
    {
        recursiveClear(*cluster.son1);
        cluster.son1 = nullptr;
    }
    if(cluster.son2 != nullptr)
    {
        recursiveClear(*cluster.son2);
        cluster.son2 = nullptr;
    }
    cluster.indices.clear();
    cluster.indices.squeeze();
}
