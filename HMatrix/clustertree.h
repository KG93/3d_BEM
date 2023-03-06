#ifndef CLUSTERTREE_H
#define CLUSTERTREE_H

#include "boundaryelements.h"
#include "LinearElements/linearboundaryelements.h"
#include "cuboid.h"
#include "global.h"

#include <concepts>

// constrain the template GeometryContainer to the types BoundaryElements, LinearBoundaryElements or QVector<VectorTriangle>
template<typename GeometryContainer>
    concept isGeometry  = std::is_same<GeometryContainer, BoundaryElements>::value
                            || std::is_same<GeometryContainer, LinearBoundaryElements>::value
                            || std::is_same<GeometryContainer, QVector<VectorTriangle>>::value;

/**
* \class Cluster
*\brief A recursive container for triangle indices.
*
* The triangle indices correspond to matrix row or column indices. A cluster contains a subset of the entire triangle index set.
*/
class Cluster
{
public:
    Cluster();
    ~Cluster(){
    }
    Cluster(QVector<long> indices, Cuboid minCuboid){
        this->indices = indices;
        this-> minCuboid = minCuboid;
    }
    Cluster(QVector<long> indices, Cuboid minCuboid, unsigned int clusterTreeDepth){
        this->indices = indices;
        this-> minCuboid = minCuboid;
        this-> clusterTreeDepth = clusterTreeDepth;
    }
    long size() const {return indices.size();} /*!< \brief Get the size of the cluster. */
    long startIndex() const {return indices.first();}; /*!< \brief Get the first index of the cluster. */

    void getRandomElementFromEachLeafNode(QVector<long> &rowIndices, const long maxNumberOfElements = 0);

    Cluster* returnCopy(Cluster* father); // returns a pointer to a copy of the entire subtree of this

    QVector<long> indices;
    unsigned int clusterTreeDepth;
    bool isRoot = false;
    bool isLeaf = false;

    /**
    * \brief The minimal axis aligned cuboid that contains the triangles (that are referred to by the triangle indices).
    */
    Cuboid minCuboid;

    /**
    * \brief The superset (supercluster) of the triangle indices.
    */
    Cluster* father = nullptr;

    /**
    * \brief The clusters son1 and son2 are a bipartition of the cluster. The cluster is split along its minCuboids' maximum elongation.
    */
    Cluster* son1 = nullptr;
    Cluster* son2 = nullptr;
};

/**
* \class ClusterTree
* \brief A clustertree contains increasingly (by tree depth) fine partitions of the triangle geometry.
*/
class ClusterTree
{
public:
    ClusterTree();

    ClusterTree(BoundaryElements* boundaryElements)
    {
        if(boundaryElements != nullptr)
        {
            boundaryElements->calculateTriangleMidPoints();
            assembleClustertree<BoundaryElements>(boundaryElements);
            boundaryElements->calculateCollocationPoints();
            boundaryElements->calculateTrianglesArea();
        }
    }

    ClusterTree(LinearBoundaryElements* linearBoundaryElements)
    {
        if(linearBoundaryElements != nullptr)
        {
            assembleClustertree<LinearBoundaryElements>(linearBoundaryElements);
        }
    }

    ClusterTree(QVector<VectorTriangle>* triangles)
    {
        if(triangles != nullptr)
        {
            for(int i=0; i<triangles->size(); i++)
            {
                (*triangles)[i].triangleMidpoint = global::midpointOfTriangle(triangles->at(i));
            }
            assembleClustertree<QVector<VectorTriangle>>(triangles);
        }
    }

    /**
     *\brief Create a copy of a cluster tree.
    */
    ClusterTree(ClusterTree &originalTree);
    ~ClusterTree()
    {
    }

//    bool isInCuboid(VectorTriangle triangle, Cuboid cuboid);
//    bool isInCuboid(LinearTriangle triangle, QVector<LinearTriangleNode> &nodes, Cuboid cuboid);
//    bool isInCuboid(LinearTriangleNode node, Cuboid cuboid);

    long getNumberOfElements(BoundaryElements* boundaryElements){ return boundaryElements->triangles.size(); }
    long getNumberOfElements(LinearBoundaryElements* linearBoundaryElements){ return linearBoundaryElements->nodes.size(); }
    long getNumberOfElements(QVector<VectorTriangle>* triangles){ return triangles->size(); }

    /**
    * \brief Return a vector of all clusters.
    */
    QVector<Cluster*> getClusters();

    /**
    * \brief Get the root cluster of the cluster tree.
    */
    Cluster* getRootCluster();

//    static const unsigned int triangleThreshold = 26;
    static constexpr unsigned int triangleThreshold = 50;

    unsigned int maxTreeDepth = 0;

    /**
    * \brief Delete all tree nodes.
    */
    void clear();

    /**
    * \brief Recursive node deletion.
    */
    void recursiveClear(Cluster* cluster);

    /**
    * \brief Recalculate the minimal cuboids for an existing clustertree. The function is used for the calculation of reflection matrices, where the cluster indices have to correspond to some preexisting partition, but the geometry hos been changed (reflected) and the mincuboids have to be updated accordingly.
    */
    template<class GeometryContainer>  requires isGeometry<GeometryContainer>
    void updateMinCuboids(GeometryContainer* container);

private:
    /**
     * \brief Build up the tree from the root node.
    */
    template<class GeometryContainer> requires isGeometry<GeometryContainer>
    void assembleClustertree(GeometryContainer* container);

    template<class GeometryContainer> requires isGeometry<GeometryContainer>
    void splitCluster( Cluster &clusterToSplit, GeometryContainer* container);

    template<class GeometryContainer>  requires isGeometry<GeometryContainer>
    void createRootCluster(GeometryContainer* container);

    Cuboid minimalCuboidForTriangles(QVector<VectorTriangle> triangles) const; /*!< Calculate the minimal axis aligned cuboid for a set of triangles. */
    Cuboid minimalCuboidForIndices(const QVector<long> &triangleIndices, const BoundaryElements *boundaryElements) const; /*!< Calculate the minimal axis aligned cuboid for a set of triangles. */
    Cuboid minimalCuboidForIndices(const QVector<long> &triangleIndices, const LinearBoundaryElements* linearBoundaryElements) const;
    Cuboid minimalCuboidForIndices(const QVector<long> &triangleIndices, const QVector<VectorTriangle>* triangles) const;

    Cuboid minimalCuboidForTriangleIndexes(QVector<long> triangleIndexes) const;
    Cuboid minimalCuboidForTriangleIndexesMidpointBased(const QVector<long> &triangleIndexes) const; /*!< Calculate the minimal axis aligned cuboid for the minpoints of a set of triangles. */
    Cuboid minimalCuboidForTriangleIndicesMidpointBased2(const QVector<long> &triangleIndexes) const;

    std::tuple<Cuboid, Cuboid> splitQuboid(Cuboid quboid) const; /*!< Split a cuboid in two halves along its axis of longest elongation. */

    bool isInCuboid(const long index, const BoundaryElements* boundaryElements, const Cuboid cuboid);
    bool isInCuboid(const long index, const LinearBoundaryElements* linearBoundaryElements, const Cuboid cuboid);
    bool isInCuboid(const long index, const QVector<VectorTriangle>* triangles, const Cuboid cuboid);

//    bool cuboidContainsTriangle(VectorTriangle triangle, Cuboid cuboid); /*!< Check whether a triangle (midpoint) is contained in a cuboid. */
//    bool cuboidContainsCuboid(Cuboid cuboid, Cuboid quboidContainer); /*!< Check whether a cuboid is contained in a second cuboid. */
    void assembleDepthFirstClusterVector(); /*!< Construct a vector of pointers to all clusters in the tree in a depth-first fashion. */
    void assembleClusterVectorRecursion(Cluster* cluster);

    template<class GeometryContainer>  requires isGeometry<GeometryContainer>
    void orderGeometryContiguously(GeometryContainer* container);

    void recursiveReorder(Cluster* cluster, long &globalTriangleIndex, BoundaryElements* boundaryElements, const BoundaryElements &boundaryElementsCopy); /*!< Reindex the triangles in 'BoundaryElements* boundaryElements' so that each cluster contains a contiguous set of triangle indices; so that a cluster corresponds to a contiguous vector segment and that cluster x cluster corresponds to a matrix block. */
    void recursiveReorder(Cluster* cluster, long &globalTriangleIndex, QVector<VectorTriangle>* triangles, const QVector<VectorTriangle> &trianglesCopy); /*!< Reindex the triangles in 'QVector<VectorTriangle>* triangles' so that each cluster contains a contiguous set of triangle indices; so that a cluster corresponds to a contiguous vector segment and that cluster x cluster corresponds to a matrix block. */
    void recursiveReorder(Cluster* cluster, long &globalNodeIndex, LinearBoundaryElements* linearBoundaryElements, const LinearBoundaryElements &boundaryElementsCopy); /*!< Reindex the triangles in 'LinearBoundaryElements* linearBoundaryElements' so that each cluster contains a contiguous set of triangle indices; so that a cluster corresponds to a contiguous vector segment and that cluster x cluster corresponds to a matrix block. */


    Cluster* rootCluster = nullptr;
    QVector<Cluster*> clusters;
    void reindexLinearTriangleNodes(LinearTriangle &triangle, const long oldIndex, const long  newIndex);
};

template<class GeometryContainer> requires isGeometry<GeometryContainer>
void ClusterTree::assembleClustertree(GeometryContainer* container)
{
    createRootCluster<GeometryContainer>(container);
    splitCluster<GeometryContainer>(*rootCluster, container);
    orderGeometryContiguously<GeometryContainer>(container);
    assembleDepthFirstClusterVector();
}

template<class GeometryContainer> requires isGeometry<GeometryContainer>
void ClusterTree::updateMinCuboids(GeometryContainer* container)
{
    this->assembleDepthFirstClusterVector();
    for(long i = 0; i < clusters.length(); i++)
    {
//        long startIndex = clusters.at(i)->indices.first();
//        long endIndex = clusters.at(i)->indices.last();
        clusters[i]->minCuboid = minimalCuboidForIndices(clusters.at(i)->indices, container);
    }
}

template<class GeometryContainer> requires isGeometry<GeometryContainer>
void ClusterTree::createRootCluster(GeometryContainer* container)
{
    if(rootCluster != nullptr)
    {
        recursiveClear(rootCluster);
    }
    clusters.clear();
    if(container == nullptr)
    {
        std::cout << "No boundary elements were passed to clusterTree before create rootCluster call!" << std::endl;
        return;
    }
    long numberOfTriangles = getNumberOfElements(container);
    if(numberOfTriangles == 0)
    {
        std::cout<<"numberOfTriangles == 0 in create rootCluster call!"<<std::endl;
        return;
    }
    QVector<long> triangleIndices = global::createContiguousIndexVector(0, numberOfTriangles-1);
    rootCluster = new Cluster(triangleIndices, minimalCuboidForIndices(triangleIndices, container), 0);
    rootCluster->isRoot = true;
}

template<class GeometryContainer> requires isGeometry<GeometryContainer>
void ClusterTree::splitCluster(Cluster &clusterToSplit, GeometryContainer* container)
{
    if(clusterToSplit.indices.length() < triangleThreshold || clusterToSplit.indices.length() <= 1)
    {
        clusterToSplit.isLeaf = true;
//        std::cout<<"cluster is leaf. Unneccessary splitCluster() call."<< std::endl;
        return;
    }
    Eigen::Vector3d cuboidSideLengths = clusterToSplit.minCuboid.maxPoint - clusterToSplit.minCuboid.minPoint;
    unsigned int longestEdgeIndex = 0;
    cuboidSideLengths.cwiseAbs().maxCoeff(&longestEdgeIndex);

    Cuboid son1Cuboid = clusterToSplit.minCuboid;
    son1Cuboid.minPoint(longestEdgeIndex) = clusterToSplit.minCuboid.minPoint(longestEdgeIndex) + 0.5 * cuboidSideLengths(longestEdgeIndex);
    Cuboid son2Cuboid = clusterToSplit.minCuboid;
    son2Cuboid.maxPoint(longestEdgeIndex) = clusterToSplit.minCuboid.minPoint(longestEdgeIndex) + 0.5 * cuboidSideLengths(longestEdgeIndex);

    QVector<long> son1TriangleIndexes;
    QVector<long> son2TriangleIndexes;
    for(int i = 0; i < clusterToSplit.indices.length(); i++)
    {
        long triangleIndex = clusterToSplit.indices.at(i);
        bool triangleInSon1Cuboid = isInCuboid(triangleIndex, container, son1Cuboid);

        if(triangleInSon1Cuboid)
        {
            son1TriangleIndexes.append(triangleIndex);
        }
        else
        {
            son2TriangleIndexes.append(triangleIndex);
        }
    }
    clusterToSplit.son1 = new Cluster(son1TriangleIndexes, minimalCuboidForIndices(son1TriangleIndexes, container), clusterToSplit.clusterTreeDepth+1);
    clusterToSplit.son1->father = &clusterToSplit;
    splitCluster(*clusterToSplit.son1, container);

    clusterToSplit.son2 = new Cluster(son2TriangleIndexes, minimalCuboidForIndices(son2TriangleIndexes, container), clusterToSplit.clusterTreeDepth+1);
    clusterToSplit.son2->father = &clusterToSplit;
    splitCluster(*clusterToSplit.son2, container);

    if(maxTreeDepth < clusterToSplit.clusterTreeDepth + 1)
    {
        maxTreeDepth = clusterToSplit.clusterTreeDepth + 1;
    }
    return;
}

template<class GeometryContainer> requires isGeometry<GeometryContainer>
void ClusterTree::orderGeometryContiguously(GeometryContainer* container)
{
    long globalIndex = 0;
    if(rootCluster == nullptr)
    {
        std::cerr << "orderTrianglesContiguously called on nullptr rootCluster." << std::endl;
        return;
    }
    GeometryContainer geometryContainerCopy;
    if(container != nullptr)
    {
        geometryContainerCopy = *container;
        recursiveReorder(rootCluster, globalIndex, container, geometryContainerCopy);
    }
}
#endif // CLUSTERTREE_H
