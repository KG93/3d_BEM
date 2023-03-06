#ifndef NODES_H
#define NODES_H

#include "node.h"
#include "global.h"
#include "meshfunctions.h"
#include "globallogstrings.h"

#include <QMatrix4x4>
#include <QHash>

/**
* \class NodesSection
* \brief This class represents a section of nodes in a script file.
*/
class NodesSection
{
public:
    NodesSection();

    /**
    * \brief Construct a new NodesSection object with the given name.
    * \param[in] name The name of the NodesSection.
    */
    NodesSection(QString name){
        this->name = name;
        Scale = {1,1,1};
    }

    /**
    * \brief Construct a new NodesSection object with the given name, MeshFiles, and MeshFileAlias.
    * \param[in] name The name of the NodesSection.
    * \param[in] MeshFiles The list of MeshFiles for the NodesSection.
    * \param[in] MeshFileAlias The list of MeshFileAliases for the NodesSection.
    */
    NodesSection(QString name, const QStringList& MeshFiles, const QStringList& MeshFileAlias){
        this->name = name;
        meshFiles = MeshFiles;
        meshFileAlias = MeshFileAlias;
        Scale = {1,1,1};
    }

    /**
    * \brief Handle a line of script for the NodesSection.
    * \param[in] line The line of script to handle.
    * \return `true` if the line was handled successfully, `false` otherwise.
    */
    bool handleScriptLine(const QString& line);
//    bool checkIdentifiersAllUnique(QVector<Node> nodes);
    /**
    * \brief Handle a line of script for the NodesSection.
    * \param[in] line The line of script to handle.
    * \return `true` if the line was handled successfully, `false` otherwise.
    */
    bool checkIdentifiersAllUnique();

    /**
    * \brief Scale the nodes in the NodesSection by the given factor.
    */
    void scaleNodes();

    /**
    * \brief Get the non-unique node identifiers in the NodesSection.
    * \return A vector of non-unique node identifiers.
    */
    QVector<quint64> getNonUniqueIdentifiers(){return nonUniqueIdentifiers;}

    /**
    * \brief The name of the NodesSection.
    */
    QString name;
    Eigen::Vector3d Scale = Eigen::Vector3d(1,1,1);

    /**
    * \brief A shifting factor for the nodes in the NodesSection.
    */
    Eigen::Vector3d Shift = Eigen::Vector3d(0,0,0);

    /**
    * \brief Rotation values for rotating the nodes in the NodesSection around the x,y and z axis respectively.
    */
    Eigen::Vector3d Rotate = Eigen::Vector3d(0,0,0);

    /**
    * \brief The combined transformation matrix.
    */
    QMatrix4x4 transformationMatrix;

    /**
    * \brief The list of node names in the NodesSection.
    */
    QList<quint64> nodeNames;
//    QVector<Eigen::Vector3d> nodeCoordinates;
//    QVector<Node> nodes;
    QMultiHash<quint64,Node> nodes;

//    static  bool lessThan( const Node & node1, const Node & node2 )
//     {
//         if(node1.index<node2.index)
//         {return true;}
//         return false;
//     }

static  bool lessThan( const int & node1, const int & node2 )
    {
        if(node1<node2)
        {return true;}
        return false;
    }
//void sortNodes(){std::sort( nodes.begin(), nodes.end(), lessThan );
//                 sorted = true;
//                }
void sortNodeNames(){std::sort( nodeNames.begin(), nodeNames.end(), lessThan );
                 sorted = true;
                }
   bool getSorted(){return sorted;}

private:

    bool handleNodesListLine(const QString& line);
    bool handleParameterSectionLine(const QString& line);
    bool validLine;
    bool readingNodesList = false;
    QStringList meshFiles;
    QStringList meshFileAlias;
    bool sorted = false;
    QVector<quint64> nonUniqueIdentifiers;
};

#endif // NODES_H
