#ifndef PROJECTFILEHANDLER_H
#define PROJECTFILEHANDLER_H

#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QObject>
#include <QTextStream>
#include <iostream>
#include <fstream>
#include <QVector>
#include <QStringView>
#include <QObject>

/**
* \class ProjectFileHandler
* \brief This class provides methods for reading, writing, and creating a project file.
*/
class ProjectFileHandler : public QObject
{
    Q_OBJECT
public:
    explicit ProjectFileHandler(QObject *parent = 0);

    /**
    * \brief Read a project file with the given file name.
    * \param[in] fileName The file name of the project file to read.
    */
    void readProjectFile(QString fileName);

    /**
    * \brief Write the current project file.
    */
    void writeProjectFile();


    /**
    * \brief Create a new project file with the given file name.
    * \param[in] fileName The file name of the new project file.
    */
    void createProjectFile(QString fileName);

    /**
    * \brief Get the name of the info file for the current project.
    * \return The name of the info file.
    */
    QString getInfoFile();

    /**
     * \brief Get the name of the solving script for the current project.
     * \return The name of the solving script.
     */
    QString getSolvingScript();

    /**
    * \brief Get the name of the direct sound script for the current project.
    * \return The name of the direct sound script.
    */
    QString getDirektSoundScript();

    /**
    * \brief Get the name of the LE script for the current project.
    * \return The name of the LE script.
    */
    QString getLEScript();

    /*!
    * \brief Get the list of mesh files for the current project.
    * \return The list of mesh files.
    */
    QStringList getMeshFileList();

    /*!
    * \brief Get the list of mesh file aliases for the current project.
    * \return
    */
    QStringList getMeshFileAliasList();

    /**
    * \brief Get the list of observation scripts for the current project.
    * \return The list of observation scripts.
    */
    QStringList getObservScriptList();

    /**
    * \brief Set the name of the solving script for the current project.
    * \param[in] Scriptname_Solving The name of the solving script.
    */
    void setSolvingScript(QString Scriptname_Solving);

    /**
    * \brief Set the list of mesh files for the current project.
    * \param[in] meshFiles The list of mesh files.
    */
    void setMeshFileList(QStringList meshFiles);

    /**
    * \brief Set the list of mesh file aliases for the current project.
    * \param[in] meshFileAlias The list of mesh file aliases.
    */
    void setMeshFileAliasList(QStringList meshFileAlias);


    /**
    * \brief Set the list of observation scripts for the current project.
    * \param[in] observationFiles The list of observation scripts.
    */
    void setObservScriptList(QStringList observationFiles);

private:    
    /**
    * \brief Make the given path absolute.
    * \param[in] path The path to make absolute.
    * \return The absolute path.
    */
    QString makePathAbsolute(QString path);

    /**
    * \brief Make the given list of paths absolute.
    * \param[in] path The list of paths to make absolute.
    * \return The list of absolute paths.
    */
    QStringList makePathListAbsolute(QStringList path);

    //! The file name of the current project file.
    QString projectFile;

    //! The name of the info file for the current project.
    QString Scriptname_InfoFile;

    //! The file name of the current project file.
    QString Scriptname_Solving;

    //! The name of the direct sound script for the current
    QString Scriptname_DirectSound;

    //! The name of the LE script for the current project.
    QString Scriptname_LEScript;

    //! The list of mesh files for the current project.
    QStringList meshFiles;

    //! The list of mesh file aliases for the current project.
    QStringList meshFileAlias;

    //! The list of observation scripts for the current project.
    QStringList observationFiles;
};

#endif // PROJECTFILEHANDLER_H
