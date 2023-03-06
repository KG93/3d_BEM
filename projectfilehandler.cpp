#include "projectfilehandler.h"

ProjectFileHandler::ProjectFileHandler(QObject *parent) : QObject(parent)
{

}

void ProjectFileHandler::readProjectFile(QString fileName)
{
    if ( fileName.isEmpty() )
    {
        std::cout <<"Empty filename in readProjectFile"<<std::endl;
        return;
    }
    meshFileAlias.clear();
    meshFiles.clear();
    observationFiles.clear();

    QFile projFile(fileName);
    if (!projFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        return;
        //    emit nachricht("Konnte die Datei nicht einlesen: "+mshFilename);
    }
    projectFile=fileName;
    QTextStream in(&projFile);
    QString line = in.readLine();
    if(line=="; ABEC Project File"){
            std::cout <<"correct .abec file : "<< fileName.toUtf8().constData()<<std::endl;
    }
//    else{            std::cout <<"incorrect .abec file : "<< fileName.toUtf8().constData()<<std::endl;
//    }

    bool meshOrObservationSection=false;
    QString string;
    while (!line.isNull())
    {
        if(!meshOrObservationSection) //hack because no peek() in QTextStream
        {
            line = in.readLine();
        }
        meshOrObservationSection=false;
        if(line=="[Project]")
        {
            string = in.readLine();
            Scriptname_InfoFile = string.mid(20);
            std::cout <<"Scriptname_InfoFile="<< Scriptname_InfoFile.toUtf8().constData()<<std::endl;

        }
        if(line=="[Solving]")
        {
            string = in.readLine();
            Scriptname_Solving = string.mid(19);
            std::cout <<"Scriptname_Solving="<< Scriptname_Solving.toUtf8().constData()<<std::endl;

        }
        if(line=="[DirectSound]")
        {
            string = in.readLine();
            Scriptname_DirectSound = string.mid(23);
            std::cout <<"Scriptname_DirectSound="<< Scriptname_DirectSound.toUtf8().constData()<<std::endl;

        }
        if(line=="[LEScript]")
        {
            string = in.readLine();
            Scriptname_LEScript = string.mid(20);
            std::cout <<"Scriptname_LEScript="<< Scriptname_LEScript.toUtf8().constData()<<std::endl;

        }
        if(line=="[Export]")
        {
//            string = in.readLine();
//             = string.mid(20);
        }
        if(line=="[Observation]")
        {
            string = in.readLine();
            int i=0;
            while(string.left(1)=="C")
            {
                observationFiles.append(string.mid(3));
                std::cout <<"C"<<i<<"="<< observationFiles[i].toUtf8().constData()<<std::endl;
                i=i+1;
                string = in.readLine();

            }
            line=string; //hack because no peek() in QTextStream
            meshOrObservationSection=true;
        }

//        QString identifier="[MeshFiles]";
        if(line=="[MeshFiles]")
        {
            string = in.readLine();
            int i=0;
            while(string.left(1)=="C")
            {
                QString meshFilePath=string.mid(3);
                int separatorPosition=meshFilePath.indexOf(",");

                meshFiles.append(meshFilePath.left(separatorPosition));
                std::cout <<"C"<<i<<"="<< meshFiles[i].toUtf8().constData()<<std::endl;
                meshFileAlias.append(meshFilePath.mid(separatorPosition+1));
                std::cout <<"Alias "<<i<<"="<< meshFileAlias[i].toUtf8().constData()<<std::endl;

                i=i+1;
                string = in.readLine();
            }
            line=string; //hack because no peek() in QTextStream
            meshOrObservationSection=true;

        }
    }
}

QString ProjectFileHandler::getInfoFile()
{
    return makePathAbsolute(Scriptname_InfoFile);
}

QString ProjectFileHandler::getSolvingScript()
{
    return makePathAbsolute(Scriptname_Solving);
}

QString ProjectFileHandler::getDirektSoundScript()
{
    return makePathAbsolute(Scriptname_DirectSound);
}

QString ProjectFileHandler::getLEScript()
{
    return makePathAbsolute(Scriptname_LEScript);

}
QStringList ProjectFileHandler::getMeshFileList()
{
    return makePathListAbsolute(meshFiles);
}
QStringList ProjectFileHandler::getMeshFileAliasList()
{
    return meshFileAlias;
}
QStringList ProjectFileHandler::getObservScriptList()
{
    return makePathListAbsolute(observationFiles);
}

void ProjectFileHandler::setSolvingScript(QString Scriptname_Solving)
{
    QFileInfo projectFileInfo=QFileInfo(projectFile);
    QDir projDir=QDir(projectFileInfo.absoluteDir());
    this-> Scriptname_Solving=projDir.relativeFilePath(Scriptname_Solving);
}

void ProjectFileHandler::setMeshFileList(QStringList meshFiles)
{
    this-> meshFiles.clear();
    QFileInfo projectFileInfo=QFileInfo(projectFile);
    QDir projDir=QDir(projectFileInfo.absoluteDir());
    for(int i=0;i<meshFiles.size();i++)
    {
        this-> meshFiles.append(projDir.relativeFilePath(meshFiles.at(i)));
    }
}

void ProjectFileHandler::setMeshFileAliasList(QStringList meshFileAlias)
{
    this-> meshFileAlias=meshFileAlias;
}

void ProjectFileHandler::setObservScriptList(QStringList observationFiles)
{
    this-> observationFiles.clear();
    QFileInfo projectFileInfo=QFileInfo(projectFile);
    QDir projDir=QDir(projectFileInfo.absoluteDir());
    for(int i=0;i<observationFiles.size();i++)
    {
        this-> observationFiles.append(projDir.relativeFilePath(observationFiles.at(i)));
    }
}

QString ProjectFileHandler::makePathAbsolute(QString path)
{
    QFileInfo projectFileInfo=QFileInfo(projectFile);
    QDir absDir=projectFileInfo.absoluteDir();
    QString absPath=(absDir.absolutePath()).append("/");
    absPath.append(path);
    return absPath;
}

QStringList ProjectFileHandler::makePathListAbsolute(QStringList path)
{
    QStringList absPathList;
    QFileInfo projectFileInfo=QFileInfo(projectFile);
    QDir absDir=projectFileInfo.absoluteDir();
    for(int i=0;i<path.length();i++)
    {
        QString absPath=(absDir.absolutePath()).append("/");
        absPathList.append(absPath.append(path.at(i)));
    }
    return absPathList;
}

void ProjectFileHandler::writeProjectFile(/*QString fileName*/)
{
    if ( projectFile.isEmpty() )
    {
        std::cout <<"Empty projectFileName in writeProjectFile"<<std::endl;
        return;
    }
    QFile file(projectFile);
    if ( file.open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Text) )
    {
        QTextStream stream( &file );
        stream << "; ABEC Project File" << Qt::endl;
        stream << "[Project]"<< Qt::endl;
        stream << "Scriptname_InfoFile="<<Scriptname_InfoFile<< Qt::endl;
        stream <<"[Solving]"<< Qt::endl;
        stream <<"Scriptname_Solving="<<Scriptname_Solving<< Qt::endl;
        stream <<"[DirectSound]"<< Qt::endl;
        stream <<"Scriptname_DirectSound="<<Scriptname_DirectSound<< Qt::endl;
        stream <<"[LEScript]"<< Qt::endl;
        stream <<"Scriptname_LEScript="<<Scriptname_LEScript<< Qt::endl;
        stream <<"[DirectSound]"<< Qt::endl;
        stream <<"Scriptname_DirectSound="<<Scriptname_DirectSound<< Qt::endl;
        stream <<"[Observation]"<< Qt::endl;
        for(int i=0;i<observationFiles.size();i++)
        {
            stream <<"C"<<i<<"="<<observationFiles.at(i)<< Qt::endl;
        }
        stream <<"[MeshFiles]"<< Qt::endl;
        for(int i=0;i<meshFiles.size();i++)
        {
            stream <<"C"<<i<<"="<<meshFiles.at(i)<<","<<meshFileAlias.at(i)<< Qt::endl;
        }
    }
}

void ProjectFileHandler::createProjectFile(QString fileName)
{
    if ( fileName.isEmpty() )
    {
        std::cout <<"Empty filename in createProjectFile"<<std::endl;
        return;
    }
    projectFile=QString(fileName);
    QFile file( fileName );
    if ( file.open(QIODevice::ReadWrite) )
    {
        QTextStream stream( &file );
        stream << "; ABEC Project File" << Qt::endl;
        stream << "[Project]"<< Qt::endl;
        stream << "Scriptname_InfoFile="<< Qt::endl;
        stream <<"[Solving]"<< Qt::endl;
        stream <<"Scriptname_Solving="<< Qt::endl;
        stream <<"[DirectSound]"<< Qt::endl;
        stream <<"Scriptname_DirectSound="<< Qt::endl;
        stream <<"[LEScript]"<< Qt::endl;
        stream <<"Scriptname_LEScript="<< Qt::endl;
        stream <<"[DirectSound]"<< Qt::endl;
        stream <<"Scriptname_DirectSound="<< Qt::endl;
        stream <<"[Observation]"<< Qt::endl;
        stream <<"[MeshFiles]"<< Qt::endl;
    }
}
