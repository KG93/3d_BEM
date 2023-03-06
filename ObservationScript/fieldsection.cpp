#include "fieldsection.h"

FieldSection::FieldSection()
{

}

bool FieldSection::handleScriptLine(const QString& currentLine)
{

    validLine=false;
    if(!readingFieldList)
    {
        validLine=handleParameterSectionLine(currentLine);
    }
    else
    {
        validLine=handleFieldListLine(currentLine);
    }
    return validLine;}

bool FieldSection::handleFieldListLine(const QString& currentLine)
{
    validLine=false;
    QString currentLineCopy=currentLine;
    QRegularExpressionMatch match;
    //find RefElements
    QRegularExpression identifier("\\s+RefElements{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);

    int matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
    if(matchPosition>=0) //RefNodes declaration in BE_Spectrum item declaration line found
    {
        QString RefElements=match.captured(1);
        global::removeQuotes(RefElements);
        std::cout <<"RefElements: "<< RefElements.toUtf8().constData()<<std::endl;
        currentLineCopy.remove(matchPosition,match.capturedLength() );
        std::cout <<"Current line with RefNodes removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
        if(RefElements.isEmpty()==false)
        {
            std::cout <<"RefElements: "<< currentLine.toUtf8().constData()<<std::endl;
            QStringList list =global::stringToStringListQuotesIntact(currentLineCopy);
            if(!list.isEmpty())
            {
                bool isNumber;
                int fieldIndex=list.first().toInt(&isNumber);
                if(isNumber)
                {
                    validLine=true;
                    boundaryFields.push_back(BoundaryField(fieldIndex,RefElements));
                }
            }
        }
        else
        {
            validLine=false;
            std::cout <<"No RefElements declared in Line: "<< currentLine.toUtf8().constData()<<" in FieldSection"<< name.toUtf8().constData()<<std::endl;
            logStrings::errorLogString.append("No RefElements declared in Line: "+ currentLine+" in FieldSection"+ name+"\r\n");
        }
    }
    else
    {
        bool isNumber;
        bool isNumber2;
        bool isNumber3;
        bool isNumber4;
        bool isNumber5;
        QStringList list =global::stringToStringListQuotesIntact(currentLineCopy);
        if(!list.isEmpty())
        {
            int fieldIndex=list.first().toInt(&isNumber);
            if(!isNumber)
            {
                validLine=true;
                return validLine;
            }
            else
            {
                if(list.size()>=2)
                {
                    int firstNode=list.at(1).toInt(&isNumber2);
                    if(isNumber2)// declaration of observation field via 4 node numbers
                    {
                        if(list.size()>=5)
                        {
                            int secondNode=list.at(2).toInt(&isNumber3);
                            int thirdNode=list.at(3).toInt(&isNumber4);
                            int fourthNode=list.at(4).toInt(&isNumber5);
                            if(!(isNumber&&isNumber2&&isNumber3&&isNumber4&&isNumber5)){validLine=false; std::cout <<"Invalid observation field declaration "<< currentLine.toUtf8().constData()<<std::endl;
                            return validLine;}
//                            QList<int> nodeIndexes={firstNode,secondNode,thirdNode,fourthNode};
//                            QSet<int> set = nodeIndexes.toSet();
                            if(!global::allIntsUnique(firstNode,secondNode,thirdNode,fourthNode)){validLine=false; std::cout <<"Invalid observation field declaration. Node Indexes should be pairwise distinct: "<< currentLine.toUtf8().constData()<<std::endl;
                                return validLine;}
                            //find RefNodes
                            identifier=QRegularExpression("\\s*RefNodes{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
                            QString localRefNodes=RefNodes;
                            matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
                            if(matchPosition>=0) //DrvGroups declaration in BE_Spectrum item declaration line found
                            {
                                localRefNodes=match.captured(1);
                                global::removeQuotes(localRefNodes);
                                std::cout <<"localRefNodes: "<< localRefNodes.toUtf8().constData()<<std::endl;
                                currentLineCopy.remove(matchPosition,match.capturedLength() );
                                std::cout <<"Current line with RefNodes removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
                            }
                            if(localRefNodes.isEmpty())
                            {
                                std::cout <<"No RefNodes declared for Field with index "<<fieldIndex<<"."<<std::endl;
                                logStrings::errorLogString.append("No RefNodes declared for Field with index "+QString::number(fieldIndex)+".\r\n");
                                validLine=false;
                                return validLine;
                            }
                            validLine=true;
                            nodesFields.push_back(NodesField(fieldIndex,firstNode,secondNode,thirdNode,fourthNode,localRefNodes));
//                            std::cout <<"Field with index "<<fieldIndex<<" declared via nodes list."<<std::endl;

                        }
                        else
                        {
                            std::cout<<"Invalid observation field declaration. Must have four nodes. Line: "<<currentLine.toUtf8().constData()<<std::endl;
                            logStrings::errorLogString.append("Invalid observation field declaration. Must have four nodes. Line: "+currentLine+".\r\n");
                        }


                    }
                    else    // declaration of observation field via meshfile
                    {
                        if(!(currentLine.contains(QRegularExpression("\\s+include([\\s]|$)",QRegularExpression::CaseInsensitiveOption)) || currentLine.contains(QRegularExpression("\\s+Mesh([\\s]|$)",QRegularExpression::CaseInsensitiveOption)) || currentLine.contains(QRegularExpression("\\s+meshfileAlias\\s*=",QRegularExpression::CaseInsensitiveOption))))
                        {
                            validLine=false;
                            std::cout<<"Invalid observation field declaration. Line: "<<currentLine.toUtf8().constData()<<std::endl;
                            logStrings::errorLogString.append("Invalid observation field declaration. Line: "+currentLine+".\r\n");
                            return validLine;
                        }
                        //find MeshFileAlias
                        identifier=QRegularExpression("\\s+MeshFileAlias{1}\\s*={1}\\s*([\"\'][^\"\']+[\"\']|\\S+)", QRegularExpression::CaseInsensitiveOption);
                        QString localMeshFileAlias=MeshFileAlias;
                        matchPosition=currentLineCopy.indexOf(identifier, 0, &match);
                        if(matchPosition>=0) //MeshFileAlias declaration in field declaration line found
                        {
                            localMeshFileAlias=match.captured(1);
                            global::removeQuotes(localMeshFileAlias);
                            std::cout <<"localMeshFileAlias: "<< localMeshFileAlias.toUtf8().constData()<<std::endl;
                            currentLineCopy.remove(matchPosition,match.capturedLength() );
                            std::cout <<"Current line with MeshFileAlias removed: "<< currentLineCopy.toUtf8().constData()<<std::endl;
                        }
                        if(localMeshFileAlias.isEmpty())
                        {
                            std::cout <<"No MeshFileAlias declared for Field with index "<<fieldIndex<<"."<<std::endl;
                            logStrings::errorLogString.append("No MeshFileAlias declared for Field with index "+QString::number(fieldIndex)+".\r\n");
                            validLine=false;
                            return validLine;
                        }
                        list =global::stringToStringListQuotesIntact(currentLineCopy);
                        int includePos=list.indexOf(QRegularExpression("^\\s*include$", QRegularExpression::CaseInsensitiveOption));
                        int excludePos=list.indexOf(QRegularExpression("^\\s*exclude$", QRegularExpression::CaseInsensitiveOption));
                        for (int i=0; i<=list.length();i++)
                        {
                            std::cout << list.value(i).toUtf8().constData()<<std::endl;
                        }
                        bool includeAll=false;
                        QStringList include;
                        QStringList exclude;

                        if( (includePos==-1) | (includePos==list.length()-1) )
                        {
                            includeAll=true;
                        }
                        else
                        {
                            for(int i=includePos+1;i<=list.length()-1;i++)
                            {
                                QString currentListElement=list.value(i);
                                if(currentListElement.contains(QRegularExpression("^\\s*MeshFileAlias",QRegularExpression::CaseInsensitiveOption)))
                                {
                                    break;
                                }
                                if(currentListElement.contains(QRegularExpression("^\\s*exclude\\s*$",QRegularExpression::CaseInsensitiveOption)))
                                {
                                    break;
                                }
                                if(currentListElement.contains(QRegularExpression("^\\s*all\\s*$",QRegularExpression::CaseInsensitiveOption)))
                                {
                                    includeAll=true;
                                    continue;
                                }
                                if(currentListElement.contains(QRegularExpression("^\\s*Mesh\\s*$",QRegularExpression::CaseInsensitiveOption)))
                                {
                                    continue;
                                }
                                include.push_back(currentListElement);
                                std::cout <<"include: "<< currentListElement.toUtf8().constData()<<std::endl;

                            }
                        }
                        if(includeAll){std::cout <<"includeAll=true"<<std::endl;}

                        if(excludePos>=1 && (excludePos<=list.length()-2))
                        {
                            for(int i=excludePos+1;i<=list.length()-1;i++)
                            {
                                QString currentListElement=list.value(i);
                                if(currentListElement.contains(QRegularExpression("^\\s*MeshFileAlias",QRegularExpression::CaseInsensitiveOption)))
                                {
                                    break;
                                }
                                if(currentListElement.contains(QRegularExpression("^\\s*include\\s*$",QRegularExpression::CaseInsensitiveOption)))
                                {
                                    break;
                                }
                                if(currentListElement.contains(QRegularExpression("^\\s*Mesh\\s*$",QRegularExpression::CaseInsensitiveOption)))
                                {
                                    continue;
                                }

                                exclude.push_back(currentListElement);
                                std::cout <<"exclude: "<< currentListElement.toUtf8().constData()<<std::endl;
                            }

                        }
                        meshFields.push_back(MeshField(fieldIndex,localMeshFileAlias,include,exclude,includeAll));
                        validLine=true;
                    }
                }
            }
        }
    }

    return validLine;
}

bool FieldSection::handleParameterSectionLine(const QString& currentLine)
{
    validLine=false;
    QRegularExpressionMatch match;
    QString value;
    QRegularExpression identifier;

    identifier=global::regExWithOneValue("RefNodes");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        RefNodes=value;
        global::removeQuotes(RefNodes);
        if(RefNodes.isEmpty()){validLine=false;}
        else{validLine=true;}
        return validLine;
    }

    identifier=global::regExWithOneValue("MeshFileAlias");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        MeshFileAlias=value;
        global::removeQuotes(MeshFileAlias);
        if(MeshFileAlias.isEmpty()){validLine=false;}
        else{validLine=true;}
        return validLine;
    }

    identifier=global::regExWithOneValue("MeshFrequency");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        MeshFrequency=global::stringWithUnitsToDouble(value,QString("Hz"),validLine);

        if(MeshFrequency<=0){validLine=false;}
        return validLine;
    }

    identifier=global::regExWithOneValue("EdgeLength");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        EdgeLength=global::stringWithUnitsToDouble(value,QString("m"),validLine);
        std::cout<<"Edgelength: "<<EdgeLength<<std::endl;
        if(EdgeLength<=0){EdgeLength=0; validLine=false;}
        return validLine;
    }

    identifier=global::regExWithOneValue("Alpha");
    if(currentLine.contains(identifier,&match))
    {
        value=match.captured(1);
        global::removeQuotesAndWhitespace(value);
        Alpha=value.toDouble(&validLine);
        if(Alpha<0){std::cout<<"Alpha="<<Alpha<< "<0. Alpha will be set to 0."<<std::endl; Alpha=0;}
        if(Alpha>1){std::cout<<"Alpha="<<Alpha<< ">1. Alpha will be set to 1."<<std::endl; Alpha=1;}
        return validLine;
    }


//    identifier=global::regExWithOneValue("scale");
//    if(currentLine.contains(identifier,&match))
//    {
//        value=match.captured(1);
//        global::removeQuotes(value);
//        bool alreadyCaptured=false;
//        QRegularExpression separator( "[ \\t,]+");
//        QStringList list =value.split(separator, Qt::SkipEmptyParts);
//        if(list.length()==1)
//        {
//            QString captured1=list.at(0);
//            double val =global::stringWithUnitsToDouble(captured1,QString("m"),validLine);
//            Scale={val,val,val};
//            alreadyCaptured=true;
//        }
//        if(list.length()==2)
//        {
//            bool valid1,valid2;
//            QString captured1=list.at(0);
//            QString captured2=list.at(1);
//            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
//            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
//            Scale={val1,val2,1};
//            validLine=(valid1 && valid2);
//            alreadyCaptured=true;
//        }
//        if(list.length()==3)
//        {
//            bool valid1,valid2,valid3;
//            QString captured1=list.at(0);
//            QString captured2=list.at(1);
//            QString captured3=list.at(2);
//            double val1 =global::stringWithUnitsToDouble(captured1,QString("m"),valid1);
//            double val2 =global::stringWithUnitsToDouble(captured2,QString("m"),valid2);
//            double val3 =global::stringWithUnitsToDouble(captured3,QString("m"),valid3);
//            Scale={val1,val2,val3};
//            validLine=(valid1 && valid2 && valid3);
//            alreadyCaptured=true;
//        }
//        if(!alreadyCaptured)
//        {
//            validLine=false;
//        }
//        if(Scale(0)==0){Scale(0)=1;}
//        if(Scale(1)==0){Scale(1)=1;}
//        if(Scale(2)==0){Scale(2)=1;}
//        std::cout<<"scale: "<<Scale(0)<<", "<<Scale(1)<<", "<<Scale(2)<<std::endl;
//        return validLine;
//    }
    QRegularExpression separator( "[ \\t,]+"); //identify any number of white spaces and tabs intermixed or commas
    QStringList list =currentLine.split(separator, Qt::SkipEmptyParts);
    if(!list.isEmpty())
    {
        bool isNumber;
        list.first().toInt(&isNumber);
        if(isNumber)
        {
            readingFieldList=true;
            validLine=handleFieldListLine(currentLine);
            return validLine;
        }
    }
    return validLine;
}

bool FieldSection::containsNoElements()
{
    return (nodesFields.isEmpty() && boundaryFields.isEmpty() && meshFields.isEmpty());
}
