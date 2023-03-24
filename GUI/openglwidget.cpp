#include "openglwidget.h"
#include <GL/glu.h>
openglWidget::openglWidget(QWidget *parent) : QOpenGLWidget(parent)
{
}

void openglWidget::initializeGL()
{
    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glClearColor(0,0,0,1);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDepthFunc(GL_LEQUAL);
    glDepthRange(0.0f, 1.0f);
//    glEnable(GL_LIGHT0);
//    glEnable(GL_LIGHTING);
//    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA,  GL_ONE_MINUS_SRC_COLOR );
    initialized = true;
}

void openglWidget::setBoundaryElements(const BoundaryElements &newBoundaryElements)
{
    boundaryElements = newBoundaryElements;
    drawFieldSolution = false;
    if(boundaryElements.soundPressure.size() != 0)
    {
        maxDbOnBoundary = std::log10(std::abs((boundaryElements.soundPressure(0)/referenceSoundpressure)))*20;
        minDbOnBoundary = std::log10(std::abs((boundaryElements.soundPressure(0)/referenceSoundpressure)))*20;
        for(int i=0;i<boundaryElements.soundPressure.size();i++)
        {
            if(abs(boundaryElements.soundPressure(i)) < referenceSoundpressure)
            {
                boundaryElements.soundPressure(i) = referenceSoundpressure;
            }
//            double tmp=std::abs((boundaryElements.phiSolution(i)));
            double tmp = std::log10(std::abs((boundaryElements.soundPressure(i)/referenceSoundpressure)))*20;
            if(tmp>maxDbOnBoundary)
            {
                maxDbOnBoundary = tmp;
            }
            if(tmp<minDbOnBoundary)
            {
                minDbOnBoundary = tmp;
            }
        }
//        std::cout<<"Max soundpressure on element: "<<maxDbOnBoundary<<std::endl;
//        std::cout<<"Min soundpressure on element: "<<minDbOnBoundary<<std::endl;
    }
    solutionColoring();
}

void openglWidget::setObservationFields(const QVector<ObservationField> &obsFields)
{
    int index = -1; // the value stays -1 if no field with solution exists
    this->observationFields = obsFields;
    drawFieldSolution = false;
    for(int i=0; i<observationFields.size(); i++)
    {
        if(observationFields.at(i).soundPressure.size() == observationFields.at(i).triangles.size() && observationFields.at(i).triangles.size() != 0)
        {
            index = i; // --> observationFields.at(index).phiSolution.size()>0
            drawFieldSolution = true;
            continue;
        }
        else if(observationFields.at(i).triangleMidPoints.size() != observationFields.at(i).triangles.size() && observationFields.at(i).triangles.size() != 0)
        {
            drawFieldSolution = false;
            return;
        }
    }
    if(drawFieldSolution)
    {
        maxDbOnField = std::log10(std::abs((observationFields.at(index).soundPressure(0)/referenceSoundpressure)))*20;
        minDbOnField = std::log10(std::abs((observationFields.at(index).soundPressure(0)/referenceSoundpressure)))*20;

        for(int j=0; j<observationFields.size(); j++)
        {
            if(observationFields.at(j).soundPressure.size() != observationFields.at(j).triangles.size())
            {
                break;
            }
            observationFields[j].triangleColors.resize(observationFields.at(j).triangleMidPoints.size());
            for(int i=0; i<observationFields.at(j).phiSolution.rows(); i++)
            {
                if(abs(observationFields.at(j).soundPressure(i)) < referenceSoundpressure)
                {
                    observationFields[j].soundPressure(i)  =referenceSoundpressure;
                }
                double tmp = std::log10(std::abs((observationFields.at(j).soundPressure(i)/referenceSoundpressure)))*20;
                if(tmp > maxDbOnField)
                {
                    maxDbOnField = tmp;
                }
                if(tmp < minDbOnField)
                {
                    minDbOnField = tmp;
                }
            }
        }
    }
    solutionColoring();
    update();
}

void openglWidget::setElementSections(const QList<ElementSection>& elementsSections)
{
    this->elementsSections = elementsSections;
    drawFieldSolution = false;
    selectedTriangle = VectorTriangle(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
}

void  openglWidget::setCenterOfMass(Eigen::Vector3d centerOfMass)
{
    centerOfMassX = centerOfMass(0);
    centerOfMassY = centerOfMass(1);
    centerOfMassZ = centerOfMass(2);
}

void openglWidget::setCameraRadius(double radius)
{
    if(radius == 0)
    {
        radius = 0.1;
    }

    cameraDistance = 2.0 * radius;
    scale = radius;
    int width = widthWidget;
    int height = heightWidget;

    resize(width-1,height-1);//neccessary to adapt scale
    resize(width,height);
}

void openglWidget::paintGL()
{
    if(!initialized)
    {
        initializeGL();
    }
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, widthWidget, heightWidget);

    setUpCamera();

    drawPointSources();
    if(drawField){drawFields();}
    drawObservationPoints();
    if(!showBoundariesColouredBySolution)
    {
        drawElementsSections();
    }
    else
    {
        drawBoundaryElements();
    }
    if(showLegend && (showBoundariesColouredBySolution || drawFieldSolution))
    {
        drawSelectedTriangleContour();
    }
    if(boundaryElements.impedancePlanes.length() >= 1)
    {
        drawImpedancePlanes();
    }
    if(global::activeProgramming) // for testing purposes
    {
//        QVector<Cluster*> clusters = clusterTree.getClusters();
//        if(clusters.length()>=1)
//        {

//        Cluster cluster = *clusters.at(12);
//                   glColor3f(1,0,0);
//                drawCuboid(clusters.at(12)->minCuboid);
//            for(int i = 0; i<clusters.length();i++)
//            {
////                float factor = (float)clustertree.getClusters().at(i).clusterTreeDepth /(float)clustertree.maxTreeDepth;
////                glColor3f(factor, 1-factor, 0);
////                if(clusters.at(i)->isLeaf)
//                {
//                    glColor3f(1,0,0);
//                    drawCuboid(clusters.at(i)->minCuboid);
//                }
////                else
//                {
////                    glColor3f(0,1,0);
////                    drawCuboid(clusters.at(i)->minCuboid);
//                }
//            }
//        }
//        else{
//            std::cout << "No cuboids in clusterTree!" << std::endl;
//        }
    }
    storeTransformationMatrices();
    storeZValueAtDoubleClick();
    if(showLegend)
    {
        drawLegend();
    }
}

void openglWidget::resizeGL(int w, int h)
{
    widthWidget = w;
    heightWidget = h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    std::cout<<"scale: "<<scale<<std::endl;
    gluPerspective(45, (float)w / h, 0.1 * scale, 100.0 * scale);
    glMatrixMode(GL_MODELVIEW);
    update();
}

void openglWidget::getMaxCut(int value)
{
    highCut = (100.0 - value) / 100.0;
    solutionColoring();
    update();
}

void openglWidget::getMinCut(int value)
{
    lowCut = (value) / 100.0;
    solutionColoring();
    update();
}

void openglWidget::getMinAndMaxCut(int min, int max)
{
    lowCut = (min) / 100.0;
    highCut = (100.0 - max) / 100.0;
    solutionColoring();
    update();
}

void openglWidget::getNewSolution(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure)
{
    boundaryElements.phiSolution = phiSolution;
    boundaryElements.dPhiSolution = dPhiSolution;
    boundaryElements.soundPressure = soundPressure;
    setBoundaryElements(boundaryElements);
    updateSelectedTriangleSoundPressure();
    update();
}

void openglWidget::getNewSolutionField(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure, QVector<Eigen::VectorXcd> phiSolutionField, QVector<Eigen::VectorXcd> soundPressureField)
{
    boundaryElements.phiSolution = phiSolution;
    boundaryElements.dPhiSolution = dPhiSolution;
    boundaryElements.soundPressure = soundPressure;
    setBoundaryElements(boundaryElements);
    for(int i=0; i<observationFields.size(); i++)
    {
        observationFields[i].phiSolution = phiSolutionField.at(i);
        observationFields[i].soundPressure = soundPressureField.at(i);
    }
    setObservationFields(observationFields);
    updateSelectedTriangleSoundPressure();
    update();
}

void openglWidget::mouseMoveEvent(QMouseEvent *event)
 {
    if (!event->buttons())
    {
        return;
    }
    if (!mousePressed)
    {
//            std::cout <<"Mouse released."<<std::endl;
        QPointF diff = event->pos() - mouseClickPosition;
        cameraAngle1 = cameraAngle1 + 0.005 * diff.x();
        cameraAngle2 = cameraAngle2 + 0.005 * diff.y();
//            std::cout <<"camera angle 1."<<cameraAngle1<<std::endl;
//            std::cout <<"camera angle 2"<<cameraAngle2<<std::endl;
        mouseClickPosition = event->pos();
        update();
        return;
    }
    else
    {
        const int thresholdValue = 10;
        QPoint diff = event->pos() - mouseClickPosition;
        if (mousePressed)
        {
            mousePressed &= diff.x() < thresholdValue; //mausGedrueckt = (mausGedrueckt & differenz.x() < schwellwert)
            mousePressed &= diff.y() < thresholdValue;
            mousePressed &= diff.x() > - thresholdValue;
            mousePressed &= diff.y() > - thresholdValue;
        }
    }
}

void openglWidget::mousePressEvent(QMouseEvent *event)
 {
    if(event->button() == Qt::RightButton)
    {
//        emit customContextMenuRequested(event->pos());
    }

    if (event->buttons() != Qt::LeftButton)
    {
        return;
    }

    mousePressed = true;
    mouseClickPosition = mouseDragPosition = event->pos();
}

void openglWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Left)
        cameraOffsetX = cameraOffsetX + 0.2;
    if (event->key() == Qt::Key_Right)
        cameraOffsetX = cameraOffsetX - 0.2;
    if (event->key() == Qt::Key_Up)
        cameraOffsetZ = cameraOffsetZ + 0.2;
    if (event->key() == Qt::Key_Down)
        cameraOffsetZ = cameraOffsetZ - 0.2;
    update();
}

void openglWidget::wheelEvent(QWheelEvent *event) // Zoomen per Mausrad
{
    if (event->angleDelta().y() > 15)
    {
//    cameraOffsetZ=cameraOffsetZ+0.2;
        cameraDistance = cameraDistance + 0.1 * scale * scalingMultiplier;
    }

    if (event->angleDelta().y() < -15)
    {
//    cameraOffsetZ=cameraOffsetZ-0.2;
        cameraDistance = cameraDistance - 0.1 * scale * scalingMultiplier;
    }
    update();
}

void  openglWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    if (event->buttons() == Qt::LeftButton)
    {
        QPointF mousePos = event->position();
        mouseXInOpenGL=mousePos.x();
        mouseYInOpenGL=heightWidget-mousePos.y();
        update();
        QVector3D windowPos = QVector3D(mousePos.x(), mousePos.y(), mousePosZBufferVal);
        QRect viewPort=QRect(QPoint(0,heightWidget), QPoint(widthWidget,0));
        QVector3D worldPos = windowPos.unproject(modelViewMatrix, projectionMatrix, viewPort);
        Eigen::Vector3d worlPosEigen={worldPos.x(), worldPos.y(), worldPos.z()};
//        std::cout<<"world coordinates of double click: "<<worldPos.x()<<" "<<worldPos.y()<<" "<<worldPos.z()<<std::endl;

        double minFieldDistance = std::numeric_limits<double>::max();
        fieldIndex = 0;
        fieldTriangleIndex = 0;
        fieldTriangleFound = false;
        if(drawField&&drawFieldSolution)
        {
            for(int i=0;i<observationFields.size();i++)
            {
        //        drawFieldSolution
                for(int j=0; j<observationFields.at(i).triangles.size(); j++)
                {
//                    VectorTriangle triangle=observationFields.at(i).triangles.at(j);
                    Eigen::Vector3d triangleMidPoint = observationFields.at(i).triangleMidPoints.at(j);
                    double distance=(worlPosEigen - triangleMidPoint).norm();
                    if(distance < minFieldDistance)
                    {
                        minFieldDistance = distance;
                        fieldIndex = i;
                        fieldTriangleIndex = j;
                        fieldTriangleFound = true;
                    }
                }
            }
        }
        double minBoundaryElementDistance = std::numeric_limits<double>::max();
        BoundaryElementTriangleIndex = 0;
        BoundaryElementTriangleFound = false;
        if(showBoundariesColouredBySolution)
        {
            for(int i = 0; i < boundaryElements.collocationPoints.size(); i++)
            {
                Eigen::Vector3d triangleMidPoint = boundaryElements.collocationPoints.at(i);
                double distance = (worlPosEigen-triangleMidPoint).norm();
                if(distance<minBoundaryElementDistance)
                {
                    minBoundaryElementDistance=distance;
                    BoundaryElementTriangleIndex=i;
                    BoundaryElementTriangleFound=true;
                }
            }
        }
        if(BoundaryElementTriangleFound || fieldTriangleFound)
        {
            if(minBoundaryElementDistance <= minFieldDistance) // triangle on boundary element selected
            {
                fieldTriangleFound = false;
                selectedTriangle=boundaryElements.triangles.at(BoundaryElementTriangleIndex);
                if(BoundaryElementTriangleIndex <= boundaryElements.soundPressure.size())
                {
                    soundPressureOnSelectedTriangle = boundaryElements.soundPressure(BoundaryElementTriangleIndex);
                }
            }
            else // triangle on field selected
            {
                BoundaryElementTriangleFound = false;
                selectedTriangle = observationFields.at(fieldIndex).triangles.at(fieldTriangleIndex);
                if(fieldTriangleIndex <= observationFields.at(fieldIndex).soundPressure.size())
                {
                    soundPressureOnSelectedTriangle = observationFields.at(fieldIndex).soundPressure(fieldTriangleIndex);
                }
            }
        }
        else
        {
            soundPressureOnSelectedTriangle = referenceSoundpressure;
        }
    }
    else
    {
        return;
    }
}

void openglWidget::setUpCamera()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if((cameraAngle2 < 0) | (cameraAngle2 > global::PI))
    {
        float eyeX;
        float eyeZ;
        float eyeY;
        if(cameraAngle2<0)
        {
            cameraAngle2=-cameraAngle2;
             eyeX = centerOfMassX + cameraDistance*std::cos(cameraAngle1+global::PI)*std::sin(cameraAngle2);
             eyeZ = centerOfMassZ + cameraDistance*std::sin(cameraAngle1+global::PI)*std::sin(cameraAngle2);
             eyeY = centerOfMassY + cameraDistance*std::cos(cameraAngle2);
        }
        else if (cameraAngle2> global::PI)
        {
            cameraAngle2=global::PI-(cameraAngle2-global::PI);
             eyeX = centerOfMassX + cameraDistance*std::cos(cameraAngle1+global::PI)*std::sin(cameraAngle2);
             eyeZ = centerOfMassZ + cameraDistance*std::sin(cameraAngle1+global::PI)*std::sin(cameraAngle2);
             eyeY = centerOfMassY + cameraDistance*std::cos(cameraAngle2);
        }
        gluLookAt(eyeX, eyeY, eyeZ, centerOfMassX, centerOfMassY, centerOfMassZ, 0, -1, 0);
    }
    else
    {
        float eyeX = centerOfMassX + cameraDistance*std::cos(cameraAngle1)*std::sin(cameraAngle2);
        float eyeZ = centerOfMassZ + cameraDistance*std::sin(cameraAngle1)*std::sin(cameraAngle2);
        float eyeY = centerOfMassY + cameraDistance*std::cos(cameraAngle2);
        gluLookAt(eyeX, eyeY, eyeZ, centerOfMassX, centerOfMassY, centerOfMassZ, 0, 1, 0);
    }
}

void openglWidget::solutionColoring()
{
    if(showPhase)
    {
        ColorGradient phaseGradient;
        phaseGradient.phaseGradient();
        float r=0, g=0, b=0;
        double animationPhaseShift = /*animate **/ framecounter / (durationOfAnimationPeriod * phaseAnimationFramerate) * 2.0 * global::PI;
//        std::cout<<"animationPhaseShift " << animationPhaseShift <<std::endl;

        if(showBoundariesColouredBySolution)
        {
            boundaryElements.triangleColors.resize(boundaryElements.triangles.size());
            for(int i=0; i<boundaryElements.soundPressure.size(); i++)
            {
                double val = std::fmod(std::atan2(std::real(boundaryElements.soundPressure(i)),std::imag(boundaryElements.soundPressure(i))) + global::PI + animationPhaseShift, 2 * global::PI);
                phaseGradient.getColorAtValue(val, r,g,b);
                boundaryElements.triangleColors[i] = Eigen::Vector4d(r,g,b,1.0);
            }
        }
        if(drawFieldSolution && drawField)
        {
            for(int j=0; j<observationFields.size(); j++)
            {
                observationFields[j].triangleColors.resize(observationFields.at(j).triangles.size());
                if(observationFields.at(j).triangles.size() != observationFields.at(j).soundPressure.size())
                {
                    break;
                }
                for(int i=0;i<observationFields.at(j).triangleColors.size();i++)
                {
                    double val = std::fmod(std::atan2(std::real(observationFields.at(j).soundPressure(i)),std::imag(observationFields.at(j).soundPressure(i))) + global::PI + animationPhaseShift, 2 * global::PI);
                    phaseGradient.getColorAtValue(val, r,g,b);
                    observationFields[j].triangleColors[i] = Eigen::Vector4d(r,g,b,1.0);
                }
            }
        }
    }
    else
    {
        if(showBoundariesColouredBySolution && drawFieldSolution && drawField)
        {
            if(useGlobalMinMaxDB)
            {
                maxSol = std::max(maxGlobalDbOnField, maxGlobalDbOnBoundary);
                minSol = std::min(minGlobalDbOnField, minGlobalDbOnBoundary);
            }
            else
            {
                maxSol = std::max(maxDbOnField,maxDbOnBoundary);
                minSol = std::min(minDbOnField,minDbOnBoundary);
            }
        }
        else if(drawFieldSolution && drawField && (!showBoundariesColouredBySolution))
        {
            if(useGlobalMinMaxDB)
            {
                maxSol = maxGlobalDbOnField;
                minSol = minGlobalDbOnField;
            }
            else
            {
                maxSol = maxDbOnField;
                minSol = minDbOnField;
            }
        }
        else
        {
            if(useGlobalMinMaxDB)
            {
                maxSol = maxGlobalDbOnBoundary;
                minSol = minGlobalDbOnBoundary;
            }
            else
            {
                maxSol=maxDbOnBoundary;
                minSol=minDbOnBoundary;
            }
        }
        double tmpMax=maxSol;
        maxSol=maxSol-highCut*(maxSol-minSol);
        minSol=lowCut*(tmpMax-minSol)+minSol;
        if(maxSol<minSol)
        {
            double tmp=maxSol;
            maxSol=minSol;
            minSol=tmp;
        }
        ColorGradient heatMapGradient;
        heatMapGradient.defaultGradient();
        float r=0, g=0, b=0;
        if(showBoundariesColouredBySolution)
        {
            boundaryElements.triangleColors.resize(boundaryElements.triangles.size());
            for(int i=0;i<boundaryElements.soundPressure.size();i++)
            {
                double val=(std::log10(std::abs((boundaryElements.soundPressure(i)/referenceSoundpressure)))*20-minSol)/(maxSol-minSol);
                heatMapGradient.getColorAtValue(val, r,g,b);
                boundaryElements.triangleColors[i]=Eigen::Vector4d(r,g,b,1.0);
                if(val < 0 || val > 1)
                {
                    boundaryElements.triangleColors[i] = 0.8*boundaryElements.triangleColors.at(i);
                }
            }
        }
        if(drawFieldSolution&&drawField)
        {
            for(int j=0;j<observationFields.size();j++)
            {
                observationFields[j].triangleColors.resize(observationFields.at(j).triangles.size());
                if(observationFields.at(j).triangles.size() != observationFields.at(j).soundPressure.size())
                {
                    break;
                }
                for(int i=0; i<observationFields.at(j).triangleColors.size(); i++)
                {
                    double val = (std::log10(std::abs((observationFields.at(j).soundPressure(i)/referenceSoundpressure)))*20-minSol)/(maxSol-minSol);
                    heatMapGradient.getColorAtValue(val, r,g,b);
                    observationFields[j].triangleColors[i] = Eigen::Vector4d(r,g,b,1.0);
                    if(val < 0 || val > 1)
                    {
                        observationFields[j].triangleColors[i] = 0.8*observationFields.at(j).triangleColors.at(i);
                    }
                }
            }
        }
    }
}

void openglWidget::drawElementsSections()
{
    for(const ElementSection& elements : elementsSections)
    {
        for(const VectorTriangle& triangle : elements.elementsVectorTriangles) // draw the element triangles
        {
        glBegin(GL_TRIANGLES);
            glColor3f(triangle.color(0),triangle.color(1),triangle.color(2));
            glVertex3f(triangle.node1(0), triangle.node1(1), triangle.node1(2));
            glVertex3f(triangle.node2(0), triangle.node2(1), triangle.node2(2));
            glVertex3f(triangle.node3(0), triangle.node3(1), triangle.node3(2));
        glEnd();
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        glBegin(GL_TRIANGLES);
        glColor3f(0.2, 0.2, 0.2);
        glVertex3f(triangle.node1(0), triangle.node1(1), triangle.node1(2));
        glVertex3f(triangle.node2(0), triangle.node2(1), triangle.node2(2));
        glVertex3f(triangle.node3(0), triangle.node3(1), triangle.node3(2));
        glEnd();
        glPolygonMode(GL_FRONT, GL_FILL);
        }
        for(const VectorQuadrilateral& quadrilateral : elements.elementsVectorQuadrilaterals) // draw the element quadrilaterals
        {
        glBegin(GL_QUADS);
            glColor3f(quadrilateral.color(0),quadrilateral.color(1),quadrilateral.color(2));
            glVertex3f(quadrilateral.node1(0), quadrilateral.node1(1), quadrilateral.node1(2));
            glVertex3f(quadrilateral.node2(0), quadrilateral.node2(1), quadrilateral.node2(2));
            glVertex3f(quadrilateral.node3(0), quadrilateral.node3(1), quadrilateral.node3(2));
            glVertex3f(quadrilateral.node4(0), quadrilateral.node4(1), quadrilateral.node4(2));
        glEnd();
        }
    }
}

void openglWidget::drawBoundaryElements()
{
    for(int i = 0; i < boundaryElements.triangles.size(); i++)
    {
        VectorTriangle triangle = boundaryElements.triangles.at(i);
        float r = boundaryElements.triangleColors.at(i)(0);
        float g = boundaryElements.triangleColors.at(i)(1);
        float b = boundaryElements.triangleColors.at(i)(2);
        glBegin(GL_TRIANGLES);
        glColor3f(r, g, b);
        glVertex3f(triangle.node1(0), triangle.node1(1), triangle.node1(2));
        glColor3f(r, g, b);
        glVertex3f(triangle.node2(0), triangle.node2(1), triangle.node2(2));
        glColor3f(r, g, b);
        glVertex3f(triangle.node3(0), triangle.node3(1), triangle.node3(2));
        glEnd();
    }
}

void openglWidget::drawFields()
{
    glDisable(GL_CULL_FACE);
    for(int i = 0; i < observationFields.size(); i++)
    {
        for(int j = 0; j < observationFields.at(i).triangles.size(); j++)
        {
            VectorTriangle triangle = observationFields.at(i).triangles.at(j);
            if(drawFieldSolution) // vizualize the solution on the field
            {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glBegin(GL_TRIANGLES);
                observationFields[i].triangleColors[j];
                float r=observationFields.at(i).triangleColors.at(j)(0);
                float g=observationFields.at(i).triangleColors.at(j)(1);
                float b=observationFields.at(i).triangleColors.at(j)(2);
                float alpha=observationFields.at(i).triangleColors.at(j)(3);
                glColor4f(r, g, b,alpha);
                glVertex3f(triangle.node1(0), triangle.node1(1), triangle.node1(2));
                glColor4f(r, g, b,alpha);
                glVertex3f(triangle.node2(0), triangle.node2(1), triangle.node2(2));
                glColor4f(r, g, b,alpha);
                glVertex3f(triangle.node3(0), triangle.node3(1), triangle.node3(2));
                glEnd();
            }
            else // draw the field contours
            {
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
                glBegin(GL_TRIANGLES);
                glColor3f(0.2, 0.2, 0.6);
                glVertex3f(triangle.node1(0), triangle.node1(1), triangle.node1(2));
                glVertex3f(triangle.node2(0), triangle.node2(1), triangle.node2(2));
                glVertex3f(triangle.node3(0), triangle.node3(1), triangle.node3(2));
                glEnd();
                glPolygonMode(GL_FRONT, GL_FILL);
            }
        }
    }
    glEnable(GL_CULL_FACE);
}

void openglWidget::drawObservationPoints()
{
    glPointSize(5);
    for(int i = 0; i < observationPoints.size(); i++)
    {
        Eigen::Vector3d coordinates = observationPoints.at(i).coordinates;
        glBegin(GL_POINTS);
        glPointSize(100);
        glColor3f(0, 0, 1);
        glVertex3f(coordinates(0), coordinates(1), coordinates(2));
        glEnd();
    }
}

void openglWidget::drawPointSources()
{
    glPointSize(5);
    for(const PointSource& pointsource : pointSources) // draw all point sources
    {
        glBegin(GL_POINTS);
        glPointSize(100);
        glColor3f(1, 0, 0);
        glVertex3f(pointsource.position(0), pointsource.position(1), pointsource.position(2));
        glEnd();
    }
}

void openglWidget::drawImpedancePlanes()
{
    Eigen::Vector3d centerOfGeometry = {centerOfMassX, centerOfMassY, centerOfMassZ};
    glDisable(GL_CULL_FACE);
    glBegin(GL_QUADS);
    for(int i=0; i<boundaryElements.impedancePlanes.length(); i++)
    {
        ImpedancePlane tmpIpedancePlane = boundaryElements.impedancePlanes.at(i);
        tmpIpedancePlane.halfSpacePlaneNormal = tmpIpedancePlane.halfSpacePlaneNormal.normalized();
        Eigen::Vector3d centerProjectedOnPlane = centerOfGeometry - (centerOfGeometry - tmpIpedancePlane.halfSpacePlanePoint).dot(tmpIpedancePlane.halfSpacePlaneNormal) * tmpIpedancePlane.halfSpacePlaneNormal;

        Eigen::Vector3d tmp = tmpIpedancePlane.halfSpacePlaneNormal;
        long minCoeff = 0;
        tmp.array().abs().minCoeff(&minCoeff);
        tmp(minCoeff) += tmp.norm();

        Eigen::Vector3d planeVector1 = tmp.cross(tmpIpedancePlane.halfSpacePlaneNormal).normalized();
        Eigen::Vector3d planeVector2 = planeVector1.cross(tmpIpedancePlane.halfSpacePlaneNormal).normalized();

        planeVector1 *= 3*scale;
        planeVector2 *= 3*scale;

        Eigen::Vector3d corner1 = centerProjectedOnPlane + planeVector1 + planeVector2;
        Eigen::Vector3d corner2 = centerProjectedOnPlane + planeVector1 - planeVector2;
        Eigen::Vector3d corner3 = centerProjectedOnPlane - planeVector1 - planeVector2;
        Eigen::Vector3d corner4 = centerProjectedOnPlane - planeVector1 + planeVector2;


        glColor4f(0.1, 0.1, 0.3, 0.4);
        glVertex3f(corner1(0), corner1(1), corner1(2));
        glVertex3f(corner2(0), corner2(1), corner2(2));
        glVertex3f(corner3(0), corner3(1), corner3(2));
        glVertex3f(corner4(0), corner4(1), corner4(2));
    }
    glEnd();
    glEnable(GL_CULL_FACE);
}

void openglWidget::drawLegend()
{
    ColorGradient heatMapGradient;    // Used to create an array of different colors.
    heatMapGradient.defaultGradient();
    float r=0, g=0, b=0;
    QVector<double> dbValues;
    QVector<QRectF> legendFieldBoundaries;
    for(int i = 0; i < heatMapGradient.getColor().size(); i++)
    {
        double value = heatMapGradient.getColor().at(i).val;
        double db = value * (maxSol - minSol) + minSol;
        dbValues.push_back(db);
    }
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, widthWidget, 0.0, heightWidget);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glLineWidth(3);
    glBegin(GL_QUADS);
    int fieldWidth = 0.05 * widthWidget;
    int fieldHeight = 0.05 * heightWidget;
    for(int i = 0; i < dbValues.size(); i++)
    {
        heatMapGradient.getColorAtValue(heatMapGradient.getColor().at(i).val,r,g,b);
        glColor3f(r, g, b);
        glVertex3f(widthWidget-i*fieldWidth, heightWidget, 0);
        glVertex3f(widthWidget-(i+1)*fieldWidth, heightWidget, 0);
        glVertex3f(widthWidget-(i+1)*fieldWidth, heightWidget-fieldHeight, 0);
        glVertex3f(widthWidget-i*fieldWidth, heightWidget-fieldHeight, 0);
        legendFieldBoundaries.push_back(QRectF(QPointF(widthWidget-(i+1)*fieldWidth, 0),QPointF(widthWidget-i*fieldWidth, fieldHeight)));
    }
    legendFieldBoundaries.push_back(QRectF(QPointF(widthWidget-(dbValues.size()+1)*fieldWidth, 0),QPointF(widthWidget-dbValues.size()*fieldWidth, fieldHeight)));
    glEnd();
    glLineWidth(1);
    pushRelevantOpenGLAttributes();
    QPainter painter;
    painter.begin(this);

    QFont font = painter.font();
    font.setPixelSize(10);
    painter.setFont(font);
    painter.setRenderHint(QPainter::Antialiasing);
    for(int i=0;i<dbValues.size();i++)
    {
        painter.drawText(legendFieldBoundaries.at(i), /*Qt::AlignLeft|*/Qt::AlignCenter,QString::number(dbValues.at(i))+"db");
    }
    QPen pen=painter.pen();
    pen.setColor(Qt::white);
    painter.setPen(pen);
    double soundPressureOnSelectedTriangleInDb=std::log10(std::abs((soundPressureOnSelectedTriangle/referenceSoundpressure)))*20;
    painter.drawText(legendFieldBoundaries.at(legendFieldBoundaries.size()-1), /*Qt::AlignLeft|*/Qt::AlignCenter,QString::number(soundPressureOnSelectedTriangleInDb)+"db");

    painter.end();
    popRelevantOpenGLAttributes();
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

void openglWidget::storeTransformationMatrices()
{
    float modelview[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, modelview);
    float projection[16];
    glGetFloatv(GL_PROJECTION_MATRIX, projection);

    modelViewMatrix = QMatrix4x4(modelview).transposed(); //OpenGl matrices are saved in column-major order, therefore ...
    projectionMatrix = QMatrix4x4(projection).transposed();
}

void openglWidget::storeZValueAtDoubleClick()
{
    GLfloat depth;
    glReadPixels(mouseXInOpenGL, mouseYInOpenGL, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
    mousePosZBufferVal = depth;
}

void openglWidget::drawSelectedTriangleContour()
{
    glDisable(GL_CULL_FACE);
//    glDepthFunc(GL_ALWAYS);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glLineWidth(3);
    glBegin(GL_TRIANGLES);
    glColor3f(0.0, 0.0, 0.0);
    glVertex3f(selectedTriangle.node1(0), selectedTriangle.node1(1), selectedTriangle.node1(2));
    glVertex3f(selectedTriangle.node2(0), selectedTriangle.node2(1), selectedTriangle.node2(2));
    glVertex3f(selectedTriangle.node3(0), selectedTriangle.node3(1), selectedTriangle.node3(2));
    glEnd();
    glPolygonMode(GL_FRONT, GL_FILL);
//    glDepthFunc(GL_LEQUAL);
    glLineWidth(1);
    glEnable(GL_CULL_FACE);
}

void openglWidget::pushRelevantOpenGLAttributes()
{
    glPushAttrib(GL_ACCUM_BUFFER_BIT);
    glPushAttrib(GL_VIEWPORT_BIT);
    glPushAttrib(GL_TRANSFORM_BIT);
    glPushAttrib(GL_POLYGON_BIT);

    glPushAttrib(GL_PIXEL_MODE_BIT);
    glPushAttrib(GL_MULTISAMPLE_BIT);
    glPushAttrib(GL_LIGHTING_BIT);
    glPushAttrib(GL_ENABLE_BIT);

    glPushAttrib(GL_DEPTH_BUFFER_BIT);
    glPushAttrib(GL_CURRENT_BIT);
    glPushAttrib(GL_COLOR_BUFFER_BIT);
}

void openglWidget::popRelevantOpenGLAttributes()
{
    glPopAttrib();
    glPopAttrib();
    glPopAttrib();
    glPopAttrib();

    glPopAttrib();
    glPopAttrib();
    glPopAttrib();
    glPopAttrib();

    glPopAttrib();
    glPopAttrib();
    glPopAttrib();
}

//void openglWidget::setupClustertree()
//{
//    boundaryElements.calculateTriangleMidPoints();
////    clusterTree.boundaryElements = &boundaryElements;
//    clusterTree.assembleClustertree(&boundaryElements);

//    std::cout<<"Clustertree populated!"<<std::endl;
//    blockClusterTree.setClusterTrees(&clusterTree, &clusterTree);
//    blockClusterTree.populateBlockClusterTree(true);
//    std::cout<<"BlockClustertree populated!"<<std::endl;
//    int numberOfTriangles = blockClusterTree.getRootBlock()->rowCluster->indices.last() - blockClusterTree.getRootBlock()->rowCluster->indices.first() + 1;
//    std::cout<<"Number of triangles: "<<numberOfTriangles<<std::endl;
//    std::cout<<"Number of blocks: "<<blockClusterTree.getMinPartition().length()<<std::endl;
//}

void openglWidget::setAnimatePhase(bool animate)
{
    this->animate = animate;
    if(animate)
    {
        connect(&animationTimer, &QTimer::timeout, this, &openglWidget::animatePhase);
        animationTimer.start(1000.0/phaseAnimationFramerate);
    }
    else
    {
        animationTimer.stop();
    }
}

void openglWidget::drawCuboid(Cuboid cuboid)
{
    glDisable(GL_CULL_FACE);
    glLineWidth(1.5);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
//    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_QUADS);
    glVertex3f(cuboid.minPoint(0), cuboid.minPoint(1), cuboid.minPoint(2));
    glVertex3f(cuboid.minPoint(0), cuboid.maxPoint(1), cuboid.minPoint(2));
    glVertex3f(cuboid.minPoint(0), cuboid.maxPoint(1), cuboid.maxPoint(2));
    glVertex3f(cuboid.minPoint(0), cuboid.minPoint(1), cuboid.maxPoint(2));

    glVertex3f(cuboid.minPoint(0), cuboid.minPoint(1), cuboid.minPoint(2));
    glVertex3f(cuboid.minPoint(0), cuboid.minPoint(1), cuboid.maxPoint(2));
    glVertex3f(cuboid.maxPoint(0), cuboid.minPoint(1), cuboid.maxPoint(2));
    glVertex3f(cuboid.maxPoint(0), cuboid.minPoint(1), cuboid.minPoint(2));

    glVertex3f(cuboid.minPoint(0), cuboid.maxPoint(1), cuboid.minPoint(2));
    glVertex3f(cuboid.minPoint(0), cuboid.maxPoint(1), cuboid.maxPoint(2));
    glVertex3f(cuboid.maxPoint(0), cuboid.maxPoint(1), cuboid.maxPoint(2));
    glVertex3f(cuboid.maxPoint(0), cuboid.maxPoint(1), cuboid.minPoint(2));

    glVertex3f(cuboid.maxPoint(0), cuboid.minPoint(1), cuboid.minPoint(2));
    glVertex3f(cuboid.maxPoint(0), cuboid.maxPoint(1), cuboid.minPoint(2));
    glVertex3f(cuboid.maxPoint(0), cuboid.maxPoint(1), cuboid.maxPoint(2));
    glVertex3f(cuboid.maxPoint(0), cuboid.minPoint(1), cuboid.maxPoint(2));

    glEnd();
    glPolygonMode(GL_FRONT, GL_FILL);
    glEnable(GL_CULL_FACE);
    glLineWidth(1);
}

void openglWidget::updateSelectedTriangleSoundPressure()
{
    if(showBoundariesColouredBySolution && BoundaryElementTriangleFound)
    {
        if(BoundaryElementTriangleIndex <= boundaryElements.soundPressure.size())
        {
            soundPressureOnSelectedTriangle = boundaryElements.soundPressure(BoundaryElementTriangleIndex);
        }
    }
    else if(drawFieldSolution && fieldTriangleFound)
    {
        if(fieldIndex <= observationFields.size())
        {
            if(fieldTriangleIndex <= observationFields.at(fieldIndex).soundPressure.size())
            {
                soundPressureOnSelectedTriangle = observationFields.at(fieldIndex).soundPressure(fieldTriangleIndex);
            }
        }
    }
    else
    {
        soundPressureOnSelectedTriangle = referenceSoundpressure;
    }
}

void openglWidget::animatePhase()
{
    if(animate && showPhase)
    {
        framecounter++;
        solutionColoring();
        update();
    }
}
