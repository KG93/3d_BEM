#ifndef OPENGLWIDGET_H

#define OPENGLWIDGET_H
#include "../HMatrix/clustertree.h"
#include "../HMatrix/hmatrix.h"
#include "../HMatrix/harithm.h"
#include "../HMatrix/hmultiply.h"
#include "../HMatrix/gmres.h"

#include "Timer.h" //Timer

#include "boundaryelements.h"
#include "vectortriangle.h"
#include "SolvingScript/elementSection.h"
#include "pointsource.h"
#include "colorgradient.h"
#include "ObservationScript/observationpoint.h"
#include "ObservationScript/observationfield.h"

#include <QWidget>
#include <QPainter>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <iostream>
#include <QMouseEvent>
#include <QTimer>

/**
* \brief The openglWidget class is a widget for rendering 3D (BEM) geometry using OpenGL.
*/
class openglWidget : public QOpenGLWidget
{
    Q_OBJECT
public:
    explicit openglWidget(QWidget *parent = nullptr);

    /**
    * \brief initializeGL initializes the OpenGL rendering context.
    */
    void initializeGL();

    /**
    * \brief resizeGL resizes the OpenGL widget.
    *
    * \param w the width of the widget.
    * \param h the height of the widget.
    */
    void resizeGL(int w, int h);

    /**
    * \brief paintGL redraws the widget.
    */
    void paintGL() override;

    /**
    * \brief setBoundaryElements sets the boundary elements to be drawn.
    *
    * \param newBoundaryElements the boundary elements to be drawn.
    */
    void setBoundaryElements(const BoundaryElements& newBoundaryElements);

    /**
    * \brief getBoundaryElements get the boundary elements to be drawn.
    */
    BoundaryElements getBoundaryElements(){return boundaryElements;}

    /**
    * \brief setElementSections sets the element sections to be drawn.
    *
    * \param elementsSections the element sections to be drawn.
    */
    void setElementSections(const QList<ElementSection>& elementsSections);

    /**
    * \brief setPointSources sets the point sources to be drawn.
    *
    * \param newPointSources the point sources to be drawn.
    */
    void setPointSources(QVector<PointSource> newPointSources){this->pointSources = newPointSources;}

    /**
    * \brief setObservationPoints sets the observation points to be drawn.
    *
    * \param observationPoints the observation points to be drawn.
    */
    void setObservationPoints(QVector<ObservationPoint> observationPoints){this->observationPoints = observationPoints;}

    /**
    * \brief setObservationFields sets the observation fields to be drawn.
    *
    * \param observationFields the observation fields to be drawn.
    */
    void setObservationFields(QVector<ObservationField> observationFields);

    /**
    * \brief setCenterOfMass sets the center of mass of the geometry. The camera rotates around the center of mass.
    *
    * \param centerOfMass the center of mass of the geometry.
    */
    void setCenterOfMass(Eigen::Vector3d centerOfMass);
    void setCameraRadius(double radius);

    /**
    * \brief Returns a boolean indicating whether the solution is currently shown on the boundaries.
    * \return A boolean indicating whether the solution is currently shown on the boundaries.
    */
    bool getShowSolutionValue(){ return showBoundariesColouredBySolution;}

    /**
    * \brief Sets the value indicating whether the solution should be shown on the boundaries.
    * \param value A boolean indicating whether the solution should be shown on the boundaries.
    */
    void setShowSolutionValue(bool value){ this->showBoundariesColouredBySolution = value; solutionColoring();}

    /**
    * \brief Returns a boolean indicating whether the phase should be shown on the boundaries.
    * \return A boolean indicating whether the phase should be shown on the boundaries.
    */
    bool getShowPhaseValue(){ return showPhase;}

    /**
    * \brief Sets the value indicating whether the phase should be shown on the boundaries.
    * \param value A boolean indicating whether the phase should be shown on the boundaries.
    */
    void setShowPhaseValue(bool value){ this->showPhase = value; solutionColoring();}

    /**
    * \brief Gets the value of whether to draw the field or not
    * \return Whether to draw the field or not
    */
    bool getDrawField(){ return drawField;}

    /**
    * \brief Sets the value of whether to draw the field or not
    * \param value The value to set whether to draw the field or not
    */
    void setDrawField(bool value){ this->drawField = value; solutionColoring();}

    /**
    * \brief Gets the value of whether to draw the legend or not
    * \return Whether to draw the legend or not
    */
    bool getDrawLegend(){ return showLegend;}

    /**
    * \brief Sets the value of whether to draw the legend or not
    * \param value The value to set whether to draw the legend or not
    */
    void setDrawLegend(bool value){ this->showLegend = value;}

//    /**
//    * \brief Sets up the clustertree.
//    */
//    void setupClustertree();

    /**
    * \brief Gets the value of whether to animate the phase or not
    * \return Whether to animate the phase or not
    */
    bool getAnimatePhase(){ return animate;}

    /**
    * \brief Sets the value of whether to animate the phase or not
    * \param animate The value to set whether to animate the phase or not
    */
    void setAnimatePhase(bool animate);

    /**
    * \brief Gets the value of whether to use the global (frequency-wise) or local sound pressure extrema for the triangle colouring
    * \return useGlobalMinMaxDB Whether to use the global (frequency-wise) or local sound pressure extrema for the triangle colouring
    */
    bool getUseGlobalMinMaxDb(){return useGlobalMinMaxDB;};

    /**
    * \brief Sets the value of whether to use the global (frequency-wise) or local sound pressure extrema for the triangle colouring
    * \param arg The value to set whether to use the global (frequency-wise) or local sound pressure extrema for the triangle colouring
    */
    void setUseGlobalMinMaxDb(bool arg){this->useGlobalMinMaxDB = arg; solutionColoring(); update();};
    void setGlobalMaxDbOnBoundary(double maxSoundPressure){ this->maxGlobalDbOnBoundary = std::log10(std::abs((maxSoundPressure/referenceSoundpressure)))*20; }
    void setGlobalMinDbOnBoundary(double minSoundPressure){ this->minGlobalDbOnBoundary = std::log10(std::abs((minSoundPressure/referenceSoundpressure)))*20; }
    void setGlobalMaxDbOnField(double maxSoundPressure){ this->maxGlobalDbOnField = std::log10(std::abs((maxSoundPressure/referenceSoundpressure)))*20; }
    void setGlobalMinDbOnField(double minSoundPressure){ this->minGlobalDbOnField = std::log10(std::abs((std::max(minSoundPressure,referenceSoundpressure)/referenceSoundpressure)))*20; }

    int widthWidget;
    int heightWidget;

    double referenceSoundpressure = /*std::sqrt(2.0)**/20*1e-6;
    double staticPressure = 101.325; //Pa at 20°C and sea level
    double airDensity = 1.2041;
    const std::complex<double> imaginaryUnit = std::complex<double>(0.0,1.0);

private:
    QMatrix4x4 modelViewMatrix;
    QMatrix4x4 projectionMatrix;
    int mouseXInOpenGL;
    int mouseYInOpenGL;
    double mousePosZBufferVal;
    VectorTriangle selectedTriangle = VectorTriangle(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());

    std::complex<double> soundPressureOnSelectedTriangle = referenceSoundpressure;
    bool fieldTriangleFound = false;
    int fieldIndex = 0;
    int fieldTriangleIndex = 0;
    bool BoundaryElementTriangleFound = false;
    int BoundaryElementTriangleIndex = 0;

    bool mousePressed;
    QPoint mouseClickPosition;
    QPoint mouseDragPosition;
    float cameraAngle1 = 0;
    float cameraAngle2 = global::PI/2.0;
    float cameraOffsetX   = 0;
    float cameraOffsetY   = 0;
    float cameraOffsetZ   = -2;
    float cameraDistance   = 3.0;
    float scale = 1;
    float scalingMultiplier = 5;
    float cameraFocus   = 0;
    float centerOfMassX = 0;
    float centerOfMassY = 0;
    float centerOfMassZ = 0;
    BoundaryElements boundaryElements;
    QVector<ObservationPoint> observationPoints;
    QVector<ObservationField> observationFields;
//    std::vector<BoundaryElements> boundaryElements;
    //test
    ClusterTree clusterTree;
    HMatrix blockClusterTree;
    //test
    QList<ElementSection> elementsSections;
    QVector<PointSource> pointSources;
    bool initialized = false;
    bool showBoundariesColouredBySolution = false;
    bool showPhase = false;
    bool animate = false;
    bool drawField = false;
    bool drawFieldSolution = false;
    bool showLegend = true;
    double maxSol,minSol;
    double highCut = 0.0;
    double lowCut = 0.0;
    double maxDb = 100;
    double minDb = 0.1;
    double maxDbOnBoundary = 100;
    double minDbOnBoundary = 0.1;
    double maxDbOnField = 100;
    double minDbOnField = 0.1;

    bool useGlobalMinMaxDB = false;
    double maxGlobalDbOnBoundary = 100;
    double minGlobalDbOnBoundary = 0.1;
    double maxGlobalDbOnField = 100;
    double minGlobalDbOnField = 0.1;

    int phaseAnimationFramerate = 60; // animation framerate in Hz
    double durationOfAnimationPeriod = 3; // duration of an animation iteration
    unsigned int framecounter = 0;
    QTimer animationTimer;

    void setUpCamera();
    void solutionColoring();
    void drawElementsSections();
    void drawBoundaryElements();
    void drawFields();
    void drawObservationPoints();
    void drawPointSources();
    void drawImpedancePlanes();
    void drawLegend();
    void storeTransformationMatrices();
    void storeZValueAtDoubleClick();
    void drawSelectedTriangleContour();
    void pushRelevantOpenGLAttributes();
    void popRelevantOpenGLAttributes();
    void drawCuboid(Cuboid cuboid);
    void updateSelectedTriangleSoundPressure();
protected:

signals:

public slots:
    void getMaxCut(int value);
    void getMinCut(int value);
    void getMinAndMaxCut(int min, int max);
    void getNewSolution(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure);
    void getNewSolutionField(Eigen::VectorXcd phiSolution, Eigen::VectorXcd dPhiSolution, Eigen::VectorXcd soundPressure, QVector<Eigen::VectorXcd> phiSolutionField, QVector<Eigen::VectorXcd> soundPressureField);

private slots:
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
    void animatePhase();
};

#endif // OPENGLWIDGET_H
