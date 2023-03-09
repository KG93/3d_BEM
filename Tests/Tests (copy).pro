QT += testlib core gui openglwidgets svg

TARGET = 3d_BEM_Tests
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

QMAKE_CXXFLAGS += -fno-fast-math #-foffload=nvptx-none #-lomptarget #-fopenacc
QMAKE_CXXFLAGS_RELEASE += -fopenmp -march=native -O3 -fno-fast-math -frounding-math
QMAKE_LFLAGS += -fopenmp -fno-stack-protector #-foffload=nvidia-ptx #-lomptarget #-fopenacc
#QMAKE_LFLAGS +=  -fopenmp -static

SOURCES += tests.cpp    \
    ../GUI/freqlistwidget.cpp \
    ../GUI/logwidget.cpp \
    ../GUI/mainwindow.cpp \
    ../GUI/openglwidget.cpp \
    ../GUI/parameterdialog.cpp \
    ../GUI/qdoublespinboxnotrailzeros.cpp \
    ../GUI/registerscriptstab.cpp \
    ../HMatrix/clustertree.cpp \
    ../HMatrix/cuboid.cpp \
    ../HMatrix/gmres.cpp \
    ../HMatrix/golubReinschSVD.cpp \
    ../HMatrix/harithm.cpp \
    ../HMatrix/hmatrix.cpp \
    ../HMatrix/hmatrixvisuals.cpp \
    ../HMatrix/hmultiply.cpp \
    ../LinearElements/linearboundaryelementssolver.cpp \
    ../LinearElements/lineartriangle.cpp \
    ../LinearElements/lineartrianglenode.cpp \
    ../SolvingScript/infinitebafflesection.cpp \
    ../bemoperatorsubroutines.cpp \
    ../impedanceplane.cpp \
    ../LinearElements/linearboundaryelements.cpp \
    ../linequadraturerules.cpp \
    ../mshreader.cpp \
    ../projectfilehandler.cpp \
    ../triangle.cpp \
    ../SolvingScript/solvingscriptreader.cpp \
    ../boundaryelements.cpp \
    ../boundaryelementsolver.cpp \
    ../vectortriangle.cpp \
    ../quadrilateral.cpp \
    ../nodesSection.cpp \
    ../SolvingScript/elementSection.cpp \
    ../SolvingScript/meshfileelement.cpp \
    ../meshfilepropertiessection.cpp \
    ../SolvingScript/controlsolversection.cpp \
    ../subdomainpropertiessection.cpp \
    ../global.cpp \
    ../trianglequadraturerules.cpp \
    ../SolvingScript/wallimpedancesection.cpp \
    ../gls.cpp \
    ../pointsource.cpp \
    ../SolvingScript/pressurepointssection.cpp \
    ../ObservationScript/observationscriptreader.cpp \
    ../ObservationScript/controlfieldsection.cpp \
    ../ObservationScript/controlspectrumsection.cpp \
    ../ObservationScript/fieldsection.cpp \
    ../ObservationScript/bespectrumsection.cpp \
    ../meshfunctions.cpp \
    ../ObservationScript/observationfield.cpp \
    ../ObservationScript/observationpoint.cpp \
    ../globallogstrings.cpp \
    ../SolvingScript/drivingsection.cpp

HEADERS  += ../GUI/colorgradient.h \
    ../GUI/freqlistwidget.h \
    ../GUI/logwidget.h \
    ../GUI/mainwindow.h \
    ../GUI/openglwidget.h \
    ../GUI/parameterdialog.h \
    ../GUI/qdoublespinboxnotrailzeros.h \
    ../GUI/registerscriptstab.h \
    ../HMatrix/clustertree.h \
    ../HMatrix/cuboid.h \
    ../HMatrix/gmres.h \
    ../HMatrix/golubReinschSVD.h \
    ../HMatrix/harithm.h \
    ../HMatrix/hmatrix.h \
    ../HMatrix/hmatrixvisuals.h \
    ../HMatrix/hmultiply.h \
    ../LinearElements/linearboundaryelementssolver.h \
    ../LinearElements/lineartriangle.h \
    ../LinearElements/lineartrianglenode.h \
    ../SolvingScript/infinitebafflesection.h \
    ../Timer.h \
    ../bemoperatorsubroutines.h \
    ../impedanceplane.h \
    ../LinearElements/linearboundaryelements.h \
    ../linequadraturerules.h \
    ../mshreader.h \
    ../projectfilehandler.h \
    ../triangle.h \
    ../SolvingScript/solvingscriptreader.h \
    ../boundaryelements.h \
    ../boundaryelementsolver.h \
    ../vectortriangle.h \
    ../quadrilateral.h \
    ../nodesSection.h \
    ../SolvingScript/elementSection.h \
    ../SolvingScript/meshfileelement.h \
    ../meshfilepropertiessection.h \
    ../SolvingScript/controlsolversection.h \
    ../global.h \
    ../subdomainpropertiessection.h \
    ../node.h \
    ../vectorquadrilateral.h \
    ../robinboundarycondition.h \
    ../trianglequadraturerules.h \
    ../SolvingScript/wallimpedancesection.h \
    ../SolvingScript/wallimpedanceelement.h \
    ../gls.h \
    ../pointsource.h \
    ../SolvingScript/pressurepointssection.h \
    ../ObservationScript/observationscriptreader.h \
    ../ObservationScript/controlfieldsection.h \
    ../ObservationScript/controlspectrumsection.h \
    ../ObservationScript/fieldsection.h \
    ../ObservationScript/bespectrumsection.h \
    ../ObservationScript/bespectrumitem.h \
    ../ObservationScript/boundaryfield.h \
    ../ObservationScript/meshfield.h \
    ../ObservationScript/nodesfield.h \
    ../meshfunctions.h \
    ../ObservationScript/observationfield.h \
    ../ObservationScript/observationpoint.h \
    ../globallogstrings.h \
    ../SolvingScript/drivingsection.h \
    ../SolvingScript/drivingelement.h

FORMS += mainwindow.ui \

CONFIG += c++20 release warn_on
LIBS += -lGL -lGLU -lgomp -fopenmp -foffload=nvptx-none #-lcuda #-lOpenCL
