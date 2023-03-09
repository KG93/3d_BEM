#ifndef HMATRIXVISUALS_H
#define HMATRIXVISUALS_H

#include <QObject>
#include <QGraphicsScene>
#include <QSvgGenerator>
#include <QFileDialog>
#include <QPainter>
#include "hmatrix.h"
#include "GUI/colorgradient.h"
#include <QImage>

/**
 * \class hMatrixVisuals
 *
 * \brief Class for visualizing H-matrices. Generates svg and png images.
 */
class hMatrixVisuals
{
public:
    hMatrixVisuals();
    static void generateSvgFromHMat(HMatrix &hMatrix); /*!< Generate SVG image of H-matrix. */
    static void matrixToGreyscaleImage(const Eigen::MatrixXcd &matrix); /*!< Generate grayscale image of H-matrix. */
    static void matrixToGreyscaleImageHigherContrast(const Eigen::MatrixXcd &matrix); /*!< Generate SVG image of H-matrix with increased contrast. */
    static void matrixToHeatMapImage(const Eigen::MatrixXcd &matrix); /*!< Generate heatmap image of H-matrix. */
    static void matrixPhaseImage(const Eigen::MatrixXcd &matrix); /*!< Generate image of the phase of the H-matrix entries. */
};

#endif // HMATRIXVISUALS_H
