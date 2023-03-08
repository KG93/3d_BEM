#include "hmatrixvisuals.h"

hMatrixVisuals::hMatrixVisuals()
{

}

void hMatrixVisuals::generateSvgFromHMat(HMatrix &hMatrix)
{
    double maxDimension = 200;

    long rows = hMatrix.rows();
    long cols = hMatrix.cols();
    double sizefactor = maxDimension / std::max(rows, cols);
    QGraphicsScene hMatVisual(0, 0 , sizefactor * cols, sizefactor * rows);
    QColor outlineColor = Qt::black;
//    QColor outlineColor = Qt::yellow;
//    QColor outlineColor = Qt::cyan;

    hMatrix.updatePartition();
    QVector<BlockCluster*> hMatrixPartition = hMatrix.getMinPartition();
    for(int blockIndex = 0; blockIndex < hMatrixPartition.size(); blockIndex++)
    {
        BlockCluster* block = hMatrixPartition.at(blockIndex);
        QRectF blockDimensions(sizefactor * block->colStartIndex(), sizefactor * block->rowStartIndex(), sizefactor * block->cols(), sizefactor * block->rows());
        QColor infillColor = Qt::black;
        float h;
        float s;
        float v;
        float a;
        infillColor.getHsvF(&h,&s,&v,&a);
//        infillColor = infillColor.lighter(block->compression() * 100);
        if(!block->isAdmissible)
        {
//            infillColor = Qt::yellow;
//            infillColor = Qt::cyan;
//            infillColor = Qt::darkGray;
//            infillColor = Qt::black;
        }
        else
        {
//            infillColor = Qt::black;
//            infillColor = Qt::yellow;
//            infillColor = Qt::yellow;

            infillColor.setHsvF(h, s, 1-std::min(1/block->compression(), 1.0), a);
        }
        hMatVisual.addRect(blockDimensions, QPen(outlineColor, 0.1), QBrush(infillColor));
    }

    QString filter = "Images (*.svg)";
    QString filename = QFileDialog::getSaveFileName(nullptr, "Save as .svg.","",filter, &filter);
    if ( filename.isEmpty() )
    {
        return;
    }
    if (!filename.endsWith(".svg"))
    {
        filename += ".svg";
    }


    QSvgGenerator svgGen;
//    svgGen.setResolution(50);
    svgGen.setFileName(filename);
    svgGen.setSize(QSize(sizefactor * cols, sizefactor * rows));
    svgGen.setViewBox(QRect(0, 0, sizefactor * cols, sizefactor * rows));
//    svgGen.setTitle(tr("SVG Generator Example Drawing"));
    QPainter painter;
    painter.begin(&svgGen);
    hMatVisual.render(&painter);
    painter.end();
}

void hMatrixVisuals::matrixToGreyscaleImage(const Eigen::MatrixXcd &matrix)
{
    long maxDimension = 10000;

    long rows = std::min(matrix.rows(), maxDimension);
    long cols = std::min(matrix.cols(), maxDimension);

    if(rows < 1 || cols < 1)
    {
        return;
    }
    QString filter = "Images (*.png)";
    QString filename = QFileDialog::getSaveFileName(nullptr, "Save as .png.","",filter, &filter);
    if ( filename.isEmpty() )
    {
        return;
    }
    if (!filename.endsWith(".png"))
    {
        filename += ".png";
    }
    int maxInt = INT_MAX;
    Eigen::MatrixXd absMatrix = matrix.cwiseAbs();
    absMatrix /= absMatrix.maxCoeff();

    Eigen::MatrixXi mat = (maxInt * absMatrix).cast<int>();
    QImage matrixImage(cols, rows, QImage::Format_Grayscale16);
    for(long r = 0; r < rows; r++)
    {
        for(long c = 0; c < cols; c++)
        {
            matrixImage.setPixel(c, r, 2 * mat(r, c));
        }
    }
    matrixImage.save(filename, nullptr, 100);
}

void hMatrixVisuals::matrixToGreyscaleImageHigherContrast(const Eigen::MatrixXcd &matrix)
{
    long maxDimension = 10000;

    long rows = std::min(matrix.rows(), maxDimension);
    long cols = std::min(matrix.cols(), maxDimension);

    if(rows < 1 || cols < 1)
    {
        return;
    }
    QString filter = "Images (*.png)";
    QString filename = QFileDialog::getSaveFileName(nullptr, "Save as .png.","",filter, &filter);
    if ( filename.isEmpty() )
    {
        return;
    }
    if (!filename.endsWith(".png"))
    {
        filename += ".png";
    }
    int maxInt = INT_MAX;
    Eigen::MatrixXd absMatrix = matrix.cwiseAbs().array().pow(1.9);
    absMatrix /= absMatrix.maxCoeff();

    Eigen::MatrixXi mat = (maxInt * absMatrix).cast<int>();
    QImage matrixImage(cols, rows, QImage::Format_Grayscale16);
    for(long r = 0; r < rows; r++)
    {
        for(long c = 0; c < cols; c++)
        {
            matrixImage.setPixel(c, r, 2 * mat(r, c));
        }
    }
    matrixImage.save(filename, nullptr, 100);
}

void hMatrixVisuals::matrixToHeatMapImage(const Eigen::MatrixXcd &matrix)
{
    long maxDimension = 10000;

    long rows = std::min(matrix.rows(), maxDimension);
    long cols = std::min(matrix.cols(), maxDimension);

    if(rows < 1 || cols < 1)
    {
        return;
    }
    QString filter = "Images (*.png)";
    QString filename = QFileDialog::getSaveFileName(nullptr, "Save as .png.","",filter, &filter);
    if ( filename.isEmpty() )
    {
        return;
    }
    if (!filename.endsWith(".png"))
    {
        filename += ".png";
    }
    ColorGradient heatMapGradient;
    heatMapGradient.defaultGradient();
//    Eigen::MatrixXd absMatrix = matrix.cwiseAbs()/*.array().pow(0.4)*/;
    Eigen::MatrixXd absMatrix = matrix.cwiseAbs();
    absMatrix /= absMatrix.minCoeff();
    absMatrix = matrix.cwiseAbs().array().log();

    absMatrix /= absMatrix.maxCoeff();

    float r=0, g=0, b=0;
    QColor color;

    QImage matrixImage(cols, rows, QImage::Format_RGB32);
    for(long row = 0; row < rows; row++)
    {
        for(long col = 0; col < cols; col++)
        {
            heatMapGradient.getColorAtValue(absMatrix(row, col), r, g, b);
            color.setRgbF(r,g,b);
            matrixImage.setPixelColor(col, row, color);
        }
    }
    matrixImage.save(filename, nullptr, 100);
}

void hMatrixVisuals::matrixPhaseImage(const Eigen::MatrixXcd &matrix)
{
    long maxDimension = 10000;

    long rows = std::min(matrix.rows(), maxDimension);
    long cols = std::min(matrix.cols(), maxDimension);

    if(rows < 1 || cols < 1)
    {
        return;
    }
    QString filter = "Images (*.png)";
    QString filename = QFileDialog::getSaveFileName(nullptr, "Save as .png.","",filter, &filter);
    if ( filename.isEmpty() )
    {
        return;
    }
    if (!filename.endsWith(".png"))
    {
        filename += ".png";
    }
    ColorGradient phaseGradient;
    phaseGradient.phaseGradient();
    float r=0, g=0, b=0;
    QColor color;

    QImage matrixImage(cols, rows, QImage::Format_RGB32);
    for(long row = 0; row < rows; row++)
    {
        for(long col = 0; col < cols; col++)
        {
            double phase = std::atan2(std::real(matrix(row, col)),std::imag(matrix(row, col)));
            phaseGradient.getColorAtValue(phase, r,g,b);
            color.setRgbF(r,g,b);
            matrixImage.setPixelColor(col, row, color);
        }
    }
    matrixImage.save(filename, nullptr, 100);
}
