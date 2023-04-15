//The  boundary element solver is implemented after the text
//"the boundary element method in acoustics" by Stephen Kirkup
//1997/2007

#ifndef BOUNDARYELEMENTSOLVER_H
#define BOUNDARYELEMENTSOLVER_H

#include "boundaryelements.h"
//#include "LinearElements/linearboundaryelements.h"
#include "trianglequadraturerules.h"
#include "linequadraturerules.h"
#include "gls.h"
#include "pointsource.h"
#include "ObservationScript/observationpoint.h"
#include "ObservationScript/observationfield.h"
#include "globallogstrings.h"
#include <eigen3/Eigen/Geometry>
#include "Timer.h"
#include "bemparameters.h"

#include "HMatrix/harithm.h"
#include "HMatrix/gmres.h"

#include <malloc.h>
#include <functional>
//#include "HMatrix/hmatrixvisuals.h"

/**
 * \brief The boundary element solver class.
 *
 * This class provides methods to solve the Helmholtz equation using the boundary
 * element method. The solver can use H-matrices to speed up the computation
 * and provides the option of H-LU preconditioning.
 *
 * Adaptive cross approximation can be used to approximate the BEM matrices via low-rank approximations and thereby reducing the computational
 * complexity from \f$ O(n^2) \f$ to \f$ O(n*\log(n)*k^2) \f$, where \f$ O(n^1) \f$ is the number of boundary elements
 * and \f$ k \f$ is the local rank of the approximation. The local rank is \f$ O(\sqrt{\kappa} * \log(1/relError)) \f$, where \f$ \kappa \f$ is the wavenumber and \f$ relError \f$ is the relative approximation error of the H-matrix.
 * The H-matrix method is therefore only efficient for low to medium frequencies.
 */
class BoundaryElementSolver : public QObject
{
    Q_OBJECT
public:
    BoundaryElementSolver();

    /**
    * \brief Get the boundary elements from the solver.
    * \return The boundary elements.
    */
    BoundaryElements getBoundaryElements(){return boundaryElements;}

    /**
    * \brief Set the boundary elements for the solver.
    * \param elements The boundary elements.
    */
    void setBoundaryElements(BoundaryElements elements){this->boundaryElements = elements;}

    /**
    * \brief Set the wavenumber for the solver.
    * The wavenumber is a complex number that is used to solve the Helmholtz equation. It is defined as the product of the wave speed and the frequency, and is given by the formula
    * \f$ k = \omega \cdot c \f$
    * where \f$ \omega \f$ is the angular frequency and \f$ c \f$ is the wave speed.
    * \param wavenumber The wavenumber to set.
    */
    void setWavenumber(std::complex<double> wavenumber){this->wavenumber = wavenumber;}

    /**
    * Set the wavenumber for the solver.
    * \param wavenumber The complex wavenumber to use in the solver.
    */
    void setFrequency(double frequency){this->frequency = frequency;}

    /**
    * Set the wavespeed for the solver.
    * \param waveSpeed The speed of the waves to use in the solver.
    */
    void setWavespeed(double waveSpeed){this->waveSpeed = waveSpeed;}

    /**
    * Set the air density for the solver.
    * \param airDensity The density of the air to use in the solver.
    */
    void setAirDensity(double airDensity){this->airDensity = airDensity;}

    /**
    * \brief Set all the (environment) parameters for the BEM solver.
    * \param bemParameters The collection of control parameters for the BEM solver.
    */
    void setBemParameters(BemParameters bemParameters){this->frequency = bemParameters.frequency;
                                                       this->waveSpeed = bemParameters.waveSpeed;
                                                       this->airDensity = bemParameters.airDensity;}

    /**
    \brief Set the point sources for the solver.
    \param newPointSources A QVector of PointSource objects representing the point sources.
    */
    void setPointSources(QVector<PointSource> newPointSources){this->pointSources = newPointSources;}

    /**
    * \brief Set observation points.
    * \param observationPoints The observation points to set.
    */
    void setObservationPoints(QVector<ObservationPoint> observationPoints){this->observationPoints = observationPoints;}

    /**
    * \brief Set observation fields.
    * \param observationFields The observation fields to set.
    */
    void setObservationFields(QVector<ObservationField> observationFields){this->observationFields = observationFields;}

    /**
    * \brief Set the order of the Gauss quadrature rule with high accuracy requirements
    * Set the order of the Gauss quadrature rule that will be used by the solver for integrations with high accuracy requirements.
    * \param order The order of the Gauss quadrature rule
    */
    void setHighOrderQuadratureRule(int order){this->orderHigh = order;}

    /**
    * \brief Set the order of the Gauss quadrature rule without high accuracy requirements
    * Set the order of the Gauss quadrature rule that will be used by the solver for integrations without high accuracy requirements.
    * \param order The order of the Gauss quadrature rule
    */
    void setRegularOrderQuadratureRule(int order){this->orderRegular = order;}

    /**
    * \brief Get whether the solver uses H-arithmatic.
    * \return True if the solver uses H-arithmatic.
    */
    bool getHSolving(){return bemSolverParameters.hMatSolving;}

    /**
    * \brief Set whether to use H-arithmetic in the solver.
    * \param val True to use H-arithmetic.
    */
    void setHSolving(bool val){this->bemSolverParameters.hMatSolving = val;}

    /**
    * \brief Get the relative approximation error of the adaptive cross approximation.
    * \return The relative approximation error of the adaptive cross approximation.
    */
    double getACArelativeError(){return bemSolverParameters.acaRelativeError;}

    /**
    * \brief Set the relative approximation error of the adaptive cross approximation.
    * \param relError The relative approximation error.
    */
    void setACArelativeError(double relError){this->bemSolverParameters.acaRelativeError = std::max(relError, 0.0);}

    /**
    * \brief Set the maximum rank of the adaptive cross approximation.
    * \param maxRank The maximum rank of the adaptive cross approximation.
    */
    void setACAMaxRank(long maxRank){this->bemSolverParameters.acaMaxRank = maxRank;}

    /**
    * \brief Set whether to use HLU preconditioning in the HMatrix solver.
    * \param val True to use HLU preconditioning.
    */
    void setUsePreconditioner(bool val){this->bemSolverParameters.usePreconditioner = val;}

    /**
    * \brief Get the (maximum) rank of the H-matrix preconditioner.
    * \return The (maximum) rank of the H-matrix preconditioner.
    */
    unsigned long getPreconditionerRank(){return bemSolverParameters.preconditionerRank;}

    /**
    * \brief Set the (maximum) rank of the H-matrix preconditioner.
    * \param rank The (maximum) rank of the H-matrix preconditioner.
    */
    void setPreconditionerRank(unsigned long rank){this->bemSolverParameters.preconditionerRank = rank;}

    /**
    * \brief Set the relative error of the H-matrix preconditioner.
    * \param relError The relative error of the H-matrix preconditioner.
    */
    void setPreconditionerRelativeError(double relError){this->bemSolverParameters.preconditionerRelError = std::max(relError, 0.0);}

    /**
    * \brief Get whether the solver uses H-arithmatic for sound field calculation.
    * \return True, if the solver uses H-arithmatic for sound field calculation.
    */
    bool getHFieldSolving(){return bemSolverParameters.hMatFieldSolving;}

    /**
    * \brief Set to true to use H-arithmatic in the solver for sound field calculation.
    * \param val Boolean value indicating whether to use H-arithmatic.
    */
    void setHFieldSolving(bool val){this->bemSolverParameters.hMatFieldSolving = val;}

    /**
    * \brief Get the relative approximation error of the adaptive cross approximation for sound field calculation.
    * \return The relative approximation error of the adaptive cross approximation for sound field calculation.
    */
    double getFieldACARelError(){return bemSolverParameters.fieldACARelativeError;}

    /**
    * \brief Set the relative approximation error of the adaptive cross approximation for sound field calculation.
    * \param relError The relative approximation error.
    */
    void setFieldACARelError(double relError){this->bemSolverParameters.fieldACARelativeError = std::max(relError, 0.0);}

    /**
    * \brief Get the maximum rank of the adaptive cross approximation for sound field calculation.
    * \return The maximum rank of the adaptive cross approximation for sound field calculation.
    */
    unsigned long getFieldACAMaxRank(){return bemSolverParameters.fieldACAMaxRank;}

    /**
    * \brief Set the maximum rank of the adaptive cross approximation for sound field calculation.
    * \param rank The maximum rank.
    */
    void setFieldACAMaxRank(unsigned long rank){this->bemSolverParameters.fieldACAMaxRank = rank;}

    /**
    * \brief Control whether to calculate the matrix norm and condition number during bem solving.
    * \param calcNorm Boolean value indicating whether to calculate the norm and condition number.
    */
    void setCalculateNormAndConditionNumber(bool calcNorm){this->bemSolverParameters.calculateNormAndCond = calcNorm;}

    /**
    * \brief Set all the control parameters for the BEM solver.
    * \param bemSolverParameters The collection of control parameters for the BEM solver.
    */
    void setBemSolverParameters(BemSolverParameters bemSolverParameters){this->bemSolverParameters = bemSolverParameters;}

    /**
    * \brief Return the observation points from the solver.
    * \return The observation points from the solver.
    */
    QVector<ObservationPoint> getObservationPoints(){return observationPoints;}

    /**
    * \brief Return the observation fields from the solver.
    * \return The observation fields from the solver.
    */
    QVector<ObservationField> getObservationFields(){return observationFields;}

    /**
    * \brief Prepare the boundary elements for the solver.
    */
    void prepareBoundaryElements();

    /**
    * \brief Start the regular solver with LU factorization.
    */
    void regularSolveWithLU();

    /**
    * \brief Start the regular solver with GMRES.
    * \param error The relative residual error of the GMRES solver.
    */
    void regularSolveWithGMRES(double error);


//    /**
//    * \brief Start the regular solver with LU factorization.
//    */
//    void regularSolveWithLULinear();


    /**
    * \brief Start the bem solver with H-matrix acceleration.
    */
    void hMatrixSolve();

    /**
    * \brief Calculate the truncated sum of reflection matrices for two parallel impedance planes.
    * The higher terms correspond to higher order reflections between the planes.
    * Unfinished (small memory leak).
    * \param reflMatPhi The reflection matrix for the potential.
    * \param reflMatDPhi The reflection matrix for the normal derivative of the potential.
    * \param maxRank The maximum rank of the adaptive cross approximation.
    * \param relativeError The relative approximation error of the adaptive cross approximation.
    * \param elements The boundary elements.
    * \param reflElements The reflected boundary elements.
    * \param clusterTree The cluster tree.
    * \param lastPlaneIndex The last plane index. Default is -1.
    */
    void calcReflectionMatrices(HMatrix &reflMatPhi, HMatrix &reflMatDPhi, const long maxRank, const double relativeError, const BoundaryElements &elements, BoundaryElements reflElements, std::shared_ptr<ClusterTree> clusterTree, int lastPlaneIndex = -1);

    /**
    * \brief Group the block-corresponding boundary elements by their normals' orientation and return a cross product of the groups' representatives by themselves.
    * This function groups the boundary elements corresponding to the given block by their normal's orientation and returns a cross product of the group representatives by themselves.
    * The returned pivot elements can be used in the partialPivotACAgivenPivotIndices() function.
    * \param[in] block BlockCluster corresponding to the boundary elements.
    * \param[in] filterConeAngle Angle (in radians) defining the cone within which the normals are considered to be oriented in the same direction.
    * \return Vector containing pairs of indices corresponding to the pivot elements.
    */
    QVector<std::pair<long,long>> getNormalFilteredPivotIndices(BlockCluster* block, double filterConeAngle = 0.44);

    /**
    * This function calculates a row vector of an H-matrix block by using the given implicit matrix function. It is used in the ACA algorithms for the low-rank assembly of the H-matrix blocks.
    * \param rowVector The row vector to be calculated.
    * \param rowIndex The index of the row of the H-matrix block.
    * \param columnStartIndex The start index of the columns of the H-matrix block.
    * \param implicitMatrix The implicit matrix function that is used to evaluate the elements of the row vector.
    */
    void calcHBlockRowVector(Eigen::RowVectorXcd &rowVector, const long rowIndex, const long columnStartIndex, std::function<std::complex<double> (long, long)> implicitMatrix);

    /**
    * This function calculates a column vector of an H-matrix block by using the given implicit matrix function. It is used in the ACA algorithms for the low-rank assembly of the H-matrix blocks.
    * \param columnVector The column vector to be calculated.
    * \param rowStartIndex The start index of the rows of the H-matrix block.
    * \param columnIndex The index of the column of the H-matrix block.
    * \param implicitMatrix The implicit matrix function that is used to evaluate the elements of the column vector.
    */
    void calcHBlockColumnVector(Eigen::VectorXcd &columnVector, const long rowStartIndex, const long columnIndex, std::function<std::complex<double> (long, long)> implicitMatrix);


    /**
    * This function returns the mathematical constant \f$ \pi \f$.
    * \return The value of \f$ \pi \f$.
    */
    double getPI() {return global::PI;}

    /**
    * This function returns a tolerance value that is used in the boundary element solver. The value should be in the magnitude of the rounding error.
    * \return A tolerance value.
    */
    double getTiny() {return global::tiny;}

    /**
    * This function returns a reference to the boundary elements used in the boundary element solver.
    * \return A reference to the boundary elements.
    */
    BoundaryElements* getBoundaryElementsReference() {return &boundaryElements;}

    /**
    * \brief Start the solver for the field calculation.
    * The function will either call calculateFieldSolutionRegular() or calculateFieldSolutionFast() (with H-matrix acceleration).
    */
    void calculateFieldSolution();

    void phiSolutionToSoundPressure(); /*!< \brief Calculate the sound pressure from the acoustic potential (phi). */

    Eigen::VectorXcd directRightHandSide;
    Eigen::VectorXcd HMatRightHandSide;
    Eigen::VectorXcd testSol1;
    Eigen::VectorXcd testSolLU;
    Eigen::VectorXcd testSolRandGMRES;
    Eigen::VectorXcd testSol2;
    Eigen::MatrixXcd testDirect;
    Eigen::MatrixXcd testDPhiDirect;
    Eigen::MatrixXcd testPhiHMat;
    Eigen::MatrixXcd testDPhiHMat;

private:
    /**
    * \brief Calculate the sound field by full assembly of the interaction matrices.
    */
    void calculateFieldSolutionRegular();

    /**
    * \brief Accelerated calculation of the sound field by exploiting a clustering stategy / h-matrix technique.
    * \param relativeError The relative error used for the H-matrix compression.
    * \param maxRank The maximum rank used for the H-matrix compression.
    */
    void calculateFieldSolutionFast(const double relativeError, const long maxRank);

    /**
    * \brief Calculate the truncated sum of reflection matrices for two parallel impedance planes. The higher terms correspond to higher order reflections between the planes.
    * \param reflMatPhi the reflection matrix for the phi BEM operator
    * \param reflMatDPhi the reflection matrix for the dPhi BEM operator
    * \param maxRank the maximum rank of the reflection matrices
    * \param relativeError the relative error used in the truncation of the reflection matrices
    * \param observationPoints the field observation points
    * \param elements the boundary elements
    * \param obsClusterTree the cluster tree for the field observation points
    * \param elementsClusterTree the cluster tree for the boundary elements
    * \param lastPlaneIndex the index of the last reflection plane (default -1, i.e. no plane)
    */
    void calcReflectionMatricesField(HMatrix &reflMatPhi, HMatrix &reflMatDPhi, const long maxRank, const double relativeError, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &elements, std::shared_ptr<ClusterTree> obsClusterTree, std::shared_ptr<ClusterTree> elementsClusterTree, int lastPlaneIndex = -1);

    /**
    * \brief Calculate the BEM operators for a field observation point and a boundary element.
    * \param observationPoint - The location of the field observation point.
    * \param boundaryTriangleIndex - The index of the boundary element triangle.
    * \param Mk - The value of the BEM operator Mk.
    * \param Lk - The value of the BEM operator Lk.
    */
    void BemOperatorField(const Eigen::Vector3d observationPoint, const int boundaryTriangleIndex, std::complex<double> &Mk, std::complex<double> &Lk, const BoundaryElements &boundaryElements);

    /**
    * \brief Calculate the phi BEM operators for a field observation point and a boundary element.
    * \param[in] rowIndex The index of the observation point for which the BEM operator is calculated.
    * \param[in] columnIndex The index of the boundary element for which the BEM operator is calculated.
    * \param[in] observationPoints List of all observation points.
    * \param[in] boundaryElements List of all boundary elements.
    * \return The value of the phi BEM operator for the given observation point and boundary triangle.
    */
    std::complex<double> implicitPhiMatrixField(const long rowIndex, const long columnIndex, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &boundaryElements);

    /**
    * \brief Calculate the dphi BEM operators for a field observation point and a boundary element.
    * \param rowIndex The index of the observation point for which the BEM operator is calculated.
    * \param columnIndex The index of the boundary element for which the BEM operator is calculated.
    * \param observationPoints List of all observation points.
    * \param boundaryElements List of all boundary elements.
    * \return The dphi BEM operator for the given observation point and boundary element.
    */
    std::complex<double> implicitDPhiMatrixField(const long rowIndex, const long columnIndex, const QVector<Eigen::Vector3d> &observationPoints, const BoundaryElements &boundaryElements);

    /**
    * \brief Calculates the pivot indices for the ACA with extra accuracy. The approximation accuracy of the ACA will be tested on all pivot indices.
    * \param block A pointer to the hmatrix block for which the pivot indices will be calculated
    * \param filterConeAngle The angle of the filter cone used to filter the pivot indices
    * \return QVector<std::pair<long,long>> A vector of pivot indices for the given block
    */
    QVector<std::pair<long,long>> getNormalFilteredPivotIndicesForField(BlockCluster* block, const double filterConeAngle = 0.44);

    /**
    * \brief Set up the gauss-puadrature weights and abscissas for the prescribed quadrature orders.
    */
    void setUpQuadratureRules();

    /**
    * \brief Boundary elements of the BEM simulation.
    */
    BoundaryElements boundaryElements;

    /**
    * \brief List of observation points.
    */
    QVector<ObservationPoint> observationPoints;

    /**
    * \brief List of observation fields.
    */
    QVector<ObservationField> observationFields;

    /**
    * \brief Wave speed for the BEM simulation.
    */
    double waveSpeed;

    /**
    * \brief The frequencies for the simulation.
    */
    double frequency;

    /**
    * \brief Wavenumber for the BEM simulation.
    */
    std::complex<double> wavenumber;
    std::complex<double> wavenumberSquared;
    std::complex<double> wavenumberSquaredHalf;
    std::complex<double> iWavenumber;

    /**
    * \brief Coupling parameter for the BEM simulation.
    */
    std::complex<double> couplingParameter;
    std::complex<double> couplingParameterHalf;

    /**
    * \brief Vector of source terms for the BEM simulation.
    */
    Eigen::VectorXcd sourceTermVector;

    // alpha * Phi + beta * dPhi = f
    /**
    * \brief Vector of alpha-coefficients from the robin boundary conditions for all the triangles.
    */
    Eigen::VectorXcd alpha;
    /**
    * \brief Vector of beta-coefficients from the robin boundary conditions for all the triangles.
    */
    Eigen::VectorXcd beta;
    /**
    * \brief Vector of right-hand sides from the robin boundary conditions for all the triangles.
    */
    Eigen::VectorXcd f;
    /**
    * \brief Stores information, whether the implicit substitutions in the BEM matrics will be made for phi or for dphi.
    */
    Eigen::VectorXi substituteDPhiWithPhi;

    /**
    * \brief List of point sources for the simulation.
    */
    QVector<PointSource> pointSources;

    /**
    * \brief Calculate the collocation BEM operators for constant triangles.
    * \param row The row index for which the BEM operators are calculated. Corresponds to the collocation point \mathbf{p}.
    * \param column The column index for which the BEM operators are calculated. Corresponds to the triangle over which is integrated.
    * \param Lk The resulting Lk operator \f$ (G(\mathbf{p},\mathbf{q})) \f$ will be stored in this variable.
    * \param Mk The resulting Mk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
    * \param Mtk The resulting Mtk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
    * \param Nk The resulting Nk \frac{\partial^2 G(\mathbf{p},\mathbf{q}) }{\partial n_q \partial n_p} operator will be stored in this variable.
    */
    inline void BemOperatorsConst(const long row, const long column, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk);

    /**
    * Calculate the singular collocation BEM operators for constant triangles. Used in the ordinary (non-accalerated) BEM and the full matrices of the H-BEM. Singularities are treated by a polar integration method.
    * \param triangleIndex The index of the triangle to calculate the operators for.
    * \param Lk The resulting Lk operator \f$ (G(p,q)) \f$ will be stored in this variable.
    * \param Mk The resulting Mk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
    * \param Mtk The resulting Mtk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
    * \param Nk The resulting Nk \frac{\partial^2 G(\mathbf{p},\mathbf{q}) }{\partial n_q \partial n_p} operator will be stored in this variable.
    */
    inline void BemOperatorsSingularPolarInt(const long triangleIndex, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk);

    /**
    * Calculate the singular collocation BEM operators for constant triangles. Used in the ordinary (non-accalerated) BEM and the full matrices of the H-BEM. Singularities are treated by sigularity subtraction.
    * \param row The index of the row of the operator matrix.
    * \param column The index of the column of the operator matrix.
    * \param Lk The resulting Lk operator \f$ (G(p,q)) \f$ will be stored in this variable.
    * \param Mk The resulting Mk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
    * \param Mtk The resulting Mtk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
    * \param Nk The resulting Nk \frac{\partial^2 G(\mathbf{p},\mathbf{q}) }{\partial n_q \partial n_p} operator will be stored in this variable.
    */
    inline void BemOperatorsSingularitySubtraction(const long row, const long column, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk);

    /**
    * \brief Calculates the nearly-singular BEM operators for collocation points near the singularity for constant triangles.
    * This function is used for evaluating the off-diagonal elements of the H-BEM matrices.
    * (near) Singularities are treated by splitting the source triangle into three sub-triangles and using polar integration with sinh transformation.
    * \param row The index of the observation triangle in the list of triangles.
    * \param column The index of the source triangle in the list of triangles.
    * \param Lk The resulting Lk operator will be stored in this variable.
    * \param Mk The resulting Mk operator will be stored in this variable.
    * \param Mtk The resulting Mtk operator will be stored in this variable.
    * \param Nk The resulting Nk operator will be stored in this variable.
    */
    inline void BemOperatorsNearSingSinh(const long row, const long column, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk);

    /**
    * \brief Calculates the nearly-singular BEM operators for collocation points near the singularity for constant triangles.
    * This method is used on the off-diagonal of the full matrices of the H-BEM. It splits the source triangle into three triangles and uses polar integration to evaluate the BEM operators.
    * \param row The index of the observation point triangle in the boundaryElements vector.
    * \param column The index of the source triangle in the boundaryElements vector.
    * \param Lk The resulting L_k operator will be written to this variable.
    * \param Mk The resulting M_k operator will be written to this variable.
    * \param Mtk The resulting M^T_k operator will be written to this variable.
    * \param Nk The resulting N_k operator will be written to this variable.
    */
    inline void BemOperatorsNearSing(const long row, const long column, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk);

    /**
    * \brief Calculate the collocation BEM operators for a reflected geometry. No special singularity treatment.
    * This function calculates the boundary element method (BEM) operators Lk, Mk, Mtk, and Nk for a given row and column index and a given reflected geometry. The reflected geometry is defined by the reflectedElements and the collocation points lie on the midpoints of boundaryElements.
    * \param row The row index.
    * \param column The column index.
    * \param Lk The Lk operator. This will be overwritten by the function.
    * \param Mk The Mk operator. This will be overwritten by the function.
    * \param Mtk The Mtk operator. This will be overwritten by the function.
    * \param Nk The Nk operator. This will be overwritten by the function.
    * \param boundaryElements The boundary elements.
    * \param reflectedElements The reflected boundary elements.
    */
    inline void BemOperatorsReflected(const long row, const long column, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements);

    /**
    *  \brief Calculates the non-singular BEM operators for the reflection phi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \param boundaryElements The boundary elements for the original geometry
    * \param reflectedElements The boundary elements for the reflected geometry
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitPhiMatrixReflected(const long rowIndex, const long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements);

    /**
    * \brief Calculates the non-singular BEM operators for the reflection dphi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \param boundaryElements The boundary elements for the original geometry
    * \param reflectedElements The boundary elements for the reflected geometry
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitDPhiMatrixReflected(const long rowIndex, const long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements);

    /**
    * \brief Calculates the non-singular BEM operators for the reflection phi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitPhiMatrix(long rowIndex, long columnIndex);

    /**
    * \brief Calculates the non-singular BEM operators for the dphi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitDPhiMatrix(long rowIndex, long columnIndex);

    /**
    * \brief Calculates the non-singular BEM operators for the reflection phi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitPhiMatrix2(long rowIndex, long columnIndex);

    /**
    * \brief Calculates the non-singular BEM operators for the dphi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitDPhiMatrix2(long rowIndex, long columnIndex);

    /**
    * \brief Calculates the non-singular BEM operators for the reflection phi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitPhiMatrixPlusReflected(long rowIndex, long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements);

    /**
    * \brief Calculates the non-singular BEM operators for the dphi-matrix. Used in the ACA method.
    * \param rowIndex Index of the row in the BEM matrix
    * \param columnIndex Index of the column in the BEM matrix
    * \return The calculated value for the BEM operator
    */
    std::complex<double> implicitDPhiMatrixPlusReflected(long rowIndex, long columnIndex, const BoundaryElements &boundaryElements, const BoundaryElements &reflectedElements);

    ////////////////////////////////////linear Elements//////////////////////////////////

//    /**
//    * \brief Calculate the collocation BEM operators for linear triangles.
//    * \param row The row index for which the BEM operators are calculated. Corresponds to the collocation point \mathbf{p}.
//    * \param column The column index for which the BEM operators are calculated. Corresponds to the triangle over which is integrated.
//    * \param Lk The resulting Lk operator \f$ (G(\mathbf{p},\mathbf{q})) \f$ will be stored in this variable.
//    * \param Mk The resulting Mk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
//    * \param Mtk The resulting Mtk operator \f$ \frac{\partial G(\mathbf{p},\mathbf{q}) }{\partial n_q} \f$ will be stored in this variable.
//    * \param Nk The resulting Nk \frac{\partial^2 G(\mathbf{p},\mathbf{q}) }{\partial n_q \partial n_p} operator will be stored in this variable.
//    */
//    void BemOperatorsLinear(const long collocationPointIndex, const long domainIndex, Eigen::Vector3d shapeFuncNodeWeights, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk, const LinearBoundaryElements &elements, const bool singularity);
//    inline void BemOperatorsLinearNode1(const long collocationPointIndex, const long domainIndex, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk, const LinearBoundaryElements &elements, const bool singularity);
//    inline void BemOperatorsLinearNode2(const long collocationPointIndex, const long domainIndex, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk, const LinearBoundaryElements &elements, const bool singularity);
//    inline void BemOperatorsLinearNode3(const long collocationPointIndex, const long domainIndex, std::complex<double> &Lk, std::complex<double> &Mk, std::complex<double> &Mtk, std::complex<double> &Nk, const LinearBoundaryElements &elements, const bool singularity);


    void calculateBoundaryConditionAndSourceTermVector(); /*!< Set up the boundary conditions and the resulting implicit substitutions. Also calculate the source term vector. */
    inline std::complex<double> sourceTerm(const PointSource source, const Eigen::Vector3d listeningPosition, const Eigen::Vector3d normal); /*!< Calculate the external sound field source term at an observation point. */
    std::complex<double> sourcePhiTerm(const PointSource source, const Eigen::Vector3d listeningPosition); /*!< Calculate only the phi term os the external sound field source at an observation point. */
    void setCouplingParameterNegative(); /*!< Set the BEM coupling parameter negative. */

    static constexpr double PI2 = 2.0 * global::PI;
    static constexpr double PI4 = 4.0 * global::PI;
    static constexpr std::complex<double> imaginaryUnit = std::complex<double>(0.0, 1.0);
    double staticPressure = 101325; //Pa at 20°C and sea level
    double airDensity = 1.2041;
    QRandomGenerator randomGenerator;

    // Set up the order of the quadrature methods
    int orderLow = 5;
    Eigen::MatrixX3d lowOrderweightsAndAbscissa = triangleQuadratureRules::weightsandAbscissa(orderLow);
    int orderRegular = 8;
    Eigen::MatrixX3d regularWeightsAndAbscissa = triangleQuadratureRules::weightsandAbscissa(orderRegular);
    int orderHigh = 8;
    Eigen::MatrixX3d highOrderweightsAndAbscissa = triangleQuadratureRules::weightsandAbscissa(orderHigh);
    int numberOfLowOrderWeightsAndAbscissa = lowOrderweightsAndAbscissa.rows();
    int numberOfRegularWeightsAndAbscissa = regularWeightsAndAbscissa.rows();
    int numberOfHighOrderWeightsAndAbscissa = highOrderweightsAndAbscissa.rows();

    int orderLineQuadrature = 5;
    Eigen::MatrixX2d weightsAndAbscissaLineQuadrature = lineQuadratureRules::weightsandAbscissa(orderLineQuadrature);
    int numberOfWeightsAndAbscissaLineQuadrature = weightsAndAbscissaLineQuadrature.rows();

    BemSolverParameters bemSolverParameters;  /*!< \brief Control parameters for the BEM solver. */

    static constexpr bool polarIntegOfSingOps = true; /*!< Use polar integration for the calculation of singular BEM operators in the regular BEM. */
    static constexpr bool avoidHMatGLS = true;
    static constexpr bool avoidGLS = true; /*!< Assemble the BEM matrices with implicit substitution and avoid the GLS routine. */
    static constexpr bool testACA = false; ///// set true only for testing; stops the solve routines after matrix assembly

//    bool negativeCoupling = false;

signals:
    void updateLog(); /*!< Update the GUI log. */

public slots:
    /**
    * \brief Start the boundary element solver.
    * The function will either call regularSolve(), regularSolveWithGMRES() or hMatrixSolve().
    */
    void calculateBoundarySolution();
    void setCoupling(BemCoupling coupling){bemSolverParameters.coupling = coupling;} /*!< Sets the different BEM coupling method to solve the non-uniqueness problem. */
};

#endif // BOUNDARYELEMENTSOLVER_H
