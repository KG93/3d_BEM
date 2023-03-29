#ifndef BEMPARAMETERS_H
#define BEMPARAMETERS_H
/**
* \brief Simulation parameters for the BEM.
*/
struct BemParameters
{
    double frequency;
    double waveSpeed;
    double airDensity;
};

/**
* \brief Encodes different BEM coupling method to solve the non-uniqueness problem.
*/
enum BemCoupling : int {NoCoupling, BurtonMillerCoupling, KirkupCoupling};

/**
* \brief Control parameters for the BEM solver.
*/
struct BemSolverParameters
{
    bool hMatSolving = true; /*!< Use the (fast) H-matrix solver. */
    double acaRelativeError = 0.01; /*!< Controls the approximation error for the adaptive cross approximation. */
    unsigned long acaMaxRank = 0; /*!< Controls the maximum rank for the adaptive cross approximation. */
    bool usePreconditioner = false; /*!< Use HLU preconditioning for the H-matrix solver. */
    unsigned long preconditionerRank = 2; /*!< Sets the (maximum) local rank of the H-LU preconditioner for the GMRES. If set to zero, the acaRelativeError tolerance is used.*/
    double preconditionerRelError = 0; /*!< Controls the accuracy during the preconditioner calculation. */
    bool hMatFieldSolving = true; /*!< Use the (fast) H-matrix solver for the field calculation. */
    double fieldACARelativeError = 0.01; /*!< Controls the approximation error for the adaptive cross approximation in the field calculation. */
    unsigned long fieldACAMaxRank = 0; /*!< Sets the (maximum) local rank of the ACA for the field calculation. If set to zero, only the fieldACARelativeError tolerance is used.*/
    bool calculateNormAndCond = false; /*!< Controls whether to estimate the norm and condition number of the dphi bem matrix. */
    BemCoupling coupling = KirkupCoupling; /*!< Encodes different BEM coupling method to solve the non-uniqueness problem. */
};


#endif // BEMPARAMETERS_H
