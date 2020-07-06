#ifndef PNFLOW_MATRIXSOLVER_H
#define PNFLOW_MATRIXSOLVER_H

#include <vector>

#include <DiffusionFlow.h>
#include <DiffusionMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class MatrixSolver {

public:

    MatrixSolver(EquationPNM &equationPNM,
                 EquationDiffusion &equationDiffusion,
                 DiffusionMath &diffusionMath,
                 DiffusionFlow &diffusionFlow,
                 const std::vector<double> &langmuirCoeff,
                 const double &matrixVolume);

    virtual  ~MatrixSolver() = default;

    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;
    DiffusionMath &diffusionMath;
    DiffusionFlow &diffusionFlow;

    void calcMatCoeffDiff();

    void calcMatCoupledCoeff();

    void calcCoupledFreeVector();

    void solveCoupledMatrix();

    std::vector<double> connCoeffDiff;
    std::vector<double> centralCoeffDiff;

};

#endif //PNFLOW_MATRIXSOLVER_H
