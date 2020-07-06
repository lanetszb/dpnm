#ifndef PNFLOW_PARAMSOUT_H
#define PNFLOW_PARAMSOUT_H

#include <vector>

#include <MatrixSolver.h>
#include <DiffusionFlow.h>
#include <DiffusionMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class ParamsOut {

public:

    explicit ParamsOut(EquationPNM &equationPNM,
                       EquationDiffusion &equationDiffusion,
                       DiffusionMath &diffusionMath,
                       DiffusionFlow &diffusionFlow,
                       MatrixSolver &matrixSolver,
                       const std::vector<double> &langmuirCoeff,
                       const double &matrixVolume);

    virtual  ~ParamsOut() = default;

    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;
    MatrixSolver &matrixSolver;
    DiffusionFlow &diffusionFlow;
    DiffusionMath &diffusionMath;

    void calcVecSum(const int &iter, const std::vector<double> &vectorToSum,
                    std::vector<double> &vectorSum, const double &mult);

    void calcMatrixMassTot();

    void calcPressInlet();

    void calcCoupledFlowParams();

    double conc_ini;

    std::vector<double> pressureAv;
    std::vector<double> pressIn;
    std::vector<double> matrixMassTotal;

    std::vector<double> totalFlowPoresOut;
    std::vector<double> totalFlowPoresIn;
    std::vector<double> totalFlowDiff;

};

#endif //PNFLOW_PARAMSOUT_H
