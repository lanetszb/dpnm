#ifndef PNFLOW_PARAMSOUT_H
#define PNFLOW_PARAMSOUT_H

#include <vector>

#include <EquationPNM.h>
#include <NetworkData.h>
#include <EquationDiffusion.h>
#include <DiffusionMath.h>
#include <IniConds.h>
#include <DiffusionFlow.h>
#include <MatrixSolver.h>

class ParamsOut {

public:

    explicit ParamsOut(NetworkData &networkData, EquationPNM &equationPNM,
                       EquationDiffusion &equationDiffusion,
                       DiffusionMath &diffusionPartMath, IniConds &iniConds,
                       DiffusionFlow &diffusionPartFlow, MatrixSolver &solver,
                       const std::vector<double> &langmuirCoeff,
                       const double &matrixVolume);

    virtual  ~ParamsOut() = default;

    NetworkData &networkData;
    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;
    DiffusionMath &diffusionMath;
    IniConds &iniConds;
    DiffusionFlow &diffusionFlow;
    MatrixSolver &matrixSolver;

    void calcPressureAv();

    void calcTotalFlowDiff();

    void calcMatrixMassTot();

    void calcPressInlet();

    void calcCoupledFlowParams();

    double conc_ini;
    int &gridBlockN;

    std::vector<double> pressureAv;
    std::vector<double> pressIn;
    std::vector<double> matrixMassTotal;

    std::vector<double> totalFlowPoresOut;
    std::vector<double> totalFlowPoresIn;
    std::vector<double> totalFlowDiff;

};

#endif //PNFLOW_PARAMSOUT_H
