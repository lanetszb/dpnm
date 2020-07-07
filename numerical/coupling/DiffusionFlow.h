#ifndef PNFLOW_DIFFUSIONFLOW_H
#define PNFLOW_DIFFUSIONFLOW_H

#include <vector>

#include <DiffusionMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>
#include <IniConds.h>

class DiffusionFlow {

public:

    DiffusionFlow(EquationPNM &equationPNM,
                  EquationDiffusion &equationDiffusion,
                  DiffusionMath &diffusionMath,
                  IniConds &iniConds,
                  const std::vector<double> &langmuirCoeff,
                  const double &matrixVolume);

    virtual  ~DiffusionFlow() = default;

    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;
    DiffusionMath &diffusionMath;
    IniConds &iniConds;

    void calcDiffFlow(std::vector<double> &diffFlowVector);

    void calcDiffFlowDeriv();

    void updateConc();

    void calcDiffPart();

    // std::vector<double> throatConc;
    // std::vector<std::vector<double>> matrixConc;

    std::vector<double> diffFlowInst;
    std::vector<double> diffFlowInstPlus;
    std::vector<double> diffFlowInstMinus;

    double dP;
    std::vector<double> flowDerivDiff;

};

#endif //PNFLOW_DIFFUSIONFLOW_H
