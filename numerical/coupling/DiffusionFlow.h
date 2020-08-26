#ifndef PNFLOW_DIFFUSIONFLOW_H
#define PNFLOW_DIFFUSIONFLOW_H

#include <vector>

#include <EquationDiffusion.h>
#include <NetworkData.h>
#include <EquationPNM.h>
#include <DiffusionMath.h>
#include <IniConds.h>

class DiffusionFlow {

public:

    DiffusionFlow(
            const std::map<std::string, std::variant<int, double>> &paramsPnm,
            NetworkData &networkData,
            EquationPNM &equationPNM,
            EquationDiffusion &equationDiffusion,
            DiffusionMath &diffusionMath,
            IniConds &iniConds,
            const std::vector<double> &langmuirCoeff,
            const double &matrixVolume);

    virtual  ~DiffusionFlow() = default;


    NetworkData &networkData;
    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;
    DiffusionMath &diffusionMath;
    IniConds &iniConds;

    void calcDiffFlow(std::vector<double> &diffFlowVector,
                      const double &dt);

    void calcDiffFlowDeriv();

    void updateConc();

    void calcDiffPart(const double &dt);

    // std::vector<double> throatConc;
    // std::vector<std::vector<double>> matrixConc;


    std::map<std::string, std::variant<int, double>> _paramsPnm;
    std::vector<double> diffFlowInst;
    std::vector<double> diffFlowInstPlus;
    std::vector<double> diffFlowInstMinus;

    double dP;
    int &gridBlockN;
    std::vector<double> flowDerivDiff;

    std::vector<std::vector<double>> concCurr;

};

#endif //PNFLOW_DIFFUSIONFLOW_H
