#ifndef PNFLOW_INICONDS_H
#define PNFLOW_INICONDS_H

#include <vector>

#include <EquationDiffusion.h>
#include <NetworkData.h>
#include <EquationPNM.h>
#include <DiffusionMath.h>

class IniConds {

public:

    explicit IniConds(
            const std::map<std::string, std::variant<int, double>> &paramsPnm,
            NetworkData &networkData,
            EquationPNM &equationPNM,
            EquationDiffusion
            &equationDiffusion,
            DiffusionMath &diffusionMath,
            const std::vector<double> &langmuirCoeff,
            const double &matrixVolume);

    virtual  ~IniConds() = default;

    NetworkData &networkData;
    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;
    DiffusionMath &diffusionMath;

    int &gridBlockN;

    std::map<std::string, std::variant<int, double>> _paramsPnm;
    std::vector<std::vector<double>> matrixConc;

    std::vector<std::vector<int>> getGamma();

    void getInletFlow();

    void setInitialCondCoupledMod();
};

#endif //PNFLOW_INICONDS_H
