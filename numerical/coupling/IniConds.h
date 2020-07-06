#ifndef PNFLOW_INICONDS_H
#define PNFLOW_INICONDS_H

#include <vector>

#include <DiffusionMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class IniConds {

public:

    explicit IniConds(EquationPNM &equationPNM,
                      EquationDiffusion &equationDiffusion,
                      DiffusionMath &diffusionMath,
                      const std::vector<double> &langmuirCoeff,
                      const double &matrixVolume);

    virtual  ~IniConds() = default;

    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;
    DiffusionMath &diffusionMath;

    std::vector<std::vector<double>> matrixConc;

    std::vector<std::vector<int>> getGamma();

    void getInletFlow();

    void setInitialCondCoupledMod();
};

#endif //PNFLOW_INICONDS_H
