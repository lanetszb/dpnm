#ifndef PNFLOW_COUPLINGINICONDS_H
#define PNFLOW_COUPLINGINICONDS_H

#include <vector>

#include <DiffusionPartMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class CouplingIniConds {

public:

    explicit CouplingIniConds(const std::vector<double> &propsPNM,
                              const std::vector<double> &propsDiffusion,
                              const std::vector<int> &throatList,
                              const std::vector<double> &throatHeight,
                              const std::vector<double> &throatLength,
                              const std::vector<double> &throatWidth,
                              const std::vector<double> &connIndIn,
                              const std::vector<double> &connIndOut,
                              const std::vector<double> &poreCoordX,
                              const std::vector<double> &poreCoordY,
                              const std::vector<double> &poreCoordZ,
                              const std::vector<double> &poreRadius,
                              const std::vector<int> &poreList,
                              const std::vector<int> &poreConns,
                              const std::vector<int> &connNumber,
                              const std::vector<int> &porePerRow,
                              const std::vector<bool> &poreLeftX,
                              const std::vector<bool> &poreRightX,
                              const std::vector<double> &hydraulicCond,
                              const std::vector<double> &langmuirCoeff,
                              const double &matrixVolume);

    virtual  ~CouplingIniConds() = default;

    DiffusionPartMath diffusionPartMath;
    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

    std::vector<std::vector<double>> matrixConc;

    std::vector<std::vector<int>> getGamma();

    void getInletFlow();

    void setInitialCondCoupledMod();
};

#endif //PNFLOW_COUPLINGINICONDS_H
