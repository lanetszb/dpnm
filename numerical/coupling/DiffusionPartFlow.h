#ifndef PNFLOW_DIFFUSIONPARTFLOW_H
#define PNFLOW_DIFFUSIONPARTFLOW_H

#include <vector>

#include <DiffusionPartMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class DiffusionPartFlow {

public:

    explicit DiffusionPartFlow(const std::vector<double> &propsPNM,
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

    virtual  ~DiffusionPartFlow() = default;

    DiffusionPartMath diffusionPartMath;
    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

    void calcDiffFlow(std::vector<double> &diffFlowVector);

    void calcDiffFlowDeriv();

    void updateConc();

    void calcDiffPart();

    std::vector<double> throatConc;
    std::vector<std::vector<double>> matrixConc;

    std::vector<double> diffFlowInst;
    std::vector<double> diffFlowInstPlus;
    std::vector<double> diffFlowInstMinus;

    double dP;
    std::vector<double> flowDerivDiff;

};
#endif //PNFLOW_DIFFUSIONPARTFLOW_H
