#ifndef PNFLOW_COUPLINGPARAMSOUT_H
#define PNFLOW_COUPLINGPARAMSOUT_H

#include <vector>

#include <CouplingMatrix.h>
#include <DiffusionPartFlow.h>
#include <DiffusionPartMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class CouplingParamsOut {

public:

    explicit CouplingParamsOut(const std::vector<double> &propsPNM,
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

    virtual  ~CouplingParamsOut() = default;

    CouplingMatrix couplingMatrix;
    DiffusionPartFlow diffusionPartFlow;
    DiffusionPartMath diffusionPartMath;
    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

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

#endif //PNFLOW_COUPLINGPARAMSOUT_H
