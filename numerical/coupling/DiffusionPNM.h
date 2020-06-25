#ifndef PNFLOW_DIFFUSIONPNM_H
#define PNFLOW_DIFFUSIONPNM_H

#include <vector>

#include <EquationDiffusion.h>
#include <EquationPNM.h>


class DiffusionPNM {

public:

    explicit DiffusionPNM(
            const std::vector<double> &propsPNM,
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
            const std::vector<double> &langmuirCoeff);

    virtual  ~DiffusionPNM() = default;

    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

    std::vector<std::vector<int>> getGamma();

    void getInletFlow();

    double calcSideLength(std::vector<double> &poreCoord);

    double calcDensConst();

    void calcRockVolume();

    void calcEffRadius();

    void calcMatrixWidth();

    void calcMatricesOmega();

    void calcMatricesVolume();

    double calcLangmConc(double pressure);

    void calcThroatConc(const double &dP);

    void calcThroatAvPress();

    void calcVecSum(const int &iter, const std::vector<double> &vectorToSum,
                    std::vector<double> &vectorSum, const double &mult);

    void calcMatrixMassTot();

    void calcPressInlet();

    void calcDiffFlow(std::vector<double> &diffFlowVector);

    void calcDiffFlowDeriv();

    void calcMatCoeffDiff();

    void cfdProcedureDiff();

    void calcDiffPart();

    void solveCoupledMatrix();

    void calcMatCoupledCoeff();

    void calcCoupledFlowParams();

    void calcCoupledFreeVector();

    void updateConc();

    void setInitialCond();

    void calcCoupledFlow();

    // Getters for Python

    const std::vector<double> getPressureAverage() const;

    const std::vector<double> getMatrixMassTotal() const;

    const std::vector<double> getTotalFlowPoresOut() const;

    const std::vector<double> getTotalFlowPoresIn() const;

    const std::vector<double> getTotalFlowDiff() const;

    const std::vector<double> getPorePressure() const;

    const std::vector<double> getInletPressure() const;


    std::vector<double> langmuirCoeff;

    double rockVolume;
    double langmConc;
    double conc_ini;
    double densityConst;


    std::vector<double> effRadius;
    std::vector<double> matrixWidth;

    std::vector<double> throatAvPress;
    std::vector<double> throatConc;

    std::vector<std::vector<double>> matrixConc;
    std::vector<std::vector<double>> matricesOmega;
    std::vector<std::vector<double>> matricesVolume;


    std::vector<double> diffFlowInst;
    std::vector<double> diffFlowInstPlus;
    std::vector<double> diffFlowInstMinus;

    double dP;

    std::vector<double> flowDerivDiff;

    std::vector<double> connCoeffDiff;
    std::vector<double> centralCoeffDiff;

    std::vector<double> pressureAv;
    std::vector<double> pressIn;
    std::vector<double> matrixMassTotal;

    std::vector<double> totalFlowPoresOut;
    std::vector<double> totalFlowPoresIn;
    std::vector<double> totalFlowDiff;
};

#endif
