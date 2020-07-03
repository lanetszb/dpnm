#include <Aggregator.h>

#include <numeric>
#include <iomanip>

Aggregator::Aggregator(const std::vector<double> &propsPNM,
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
                       const double &matrixVolume) :

        equationPNM(propsPNM, throatList, throatHeight, throatLength,
                    throatWidth, connIndIn, connIndOut, poreCoordX, poreCoordY,
                    poreCoordZ, poreRadius, poreList, poreConns, connNumber,
                    porePerRow, poreLeftX, poreRightX, hydraulicCond),

        equationDiffusion(propsDiffusion),

        diffusionPartMath(propsPNM, propsDiffusion, throatList, throatHeight,
                          throatLength, throatWidth, connIndIn, connIndOut,
                          poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                          poreList, poreConns, connNumber, porePerRow,
                          poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                          matrixVolume),

        diffusionPartFlow(propsPNM, propsDiffusion, throatList, throatHeight,
                          throatLength, throatWidth, connIndIn, connIndOut,
                          poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                          poreList, poreConns, connNumber, porePerRow,
                          poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                          matrixVolume),

        couplingIniConds(propsPNM, propsDiffusion, throatList, throatHeight,
                         throatLength, throatWidth, connIndIn, connIndOut,
                         poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                         poreList, poreConns, connNumber, porePerRow,
                         poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                         matrixVolume),

        couplingMatrix(propsPNM, propsDiffusion, throatList, throatHeight,
                       throatLength, throatWidth, connIndIn, connIndOut,
                       poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                       poreList, poreConns, connNumber, porePerRow,
                       poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                       matrixVolume),

        couplingParamsOut(propsPNM, propsDiffusion, throatList, throatHeight,
                          throatLength, throatWidth, connIndIn, connIndOut,
                          poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                          poreList, poreConns, connNumber, porePerRow,
                          poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                          matrixVolume) {}

void Aggregator::calcCoupledFlow() {

    diffusionPartFlow.calcDiffPart();
    couplingMatrix.solveCoupledMatrix();
    couplingParamsOut.calcCoupledFlowParams();
}

void Aggregator::cfdProcedurePnmDiff() {

    equationPNM.setInitialCondPurePnm();
    std::vector<std::vector<int>> gammaByPressureSaved = couplingIniConds.getGamma();
    couplingIniConds.getInletFlow();

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        for (int j = 0; j < equationPNM.gammaPnm[i].size(); j++)
            equationPNM.gammaPnm[i][j] = gammaByPressureSaved[i][j];

    couplingIniConds.setInitialCondCoupledMod();
    couplingParamsOut.calcVecSum(equationPNM.networkData.poreN,
                                 equationPNM.pressure,
                                 couplingParamsOut.pressureAv,
                                 1.0 / equationPNM.networkData.poreN);
    couplingParamsOut.calcMatrixMassTot();

    couplingParamsOut.totalFlowDiff.emplace_back(0);
    couplingParamsOut.totalFlowPoresOut.emplace_back(
            equationPNM.totFlowRate * diffusionPartMath.densityConst);
    couplingParamsOut.totalFlowPoresIn.emplace_back(
            equationPNM.totFlowRate * diffusionPartMath.densityConst);

    couplingParamsOut.calcPressInlet();

    // TODO: Remove castyl
    for (double t = equationDiffusion.propsDiffusion.timeStep;
         t < equationDiffusion.propsDiffusion.time * (1. + 1.e-3);
         t += equationDiffusion.propsDiffusion.timeStep) {

        calcCoupledFlow();
        couplingParamsOut.calcPressInlet();
        couplingParamsOut.calcVecSum(equationPNM.networkData.poreN,
                                     equationPNM.pressure,
                                     couplingParamsOut.pressureAv,
                                     1.0 / equationPNM.networkData.poreN);
        couplingParamsOut.calcMatrixMassTot();
    }
}

// Getters for Python
const std::vector<double> Aggregator::getPressureAverage() const {
    return couplingParamsOut.pressureAv;
}

const std::vector<double> Aggregator::getMatrixMassTotal() const {
    return couplingParamsOut.matrixMassTotal;
}

const std::vector<double> Aggregator::getTotalFlowPoresOut() const {
    return couplingParamsOut.totalFlowPoresOut;
}

const std::vector<double> Aggregator::getTotalFlowPoresIn() const {
    return couplingParamsOut.totalFlowPoresIn;
}

const std::vector<double> Aggregator::getTotalFlowDiff() const {
    return couplingParamsOut.totalFlowDiff;
}

const std::vector<double> Aggregator::getInletPressure() const {
    return couplingParamsOut.pressIn;
}

const std::vector<double> Aggregator::getPorePressure() const {
    return equationPNM.pressure;
}