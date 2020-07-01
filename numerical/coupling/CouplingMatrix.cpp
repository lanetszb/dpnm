#include <CouplingMatrix.h>

#include <numeric>
#include <iomanip>

CouplingMatrix::CouplingMatrix(const std::vector<double> &propsPNM,
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

        diffusionPartFlow(propsPNM, propsDiffusion, throatList, throatHeight,
                          throatLength, throatWidth, connIndIn, connIndOut,
                          poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                          poreList, poreConns, connNumber, porePerRow,
                          poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                          matrixVolume),

        diffusionPartMath(propsPNM, propsDiffusion, throatList, throatHeight,
                          throatLength, throatWidth, connIndIn, connIndOut,
                          poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                          poreList, poreConns, connNumber, porePerRow,
                          poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                          matrixVolume),

        equationPNM(propsPNM, throatList, throatHeight, throatLength,
                    throatWidth, connIndIn, connIndOut, poreCoordX, poreCoordY,
                    poreCoordZ, poreRadius, poreList, poreConns, connNumber,
                    porePerRow, poreLeftX, poreRightX, hydraulicCond),

        equationDiffusion(propsDiffusion) {}


void CouplingMatrix::calcMatCoeffDiff() {

    auto flowDerivDiff = diffusionPartFlow.flowDerivDiff;

    for (int i = 0; i < flowDerivDiff.size(); i++)
        connCoeffDiff[i] = 0;
    // connCoeffDiff[i] = flowDerivDiff[i] / 2;


    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++)

            coeffSum += equationPNM.gammaPnm[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] / 2;

        // TODO: understand plus or minus sign and properly name as derivDiff
        centralCoeffDiff[i] = 0;
        // centralCoeffDiff[i] = -1 * coeffSum;
        // centralCoeffDiff[i] = coeffSum;
    }
}

// potential error
void CouplingMatrix::calcMatCoupledCoeff() {

    equationPNM.calcMatCoeff();

    // for (int i = 0; i < equationPNM.connCoeff.size(); i++)
    //    equationPNM.connCoeff[i] += connCoeffDiff[i];

    for (int i = 0; i < equationPNM.centralCoeff.size(); i++)
        equationPNM.centralCoeff[i] += centralCoeffDiff[i];
}

void CouplingMatrix::calcCoupledFreeVector() {


    std::vector<double> porFlowDiff(equationPNM.networkData.poreN, 0);
    std::vector<double> porFlowDiffDer(equationPNM.networkData.poreN, 0);

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;

        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            coeffSum += equationPNM.gammaPnm[i][j] *
                        diffusionPartFlow.diffFlowInst[equationPNM.porConns[i][j]];
        }
        porFlowDiff[i] = coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            coeffSum -= equationPNM.gammaPnm[i][j] *
                        diffusionPartFlow.flowDerivDiff[equationPNM.porConns[i][j]] *
                        diffusionPartMath.throatAvPress[equationPNM.porConns[i][j]];
        }
        porFlowDiffDer[i] = 0;
        // porFlowDiffDer[i] = coeffSum;
    }

    equationPNM.calculateFreeVector("mixed",
                                    equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        // TODO: understand plus or minus sign
        if (!equationPNM.networkData.poreRightX[i])
            equationPNM.freeVector[i] += porFlowDiff[i];
    // equationPNM.freeVector[i] += porFlowDiff[i] + porFlowDiffDer[i];
    // equationPNM.freeVector[i] += -1 * (porFlowDiff[i] - porFlowDiffDer[i]);
}

void CouplingMatrix::solveCoupledMatrix() {

    calcMatCoeffDiff();
    calcMatCoupledCoeff();

    equationPNM.calculateMatrix(
            equationPNM.connCoeff,
            equationPNM.centralCoeff,
            equationPNM.networkData.poreRightX,
            equationPNM.gammaPnm,
            connCoeffDiff);

    calcCoupledFreeVector();
    equationPNM.calculateGuessVector();
    // TODO: the solution with guess vector should be preffered later
    equationPNM.calculatePress(1);
}
