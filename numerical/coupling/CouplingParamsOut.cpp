#include <CouplingParamsOut.h>

#include <numeric>
#include <iomanip>

CouplingParamsOut::CouplingParamsOut(const std::vector<double> &propsPNM,
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

        couplingMatrix(propsPNM, propsDiffusion, throatList, throatHeight,
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

        equationDiffusion(propsDiffusion),

        conc_ini(equationDiffusion.propsDiffusion.concIni) {}


void CouplingParamsOut::calcMatrixMassTot() {

    double sum = 0;
    std::vector<double> matrixMass;
    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
        sum = 0;
        for (int j = 0;
             j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            sum += conc_ini *
                   equationDiffusion.localDiffusion.volCartes[j];
        }
        matrixMass.emplace_back(sum);
    }

    matrixMassTotal.emplace_back(
            accumulate(matrixMass.begin(), matrixMass.end(), 0.0));
}

void CouplingParamsOut::calcVecSum(const int &iter,
                                   const std::vector<double> &vectorToSum,
                                   std::vector<double> &vectorSum,
                                   const double &mult) {

    double sum = 0;
    for (int i = 0; i < iter; i++)
        sum += vectorToSum[i];

    vectorSum.emplace_back(sum * mult);
}

void CouplingParamsOut::calcPressInlet() {

    // TODO: Should be already known, rather than calculated every time
    double boundPoresLeftSize = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i])
            boundPoresLeftSize += 1;

    double pressInlet = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++) {
        if (equationPNM.networkData.poreLeftX[i])
            pressInlet += equationPNM.pressure[i];
    }

    pressIn.emplace_back(pressInlet / boundPoresLeftSize);
}

void CouplingParamsOut::calcCoupledFlowParams() {

    // Total flow from diffusion
    calcVecSum(equationPNM.networkData.throatN, diffusionPartFlow.diffFlowInst,
               totalFlowDiff, diffusionPartMath.densityConst);
    // diffFlowThroat += diffFlowInst[i] - flowDerivDiff[i] * throatAvPress[i];

    // equationPNM.getGammaByPressure();

    equationPNM.calcThrFlowRate();
    equationPNM.calcPorFlowRate();
    equationPNM.calcTotFlow(equationPNM.networkData.poreLeftX);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate *
                                  diffusionPartMath.densityConst);
}