#include <DiffusionPartFlow.h>

#include <numeric>
#include <iomanip>


DiffusionPartFlow::DiffusionPartFlow(const std::vector<double> &propsPNM,
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

        diffFlowInst(equationPNM.networkData.throatN, 0),
        diffFlowInstPlus(equationPNM.networkData.throatN, 0),
        diffFlowInstMinus(equationPNM.networkData.throatN, 0),
        flowDerivDiff(equationPNM.networkData.throatN, 0),
        dP(0) {}

void DiffusionPartFlow::calcDiffFlow(std::vector<double> &diffFlowVector) {

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {

        auto gridBlockN = equationDiffusion.propsDiffusion.gridBlockN;

        for (int j = 0; j < gridBlockN; j++) {
            equationDiffusion.conc[0][j] = matrixConc[i][j];
            equationDiffusion.conc[1][j] = matrixConc[i][j];
        }

        equationDiffusion.cfdProcedureOneStep(
                diffusionPartMath.throatConc[i],
                equationPNM.networkData.throatRadius[i],
                diffusionPartMath.matrixWidth[i],
                equationPNM.networkData.throatLength[i],
                diffusionPartMath.matricesVolume[i],
                diffusionPartMath.matricesOmega[i]);

        // diffFlowVector[i] = equationDiffusion.flowRate / densityConst;

        double flowSum = 0;
        for (int j = 0; j < gridBlockN; j++) {
            auto conc_curr = equationDiffusion.conc[equationDiffusion.iCurr][j];
            auto conc_prev = equationDiffusion.conc[equationDiffusion.iPrev][j];

            flowSum += -1 * (conc_curr - conc_prev) *
                       equationDiffusion.localDiffusion.volCartes[j] /
                       equationDiffusion.propsDiffusion.timeStep;
        }

        diffFlowVector[i] = flowSum / diffusionPartMath.densityConst;

        for (int j = 0; j < gridBlockN; j++) {
            auto conc_prev = equationDiffusion.conc[equationDiffusion.iPrev][j];
            matrixConc[i][j] = conc_prev;
        }
    }
}

void DiffusionPartFlow::calcDiffFlowDeriv() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        flowDerivDiff[i] = ((diffFlowInstPlus[i] - diffFlowInstMinus[i]) / dP);
}

void DiffusionPartFlow::updateConc() {

    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < throatN; i++) {

        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            auto conc_curr = equationDiffusion.conc[equationDiffusion.iCurr][j];
            matrixConc[i][j] = conc_curr;
        }
    }
}

void DiffusionPartFlow::calcDiffPart() {

    diffusionPartMath.calcThroatAvPress();

    // Enhance and rethink later
    dP = 0;
    diffusionPartMath.calcThroatConc(dP);
    calcDiffFlow(diffFlowInst);

    dP = (equationPNM.pIn - equationPNM.pOut) / 10.e+5;
    diffusionPartMath.calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowInstPlus);

    diffusionPartMath.calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowInstMinus);

    calcDiffFlowDeriv();
    updateConc();
}