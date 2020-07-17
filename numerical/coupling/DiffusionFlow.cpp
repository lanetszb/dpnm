#include <DiffusionFlow.h>

#include <numeric>
#include <iomanip>


DiffusionFlow::DiffusionFlow(NetworkData &networkData,
                             EquationPNM &equationPNM,
                             EquationDiffusion &equationDiffusion,
                             DiffusionMath &diffusionMath,
                             IniConds &iniConds,
                             const std::vector<double> &langmuirCoeff,
                             const double &matrixVolume) :

        networkData(networkData),
        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),
        iniConds(iniConds),

        gridBlockN(equationDiffusion.propsDiffusion.gridBlockN),
        diffFlowInst(networkData.throatN, 0),
        diffFlowInstPlus(networkData.throatN, 0),
        diffFlowInstMinus(networkData.throatN, 0),
        flowDerivDiff(networkData.throatN, 0),
        dP(0) {}

void DiffusionFlow::calcDiffFlow(std::vector<double> &diffFlowVector,
                                 const double &dt) {

    for (int i = 0; i < networkData.throatN; i++) {

        for (int j = 0; j < gridBlockN; j++) {
            equationDiffusion.conc[0][j] = iniConds.matrixConc[i][j];
            equationDiffusion.conc[1][j] = iniConds.matrixConc[i][j];
        }

        equationDiffusion.cfdProcedureOneStep("default",
                                              diffusionMath.throatConc[i],
                                              networkData.throatRadius[i],
                                              diffusionMath.matrixWidth[i],
                                              networkData.throatLength[i],
                                              diffusionMath.matricesVolume[i],
                                              diffusionMath.matricesOmega[i],
                                              dt);

        // diffFlowVector[i] = equationDiffusion.flowRate / densityConst;

        double flowSum = 0;
        for (int j = 0; j < gridBlockN; j++) {
            auto conc_curr = equationDiffusion.conc[equationDiffusion.iCurr][j];
            auto conc_prev = equationDiffusion.conc[equationDiffusion.iPrev][j];

            flowSum += -1 * (conc_curr - conc_prev) *
                       equationDiffusion.localDiffusion.volCartes[j] /
                       equationDiffusion.propsDiffusion.timeStep;
        }

        diffFlowVector[i] = flowSum / diffusionMath.densityConst;

        for (int j = 0; j < gridBlockN; j++) {
            auto conc_prev = equationDiffusion.conc[equationDiffusion.iPrev][j];
            iniConds.matrixConc[i][j] = conc_prev;
        }
    }
}

void DiffusionFlow::calcDiffFlowDeriv() {

    for (int i = 0; i < networkData.throatN; i++)
        flowDerivDiff[i] = ((diffFlowInstPlus[i] - diffFlowInstMinus[i]) / dP);
}

void DiffusionFlow::updateConc() {

    for (int i = 0; i < networkData.throatN; i++) {

        for (int j = 0; j < gridBlockN; j++) {
            auto conc_curr = equationDiffusion.conc[equationDiffusion.iCurr][j];
            iniConds.matrixConc[i][j] = conc_curr;
        }
    }
}

void DiffusionFlow::calcDiffPart(const double &dt) {

    diffusionMath.calcThroatAvPress();

    // Enhance and rethink later
    dP = 0;
    diffusionMath.calcThroatConc(dP);
    calcDiffFlow(diffFlowInst, dt);

    dP = (equationPNM.pIn - equationPNM.pOut) / 10.e+5;
    diffusionMath.calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowInstPlus, dt);

    diffusionMath.calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowInstMinus, dt);

    calcDiffFlowDeriv();
    updateConc();
}