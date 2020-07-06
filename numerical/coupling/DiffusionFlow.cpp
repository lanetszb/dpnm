#include <DiffusionFlow.h>

#include <numeric>
#include <iomanip>


DiffusionFlow::DiffusionFlow(EquationPNM &equationPNM,
                             EquationDiffusion &equationDiffusion,
                             DiffusionMath &diffusionMath,
                             IniConds &iniConds,
                             const std::vector<double> &langmuirCoeff,
                             const double &matrixVolume) :

        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),
        iniConds(iniConds),

        diffFlowInst(equationPNM.networkData.throatN, 0),
        diffFlowInstPlus(equationPNM.networkData.throatN, 0),
        diffFlowInstMinus(equationPNM.networkData.throatN, 0),
        flowDerivDiff(equationPNM.networkData.throatN, 0),
        dP(0) {}

void DiffusionFlow::calcDiffFlow(std::vector<double> &diffFlowVector) {

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {

        auto gridBlockN = equationDiffusion.propsDiffusion.gridBlockN;

        for (int j = 0; j < gridBlockN; j++) {
            equationDiffusion.conc[0][j] = iniConds.matrixConc[i][j];
            equationDiffusion.conc[1][j] = iniConds.matrixConc[i][j];
        }

        equationDiffusion.cfdProcedureOneStep(
                diffusionMath.throatConc[i],
                equationPNM.networkData.throatRadius[i],
                diffusionMath.matrixWidth[i],
                equationPNM.networkData.throatLength[i],
                diffusionMath.matricesVolume[i],
                diffusionMath.matricesOmega[i]);

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

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        flowDerivDiff[i] = ((diffFlowInstPlus[i] - diffFlowInstMinus[i]) / dP);
}

void DiffusionFlow::updateConc() {

    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < throatN; i++) {

        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            auto conc_curr = equationDiffusion.conc[equationDiffusion.iCurr][j];
            iniConds.matrixConc[i][j] = conc_curr;
        }
    }
}

void DiffusionFlow::calcDiffPart() {

    diffusionMath.calcThroatAvPress();

    // Enhance and rethink later
    dP = 0;
    diffusionMath.calcThroatConc(dP);
    calcDiffFlow(diffFlowInst);

    dP = (equationPNM.pIn - equationPNM.pOut) / 10.e+5;
    diffusionMath.calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowInstPlus);

    diffusionMath.calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowInstMinus);

    calcDiffFlowDeriv();
    updateConc();
}