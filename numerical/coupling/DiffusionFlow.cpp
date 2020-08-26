#include <DiffusionFlow.h>

#include <numeric>
#include <iomanip>


DiffusionFlow::DiffusionFlow(
        const std::map<std::string, std::variant<int, double>> &paramsPnm,
        NetworkData &networkData,
        EquationPNM &equationPNM,
        EquationDiffusion &equationDiffusion,
        DiffusionMath &diffusionMath,
        IniConds &iniConds,
        const std::vector<double> &langmuirCoeff,
        const double &matrixVolume) :

        _paramsPnm(paramsPnm),
        networkData(networkData),
        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),
        iniConds(iniConds),

        gridBlockN(equationDiffusion.propsDiffusion.gridBlockN),
        diffFlowInst(networkData.fracturesN, 0),
        diffFlowInstPlus(networkData.fracturesN, 0),
        diffFlowInstMinus(networkData.fracturesN, 0),
        flowDerivDiff(networkData.fracturesN, 0),
        concCurr(networkData.fracturesN, std::vector<double>(gridBlockN, 0)),
        dP(0) {}

void DiffusionFlow::calcDiffFlow(std::vector<double> &diffFlowVector,
                                 const double &dt) {

    for (int i = 0; i < networkData.fracturesN; i++) {

        equationDiffusion.conc[0] = iniConds.matrixConc[i];
        equationDiffusion.conc[1] = iniConds.matrixConc[i];


        equationDiffusion.cfdProcedureOneStep("default",
                                              diffusionMath.throatConc[i],
                                              networkData.fracturesHeights[i],
                                              diffusionMath.matrixWidth[i],
                                              networkData.fracturesLengths[i],
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

        iniConds.matrixConc[i] = equationDiffusion.conc[equationDiffusion.iPrev];
        concCurr[i] = equationDiffusion.conc[equationDiffusion.iCurr];
    }
}

void DiffusionFlow::calcDiffFlowDeriv() {

    for (int i = 0; i < networkData.fracturesN; i++)
        flowDerivDiff[i] = ((diffFlowInstPlus[i] - diffFlowInstMinus[i]) / dP);
}

void DiffusionFlow::updateConc() {
    iniConds.matrixConc = concCurr;
}

void DiffusionFlow::calcDiffPart(const double &dt) {

    diffusionMath.calcThroatAvPress();

    // Enhance and rethink later
    dP = 0;
    diffusionMath.calcThroatConc(dP);
    calcDiffFlow(diffFlowInst, dt);

    auto &pressIn = std::get<double>(_paramsPnm["pressIn"]);
    auto &pressOut = std::get<double>(_paramsPnm["pressOut"]);

    dP = (pressIn - pressOut) / 10.e+5;
    diffusionMath.calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowInstPlus, dt);

    diffusionMath.calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowInstMinus, dt);

    calcDiffFlowDeriv();
    updateConc();
}