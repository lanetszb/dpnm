#include <ParamsOut.h>

#include <numeric>
#include <iomanip>

ParamsOut::ParamsOut(NetworkData &networkData, EquationPNM &equationPNM,
                     EquationDiffusion &equationDiffusion,
                     DiffusionMath &diffusionPartMath, IniConds &iniConds,
                     DiffusionFlow &diffusionPartFlow, MatrixSolver &solver,
                     const std::vector<double> &langmuirCoeff,
                     const double &matrixVolume) :

        networkData(networkData),
        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionPartMath),
        iniConds(iniConds),
        diffusionFlow(diffusionPartFlow),
        matrixSolver(solver),

        gridBlockN(equationDiffusion.propsDiffusion.gridBlockN),
        conc_ini(equationDiffusion.propsDiffusion.concIni) {}


void ParamsOut::calcMatrixMassTot() {

    double totalMass;

    std::vector<double> matrixMass;
    for (int i = 0; i < networkData.throatN; i++) {
        totalMass = 0;
        for (int j = 0; j < gridBlockN; j++) {
            totalMass += iniConds.matrixConc[i][j] *
                         equationDiffusion.localDiffusion.volCartes[j];
        }
        matrixMass.emplace_back(totalMass);
    }

    matrixMassTotal.emplace_back(
            accumulate(matrixMass.begin(), matrixMass.end(), 0.0));
}

void ParamsOut::calcTotalFlowDiff() {

    double diffusFlow;
    for (int i = 0; i < networkData.throatN; i++)
        diffusFlow += diffusionFlow.diffFlowInst[i];

    totalFlowDiff.emplace_back(diffusFlow * diffusionMath.densityConst);

}

void ParamsOut::calcPressureAv() {

    double totalPorePress;
    for (int i = 0; i < networkData.poreN; i++)
        totalPorePress += equationPNM.pressure[i];

    pressureAv.emplace_back(totalPorePress / networkData.poreN);
}

void ParamsOut::calcPressInlet() {

    // TODO: Should be already known, rather than calculated every time
    double boundPoresLeftSize = 0;
    for (int i = 0; i < networkData.poreN; i++)
        if (networkData.poreLeftX[i])
            boundPoresLeftSize += 1;

    double pressInlet = 0;
    for (int i = 0; i < networkData.poreN; i++) {
        if (networkData.poreLeftX[i])
            pressInlet += equationPNM.pressure[i];
    }

    pressIn.emplace_back(pressInlet / boundPoresLeftSize);
}

void ParamsOut::calcCoupledFlowParams() {

    // Total flow from diffusion
    calcTotalFlowDiff();
    // diffFlowThroat += diffFlowInst[i] - flowDerivDiff[i] * throatAvPress[i];

    // equationPNM.getGammaByPressure();

    equationPNM.calcThrFlowRate();
    equationPNM.calcPorFlowRate();
    equationPNM.calcTotFlow(networkData.poreLeftX);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate *
                                  diffusionMath.densityConst);

    calcPressInlet();
    calcPressureAv();
    calcMatrixMassTot();
}