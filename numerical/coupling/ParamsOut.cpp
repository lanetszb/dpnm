#include <ParamsOut.h>

#include <numeric>
#include <iomanip>

ParamsOut::ParamsOut(EquationPNM &equationPNM,
                     EquationDiffusion &equationDiffusion,
                     DiffusionMath &diffusionPartMath,
                     IniConds &iniConds,
                     DiffusionFlow &diffusionPartFlow,
                     MatrixSolver &solver,
                     const std::vector<double> &langmuirCoeff,
                     const double &matrixVolume) :

        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionPartMath),
        iniConds(iniConds),
        diffusionFlow(diffusionPartFlow),
        matrixSolver(solver),
        conc_ini(equationDiffusion.propsDiffusion.concIni) {}


void ParamsOut::calcMatrixMassTot() {

    double sum = 0;
    std::vector<double> matrixMass;
    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
        sum = 0;
        for (int j = 0;
             j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            sum += iniConds.matrixConc[i][j] *
                   equationDiffusion.localDiffusion.volCartes[j];
        }
        matrixMass.emplace_back(sum);
    }

    matrixMassTotal.emplace_back(
            accumulate(matrixMass.begin(), matrixMass.end(), 0.0));
}

void ParamsOut::calcVecSum(const int &iter,
                           const std::vector<double> &vectorToSum,
                           std::vector<double> &vectorSum,
                           const double &mult) {

    double sum = 0;
    for (int i = 0; i < iter; i++)
        sum += vectorToSum[i];

    vectorSum.emplace_back(sum * mult);
}

void ParamsOut::calcPressInlet() {

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

void ParamsOut::calcCoupledFlowParams() {

    // Total flow from diffusion
    calcVecSum(equationPNM.networkData.throatN, diffusionFlow.diffFlowInst,
               totalFlowDiff, diffusionMath.densityConst);
    // diffFlowThroat += diffFlowInst[i] - flowDerivDiff[i] * throatAvPress[i];

    // equationPNM.getGammaByPressure();

    equationPNM.calcThrFlowRate();
    equationPNM.calcPorFlowRate();
    equationPNM.calcTotFlow(equationPNM.networkData.poreLeftX);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate *
                                  diffusionMath.densityConst);
}