#include <MatrixSolver.h>

#include <numeric>

MatrixSolver::MatrixSolver(EquationPNM &equationPNM,
                           EquationDiffusion &equationDiffusion,
                           DiffusionMath &diffusionMath,
                           DiffusionFlow &diffusionFlow,
                           const std::vector<double> &langmuirCoeff,
                           const double &matrixVolume) :

        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),
        diffusionFlow(diffusionFlow) {}

#include <iomanip>

void MatrixSolver::calcMatCoeffDiff() {

    auto flowDerivDiff = diffusionFlow.flowDerivDiff;

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
void MatrixSolver::calcMatCoupledCoeff() {

    equationPNM.calcMatCoeff();

    // for (int i = 0; i < equationPNM.connCoeff.size(); i++)
    //    equationPNM.connCoeff[i] += connCoeffDiff[i];

    for (int i = 0; i < equationPNM.centralCoeff.size(); i++)
        equationPNM.centralCoeff[i] += centralCoeffDiff[i];
}

void MatrixSolver::calcCoupledFreeVector() {


    std::vector<double> porFlowDiff(equationPNM.networkData.poreN, 0);
    std::vector<double> porFlowDiffDer(equationPNM.networkData.poreN, 0);

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;

        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            coeffSum += equationPNM.gammaPnm[i][j] *
                        diffusionFlow.diffFlowInst[equationPNM.porConns[i][j]];
        }
        porFlowDiff[i] = coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            coeffSum -= equationPNM.gammaPnm[i][j] *
                        diffusionFlow.flowDerivDiff[equationPNM.porConns[i][j]] *
                        diffusionMath.throatAvPress[equationPNM.porConns[i][j]];
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

void MatrixSolver::solveCoupledMatrix() {

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
