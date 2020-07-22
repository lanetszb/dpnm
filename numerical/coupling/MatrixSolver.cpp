#include <MatrixSolver.h>

MatrixSolver::MatrixSolver(NetworkData &networkData, EquationPNM &equationPNM,
                           EquationDiffusion &equationDiffusion,
                           DiffusionMath &diffusionMath,
                           DiffusionFlow &diffusionFlow,
                           const std::vector<double> &langmuirCoeff,
                           const double &matrixVolume) :

        networkData(networkData),
        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),
        diffusionFlow(diffusionFlow),

        connCoeffDiff(networkData.fracturesN, 0),
        centralCoeffDiff(networkData.poreN, 0) {}

void MatrixSolver::calcMatCoeffDiff() {
    auto flowDerivDiff = diffusionFlow.flowDerivDiff;

    for (int i = 0; i < flowDerivDiff.size(); i++) {
        connCoeffDiff[i] = 0;
    }
    // connCoeffDiff[i] = flowDerivDiff[i] / 2;


    for (int i = 0; i < networkData.por2thrConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++)

            coeffSum += equationPNM.gammaPnm[i][j] *
                        flowDerivDiff[networkData.por2thrConns[i][j]] / 2;

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

    std::vector<double> porFlowDiff(networkData.poreN, 0);
    std::vector<double> porFlowDiffDer(networkData.poreN, 0);
    double coeffSum;

    for (int i = 0; i < networkData.por2thrConns.size(); i++) {
        coeffSum = 0;
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++) {
            coeffSum += equationPNM.gammaPnm[i][j] *
                        diffusionFlow.diffFlowInst[networkData.por2thrConns[i][j]];
        }
        porFlowDiff[i] = coeffSum;
    }

    for (int i = 0; i < networkData.por2thrConns.size(); i++) {
        coeffSum = 0;
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++) {
            coeffSum -= equationPNM.gammaPnm[i][j] *
                        diffusionFlow.flowDerivDiff[networkData.por2thrConns[i][j]] *
                        diffusionMath.throatAvPress[networkData.por2thrConns[i][j]];
        }
        porFlowDiffDer[i] = 0;
        // porFlowDiffDer[i] = coeffSum;
    }

    equationPNM.calculateFreeVector("mixed",
                                    equationPNM.pIn,
                                    equationPNM.pOut);

    for (int i = 0; i < networkData.poreN; i++)
        // TODO: understand plus or minus sign
        if (!networkData.poreOutlet[i])
            equationPNM.freeVector[i] += porFlowDiff[i];
    // equationPNM.freeVector[i] += porFlowDiff[i] + porFlowDiffDer[i];
    // equationPNM.freeVector[i] += -1 * (porFlowDiff[i] - porFlowDiffDer[i]);
}

void MatrixSolver::solveCoupledMatrix() {

    calcMatCoeffDiff();
    calcMatCoupledCoeff();

    equationPNM.calculateMatrix(equationPNM.connCoeff,
                                equationPNM.centralCoeff,
                                networkData.poreOutlet,
                                equationPNM.gammaPnm,
                                connCoeffDiff);

    calcCoupledFreeVector();
    equationPNM.calculateGuessVector();
    equationPNM.calculatePress(equationPNM.solverMethod);
}
