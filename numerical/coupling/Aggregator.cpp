#include <Aggregator.h>
#include <cmath>


Aggregator::Aggregator(PropsPnm &propsPnm, NetworkData &networkData,
                       const std::vector<double> &propsDiffusion,
                       const std::vector<double> &langmuirCoeff,
                       const double &matrixVolume,
                       const std::string &solverMethod) :

        propsPnm(propsPnm),
        networkData(networkData),

        equationPNM(propsPnm, networkData, solverMethod),
        equationDiffusion(propsDiffusion),
        diffusionMath(propsPnm, networkData, equationPNM,
                      equationDiffusion, langmuirCoeff, matrixVolume),
        iniConds(networkData, equationPNM, equationDiffusion,
                 diffusionMath, langmuirCoeff, matrixVolume),
        diffusionFlow(networkData, equationPNM, equationDiffusion,
                      diffusionMath, iniConds, langmuirCoeff, matrixVolume),
        matrixSolver(networkData, equationPNM, equationDiffusion,
                     diffusionMath,
                     diffusionFlow, langmuirCoeff, matrixVolume),
        paramsOut(networkData, equationPNM, equationDiffusion,
                  diffusionMath,
                  iniConds, diffusionFlow, matrixSolver, langmuirCoeff,
                  matrixVolume) {}

Aggregator::Aggregator(const std::vector<double> &propsVector,
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
                       const double &matrixVolume,
                       const std::string &solverMethod) :
        Aggregator(*(new PropsPnm(propsVector)),
                   *(new NetworkData(throatList, throatHeight, throatLength,
                                     throatWidth, connIndIn, connIndOut,
                                     poreCoordX,
                                     poreCoordY,
                                     poreCoordZ, poreRadius, poreList,
                                     poreConns,
                                     connNumber,
                                     porePerRow, poreLeftX, poreRightX,
                                     hydraulicCond)), propsDiffusion,
                   langmuirCoeff, matrixVolume, solverMethod) {}

void Aggregator::calcCoupledFlow(const double &dt) {

    diffusionFlow.calcDiffPart(dt);
    matrixSolver.solveCoupledMatrix();
    paramsOut.calcCoupledFlowParams();
}

void Aggregator::cfdProcedurePnmDiff() {

    equationPNM.setInitialCondPurePnm();
    iniConds.setInitialCondCoupledMod();
    paramsOut.calcCoupledFlowParams();
    equationDiffusion.calcTimeVector();

    for (int i = 0; i <= equationDiffusion.timeStepsVec.size() - 1; i++) {
        calcCoupledFlow(equationDiffusion.timeStepsVec[i]);
    }
}

// Getters for Python
const std::vector<double> Aggregator::getPressureAverage() const {
    return paramsOut.pressureAv;
}

const std::vector<double> Aggregator::getMatrixMassTotal() const {
    return paramsOut.matrixMassTotal;
}

const std::vector<double> Aggregator::getTotalFlowPoresOut() const {
    return paramsOut.totalFlowPoresOut;
}

const std::vector<double> Aggregator::getTotalFlowPoresIn() const {
    return paramsOut.totalFlowPoresIn;
}

const std::vector<double> Aggregator::getTotalFlowDiff() const {
    return paramsOut.totalFlowDiff;
}

const std::vector<double> Aggregator::getInletPressure() const {
    return paramsOut.pressIn;
}

const std::vector<double> Aggregator::getPorePressure() const {
    return equationPNM.pressure;
}

const std::vector<double> Aggregator::getTimeStepsVec() const {
    return equationDiffusion.timeStepsVec;
}