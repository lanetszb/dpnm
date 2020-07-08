#include <Aggregator.h>
#include <cmath>

Aggregator::Aggregator(const std::vector<double> &propsPNM,
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

        equationPNM(propsPNM, throatList, throatHeight, throatLength,
                    throatWidth, connIndIn, connIndOut, poreCoordX, poreCoordY,
                    poreCoordZ, poreRadius, poreList, poreConns, connNumber,
                    porePerRow, poreLeftX, poreRightX, hydraulicCond,
                    solverMethod),

        equationDiffusion(propsDiffusion),

        diffusionMath(equationPNM, equationDiffusion, langmuirCoeff,
                      matrixVolume),

        iniConds(equationPNM, equationDiffusion, diffusionMath,
                 langmuirCoeff, matrixVolume),

        diffusionFlow(equationPNM, equationDiffusion, diffusionMath,
                      iniConds, langmuirCoeff, matrixVolume),


        matrixSolver(equationPNM, equationDiffusion, diffusionMath,
                     diffusionFlow, langmuirCoeff, matrixVolume),

        paramsOut(equationPNM, equationDiffusion, diffusionMath,
                  iniConds, diffusionFlow, matrixSolver, langmuirCoeff,
                  matrixVolume) {}

void Aggregator::calcCoupledFlow() {

    diffusionFlow.calcDiffPart();
    matrixSolver.solveCoupledMatrix();
    paramsOut.calcCoupledFlowParams();
}

void Aggregator::cfdProcedurePnmDiff() {

    equationPNM.setInitialCondPurePnm();
    iniConds.setInitialCondCoupledMod();
    paramsOut.calcCoupledFlowParams();

    // TODO: Remove castyl

    /*auto time = equationDiffusion.propsDiffusion.time;
    auto configTimeStep = equationDiffusion.propsDiffusion.timeStep;
    double fullStepsN;
    auto lastStep = std::modf(time, &fullStepsN);
    auto timeSteps = std::vector<double>(fullStepsN, configTimeStep);
    if (lastStep > 0)
      timeSteps.push_back(lastStep);
    for(auto &&timeStep : timeSteps)
      std::cout << timeStep<<std::endl;*/

    for (double t = equationDiffusion.propsDiffusion.timeStep;
         t < equationDiffusion.propsDiffusion.time * (1. + 1.e-3);
         t += equationDiffusion.propsDiffusion.timeStep) {

        calcCoupledFlow();
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