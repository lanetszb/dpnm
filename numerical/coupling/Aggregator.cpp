#include <Aggregator.h>
#include <cmath>


Aggregator::Aggregator(
        const double &matrixVolume,
        const std::string &solverMethod,
        const std::vector<double> &langmuirCoeff,
        const std::map<std::string, std::variant<int, double>> &params,
        const std::map<std::string, std::variant<int, double>> &paramsPnm,
        const std::map<std::string, std::variant<std::vector<bool>,
                std::vector<int>, std::vector<double>>> &paramsNetwork) :


        networkData(paramsNetwork),
        propsDiffusion(params),
        equationPNM(paramsPnm, networkData, solverMethod),
        equationDiffusion(propsDiffusion),
        diffusionMath(paramsPnm, networkData, equationPNM,
                      equationDiffusion, langmuirCoeff, matrixVolume),
        iniConds(paramsPnm, networkData, equationPNM, equationDiffusion,
                 diffusionMath, langmuirCoeff, matrixVolume),
        diffusionFlow(paramsPnm, networkData, equationPNM,
                      equationDiffusion, diffusionMath, iniConds,
                      langmuirCoeff, matrixVolume),
        matrixSolver(paramsPnm, networkData, equationPNM,
                     equationDiffusion, diffusionMath, diffusionFlow,
                     langmuirCoeff, matrixVolume),
        paramsOut(networkData, equationPNM, equationDiffusion,
                  diffusionMath, iniConds, diffusionFlow,
                  matrixSolver, langmuirCoeff, matrixVolume) {}

//Aggregator::Aggregator(
//        const std::map<std::string, std::variant<std::vector<bool>,
//                std::vector<int>, std::vector<double>>> &paramsNetwork,
//        const std::map<std::string, std::variant<int, double>> &params,
//        const std::vector<double> &langmuirCoeff,
//        const double &matrixVolume,
//        const std::string &solverMethod) :
//
//        Aggregator(*(new NetworkData(paramsNetwork)),
//                   *(new PropsDiffusion(params)),
//                   langmuirCoeff, matrixVolume,
//                   solverMethod) {}

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