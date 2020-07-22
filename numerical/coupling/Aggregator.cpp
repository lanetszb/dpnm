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
                       const std::vector<int> &fracturesList,
                       const std::vector<double> &fracturesHeights,
                       const std::vector<double> &fracturesLengths,
                       const std::vector<double> &fracturesWidths,
                       const std::vector<double> &fracsConnIndIn,
                       const std::vector<double> &fracsConnIndOut,
                       const std::vector<double> &poresCoordsX,
                       const std::vector<double> &poresCoordsY,
                       const std::vector<double> &poresCoordsZ,
                       const std::vector<double> &poresRadii,
                       const std::vector<int> &poresList,
                       const std::vector<bool> &poresInlet,
                       const std::vector<bool> &poresOutlet,
                       const std::vector<double> &hydraulicCond,
                       const std::vector<double> &langmuirCoeffs,
                       const double &matrixVolume,
                       const std::string &solverMethod) :
        Aggregator(*(new PropsPnm(propsVector)),
                   *(new NetworkData(fracturesList, fracturesHeights,
                                     fracturesLengths, fracturesWidths,
                                     fracsConnIndIn, fracsConnIndOut,
                                     poresCoordsX, poresCoordsY, poresCoordsZ,
                                     poresRadii, poresList, poresInlet,
                                     poresOutlet, hydraulicCond)),
                   propsDiffusion, langmuirCoeffs, matrixVolume, solverMethod) {}

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