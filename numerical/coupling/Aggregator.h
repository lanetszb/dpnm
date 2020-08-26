#ifndef PNFLOW_AGGREGATOR_H
#define PNFLOW_AGGREGATOR_H

#include <vector>

#include <NetworkData.h>
#include <EquationPNM.h>
#include <PropsDiffusion.h>
#include <EquationDiffusion.h>
#include <IniConds.h>
#include <DiffusionFlow.h>
#include <MatrixSolver.h>
#include <ParamsOut.h>


class Aggregator {

public:

    explicit Aggregator(const double &matrixVolume,
                        const std::string &solverMethod,
                        const std::vector<double> &langmuirCoeff,
                        const std::map<std::string, std::variant<int, double>> &params,
                        const std::map<std::string, std::variant<int, double>> &paramsPnm,
                        const std::map<std::string, std::variant<std::vector<bool>,
                                std::vector<int>, std::vector<double>>> &paramsNetwork);


//    explicit Aggregator(
//            const std::map<std::string, std::variant<std::vector<bool>,
//                    std::vector<int>, std::vector<double>>> &paramsNetwork,
//            const std::map<std::string, std::variant<int, double>> &params,
//            const std::vector<double> &langmuirCoeffs,
//            const double &matrixVolume,
//            const std::string &solverMethod);


    virtual  ~Aggregator() = default;

    NetworkData networkData;
    PropsDiffusion propsDiffusion;

    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

    DiffusionMath diffusionMath;
    IniConds iniConds;
    DiffusionFlow diffusionFlow;
    MatrixSolver matrixSolver;
    ParamsOut paramsOut;

    void calcCoupledFlow(const double &dt);

    void cfdProcedurePnmDiff();

    // Getters for Python
    const std::vector<double> getPressureAverage() const;

    const std::vector<double> getMatrixMassTotal() const;

    const std::vector<double> getTotalFlowPoresOut() const;

    const std::vector<double> getTotalFlowPoresIn() const;

    const std::vector<double> getTotalFlowDiff() const;

    const std::vector<double> getPorePressure() const;

    const std::vector<double> getInletPressure() const;

    const std::vector<double> getTimeStepsVec() const;
};

#endif //PNFLOW_AGGREGATOR_H
