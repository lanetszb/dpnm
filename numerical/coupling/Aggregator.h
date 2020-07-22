#ifndef PNFLOW_AGGREGATOR_H
#define PNFLOW_AGGREGATOR_H

#include <vector>

#include <NetworkData.h>
#include <PropsPnm.h>
#include <EquationPNM.h>
#include <EquationDiffusion.h>
#include <IniConds.h>
#include <DiffusionFlow.h>
#include <MatrixSolver.h>
#include <ParamsOut.h>


class Aggregator {

public:

    explicit Aggregator(PropsPnm &propsPnm, NetworkData &networkData,
                        const std::vector<double> &propsDiffusion,
                        const std::vector<double> &langmuirCoeff,
                        const double &matrixVolume,
                        const std::string &solverMethod);

    explicit Aggregator(const std::vector<double> &propsVector,
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
                        const std::string &solverMethod);


    virtual  ~Aggregator() = default;

    PropsPnm &propsPnm;
    NetworkData &networkData;

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
