#ifndef PNFLOW_AGGREGATOR_H
#define PNFLOW_AGGREGATOR_H

#include <vector>

#include <CouplingParamsOut.h>
#include <CouplingMatrix.h>
#include <DiffusionPartFlow.h>
#include <CouplingIniConds.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class Aggregator {

public:

    explicit Aggregator(const std::vector<double> &propsPNM,
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
                        const double &matrixVolume);

    virtual  ~Aggregator() = default;

    CouplingParamsOut couplingParamsOut;
    CouplingMatrix couplingMatrix;
    CouplingIniConds couplingIniConds;
    DiffusionPartFlow diffusionPartFlow;
    DiffusionPartMath diffusionPartMath;
    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

    void calcCoupledFlow();

    void cfdProcedurePnmDiff();

    // Getters for Python
    const std::vector<double> getPressureAverage() const;

    const std::vector<double> getMatrixMassTotal() const;

    const std::vector<double> getTotalFlowPoresOut() const;

    const std::vector<double> getTotalFlowPoresIn() const;

    const std::vector<double> getTotalFlowDiff() const;

    const std::vector<double> getPorePressure() const;

    const std::vector<double> getInletPressure() const;
};

#endif //PNFLOW_AGGREGATOR_H
