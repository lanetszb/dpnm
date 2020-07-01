#ifndef PNFLOW_COUPLINGMATRIX_H
#define PNFLOW_COUPLINGMATRIX_H

#include <vector>

#include <DiffusionPartFlow.h>
#include <DiffusionPartMath.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>

class CouplingMatrix {

public:

    explicit CouplingMatrix(const std::vector<double> &propsPNM,
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

    virtual  ~CouplingMatrix() = default;

    DiffusionPartFlow diffusionPartFlow;
    DiffusionPartMath diffusionPartMath;
    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

    void calcMatCoeffDiff();

    void calcMatCoupledCoeff();

    void calcCoupledFreeVector();

    void solveCoupledMatrix();

    std::vector<double> connCoeffDiff;
    std::vector<double> centralCoeffDiff;

};

#endif //PNFLOW_COUPLINGMATRIX_H
