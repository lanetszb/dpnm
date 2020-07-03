#ifndef PNFLOW_DIFFUSIONPART_H
#define PNFLOW_DIFFUSIONPART_H

#include <vector>

#include <EquationDiffusion.h>
#include <EquationPNM.h>

class DiffusionPartMath {

public:

    explicit DiffusionPartMath(const std::vector<double> &propsPNM,
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

    virtual  ~DiffusionPartMath() = default;

    // ToDo: EquationPNM and EquationDiffusion should be references
    EquationPNM equationPNM;
    EquationDiffusion equationDiffusion;

    double calcSideLength(std::vector<double> &poreCoord);

    double calcDensConst();

    void calcRockVolume();

    void calcEffRadius();

    void calcMatrixWidth();

    void calcMatricesOmega();

    void calcMatricesVolume();

    double calcLangmConc(double pressure);

    void calcThroatConc(const double &dP);

    void calcThroatAvPress();

    double matrixVolume;
    double langmConc;
    double conc_ini;
    double densityConst;

    std::vector<double> langmuirCoeff;
    std::vector<double> effRadius;
    std::vector<double> matrixWidth;

    std::vector<double> throatAvPress;
    std::vector<double> throatConc;

    std::vector<std::vector<double>> matricesOmega;
    std::vector<std::vector<double>> matricesVolume;
};

#endif //PNFLOW_DIFFUSIONPART_H


