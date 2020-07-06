#ifndef PNFLOW_DIFFUSIONPART_H
#define PNFLOW_DIFFUSIONPART_H

#include <vector>

#include <EquationDiffusion.h>
#include <EquationPNM.h>

class DiffusionMath {

public:

    DiffusionMath(EquationPNM &equationPNM,
                  EquationDiffusion &equationDiffusion,
                  const std::vector<double> &langmuirCoeff,
                  const double &matrixVolume);

    virtual  ~DiffusionMath() = default;

    // ToDo: EquationPNM and EquationDiffusion should be references
    EquationPNM &equationPNM;
    EquationDiffusion &equationDiffusion;

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


