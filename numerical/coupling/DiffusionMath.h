#ifndef PNFLOW_DIFFUSIONPART_H
#define PNFLOW_DIFFUSIONPART_H

#include <vector>

#include <EquationDiffusion.h>
#include <EquationPNM.h>
#include <NetworkData.h>

class DiffusionMath {

public:

    DiffusionMath(
            const std::map<std::string, std::variant<int, double>> &paramsPnm,
            NetworkData &networkData,
            EquationPNM &equationPNM,
            EquationDiffusion &equationDiffusion,
            const std::vector<double> &langmuirCoeff,
            const double &matrixVolume);

    virtual  ~DiffusionMath() = default;

    NetworkData &networkData;
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

    std::map<std::string, std::variant<int, double>> _paramsPnm;

    double matrixVolume;
    double langmConc;
    double conc_ini;
    double densityConst;
    int &gridBlockN;

    std::vector<double> langmuirCoeff;
    std::vector<double> effRadius;
    std::vector<double> matrixWidth;

    std::vector<double> throatAvPress;
    std::vector<double> throatConc;

    std::vector<std::vector<double>> matricesOmega;
    std::vector<std::vector<double>> matricesVolume;
};

#endif //PNFLOW_DIFFUSIONPART_H


