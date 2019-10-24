#ifndef PNFLOW_DIFFUSIONPNM_H
#define PNFLOW_DIFFUSIONPNM_H

#include <vector>

#include <Equation.h>
#include <EquationPNM.h>


class DiffusionPNM {

public:

    explicit DiffusionPNM(
            const std::vector<double> &propsPNM,
            const std::vector<int> &throat_list,
            const std::vector<double> &throat_radius,
            const std::vector<double> &throat_length,
            const std::vector<double> &conn_ind_in,
            const std::vector<double> &conn_ind_out,
            const std::vector<double> &pore_coord_x,
            const std::vector<double> &pore_coord_y,
            const std::vector<double> &pore_coord_z,
            const std::vector<double> &pore_radius,
            const std::vector<int> &pore_list,
            const std::vector<int> &pore_conns,
            const std::vector<int> &conn_number,
            const std::vector<int> &pore_per_row,
            const std::vector<double> &propsDiffusion,
            const std::vector<double> &langmuirCoeff);

    virtual  ~DiffusionPNM() = default;

    EquationPNM equationPNM;
    Equation equation;

    double calcSideLength(std::vector<double> &poreCoord);

    void calcRockVolume();

    void calcEffRadius();

    double calcLangmConc(double pressure);

    void calcThroatConc(const double &dP);

    void calcThroatAvPress();

    void calcDiffFlow(std::vector<double> &diffFLow);

    void calcDiffFlowDeriv();

    void calcMatCoeffDiff();

    void cfdProcedure();

    void calcMatrixJointCoeff();


    double rockVolume;
    double langmConc;
    double conc_ini;

    std::vector<double> effRadius;
    std::vector<double> throatAvPress;
    std::vector<double> throatConc;

    std::vector<double> diffFlow;
    std::vector<double> diffFlowPlus;
    std::vector<double> diffFlowMinus;

    double dP;

    std::vector<double> flowDerivDiff;

    std::vector<double> connCoeffDiff;
    std::vector<double> centralCoeffDiff;
};

#endif
