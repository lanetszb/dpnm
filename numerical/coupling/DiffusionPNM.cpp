#include <DiffusionPNM.h>

DiffusionPNM::DiffusionPNM(
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
        const std::vector<double> &langmuirCoeff) :


        equationPNM(propsPNM, throat_list, throat_radius, throat_length,
                    conn_ind_in, conn_ind_out, pore_coord_x, pore_coord_y,
                    pore_coord_z, pore_radius, pore_list, pore_conns,
                    conn_number, pore_per_row),
        equation(propsDiffusion, langmuirCoeff),
        effRadius(equationPNM.networkData.throatN, 0),
        throatAvPress(equationPNM.networkData.throatN, 0),
        throatConc(equationPNM.networkData.throatN, 0),
        conc_ini(121) {

//    calcRockVolume();
//    calcEffRadius();

//    equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
//                                    equationPNM.propsPNM.pressOut);
//
//    equationPNM.calcThroatConns();
//    calcThroatAvPress();
//
//    for (int i = 0; i < equationPNM.networkData.throatN; i++)
//        std::cout << throatAvPress[i] << std::endl;
//    std::cout << std::endl;
//
//
//    calcThroatConc();
//
//    for (int i = 0; i < equationPNM.networkData.throatN; i++)
//        std::cout << throatConc[i] << std::endl;
//    std::cout << std::endl;

}

double DiffusionPNM::calcSideLength(std::vector<double> &poreCoord) {

    auto min = std::min_element(std::begin(poreCoord),
                                std::end(poreCoord));
    auto max = std::max_element(std::begin(poreCoord),
                                std::end(poreCoord));

    return *max - *min;
}


void DiffusionPNM::calcRockVolume() {

    auto lengthX = calcSideLength(equationPNM.networkData.poreCoordX);
    auto lengthY = calcSideLength(equationPNM.networkData.poreCoordY);
    auto lengthZ = calcSideLength(equationPNM.networkData.poreCoordZ);

    rockVolume = lengthX * lengthY * lengthZ;
}

void DiffusionPNM::calcEffRadius() {

    for (int i = 0; i < effRadius.size(); i++)
        effRadius[i] = rockVolume / equationPNM.networkData.throatN;
}

double DiffusionPNM::calcLangmConc(double pressure) {


    langmConc =
            equation.props.langmuirCoeff[0] * pressure * pressure * pressure *
            pressure +
            equation.props.langmuirCoeff[1] * pressure * pressure * pressure +
            equation.props.langmuirCoeff[2] * pressure * pressure +
            equation.props.langmuirCoeff[3] * pressure +
            equation.props.langmuirCoeff[4];

    return langmConc;
}

void DiffusionPNM::calcThroatAvPress() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        throatAvPress[i] =
                (equationPNM.pressure[equationPNM.throatConns[i].first]
                 + equationPNM.pressure[equationPNM.throatConns[i].second]) / 2;

}

void DiffusionPNM::calcThroatConc() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        throatConc[i] = calcLangmConc(throatAvPress[i]);
}

void DiffusionPNM::calcDiffFlow() {

    calcRockVolume();
    calcEffRadius();

    equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    equationPNM.calcThroatConns();
    calcThroatAvPress();
    calcThroatConc();



}

void DiffusionPNM::calcDiffFlowDeriv() {

    double dP = 100;
}













