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
        conc_ini(4),
        dP(0) {

    cfdProcedure();
    calcMatCoeffDiff();
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

void DiffusionPNM::calcThroatConc(const double &dP) {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        throatConc[i] = calcLangmConc(throatAvPress[i] + dP);
}

void DiffusionPNM::calcDiffFlow(std::vector<double> &diffFLow) {

    equation.calcConcIni(conc_ini);
    std::vector<std::vector<double>> copyConcIni;

    copyConcIni.push_back(equation.conc[0]);
    copyConcIni.push_back(equation.conc[1]);


    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
        equation.cfdProcedure(throatConc[i],
                              equationPNM.networkData.throatRadius[i],
                              effRadius[i],
                              equationPNM.networkData.throatLength[i]);

        diffFLow.emplace_back(equation.flowRate);

        for (int i = 0; i < equation.props.gridBlockN; i++) {
            equation.conc[0][i] = copyConcIni[0][i];
            equation.conc[1][i] = copyConcIni[1][i];
        }
    }
}

void DiffusionPNM::calcDiffFlowDeriv() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        flowDerivDiff.emplace_back((diffFlowPlus[i] - diffFlowMinus[i]) / dP);
}

void DiffusionPNM::cfdProcedure() {

    calcRockVolume();
    calcEffRadius();

    equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    equationPNM.calcThroatConns();
    calcThroatAvPress();

    // Enhance and rethink later

    calcThroatConc(dP);
    calcDiffFlow(diffFlow);

    dP = 100;
    calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowPlus);

    calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowMinus);

    calcDiffFlowDeriv();
    // +++++++++++++++++++++++++++++
}

void DiffusionPNM::calcMatCoeffDiff() {

    equationPNM.calcPorConns();

    for (int i = 0; i < flowDerivDiff.size(); i++)
        connCoeffDiff.emplace_back(flowDerivDiff[i] / 2);


    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++)
            coeffSum += flowDerivDiff[equationPNM.porConns[i][j]] / 2;
        centralCoeffDiff.emplace_back(coeffSum);
    }

//    std::cout << std::endl;
//
//    for (int i = 0; i < equationPNM.porConns.size(); i++)
//        std::cout << centralCoeffDiff[i] << std::endl;
}














