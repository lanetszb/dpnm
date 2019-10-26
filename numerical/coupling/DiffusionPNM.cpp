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
        matrixConc(equationPNM.networkData.throatN,
                   std::vector<double>(equation.props.gridBlockN, 0)),
        conc_ini(4),
        dP(0) {


//    for (int i = 0; i < equationPNM.porConns.size(); i++) {
//        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {
//            std::cout << equationPNM.porConns[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }

//    std::cout << std::endl;
//
//    for (int i = 0; i < equationPNM.networkData.throatN; i++)
//        std::cout << diffFlow[i] << std::endl;
//
//    std::cout << std::endl;

    cfdProcedureDiff();
//    cfdProcedureDiff();
//    cfdProcedureDiff();
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

//    equation.calcConcIni(conc_ini);
//
//    for (int i = 0; i < equationPNM.networkData.throatN; i++)
//        for (int j = 0; j < equation.props.gridBlockN; j++)
//            matrixConc[i][j] = equation.conc[equation.iPrev][j];

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {

        for (int j = 0; j < equation.props.gridBlockN; j++) {
            equation.conc[0][j] = matrixConc[i][j];
            equation.conc[1][j] = matrixConc[i][j];
        }

        equation.cfdProcedure(throatConc[i],
                              equationPNM.networkData.throatRadius[i],
                              effRadius[i],
                              equationPNM.networkData.throatLength[i]);

        for (int j = 0; j < equation.props.gridBlockN; j++) {
            matrixConc[i][j] = equation.conc[equation.iCurr][j];
        }

        diffFLow.emplace_back(equation.flowRate);

        for (int j = 0; j < equation.props.gridBlockN; j++)
            matrixConc[i][j] = equation.conc[equation.iPrev][j];
    }

//    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
//        std::cout << i << std::endl;
//        for (int j = 0; j < equation.props.gridBlockN; j++) {
//            std::cout << matrixConc[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
}

void DiffusionPNM::updateConc() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {

        for (int j = 0; j < equation.props.gridBlockN; j++) {
            equation.conc[0][j] = matrixConc[i][j];
            equation.conc[1][j] = matrixConc[i][j];
        }

        equation.cfdProcedure(throatConc[i],
                              equationPNM.networkData.throatRadius[i],
                              effRadius[i],
                              equationPNM.networkData.throatLength[i]);

        for (int j = 0; j < equation.props.gridBlockN; j++) {
            matrixConc[i][j] = equation.conc[equation.iCurr][j];
        }
    }
}

void DiffusionPNM::calcDiffFlowDeriv() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        flowDerivDiff.emplace_back(
                (diffFlowPlus[i] - diffFlowMinus[i]) / dP);
}

void DiffusionPNM::calcMatCoeffDiff() {

//    for (int i = 0; i < equationPNM.porConns.size(); i++)
//        std::cout << flowDerivDiff[i] << std::endl;

/*    std::cout << std::endl;*/

// to be rewritten

    for (int i = 0; i < flowDerivDiff.size(); i++)
        connCoeffDiff.emplace_back(flowDerivDiff[i] / 2);


    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++)
            coeffSum += flowDerivDiff[equationPNM.porConns[i][j]] / 2;
        centralCoeffDiff.emplace_back(coeffSum);
    }
}

void DiffusionPNM::calcMatCoupledCoeff() {

    equationPNM.calcMatCoeff();

    for (int i = 0; i < equationPNM.connCoeff.size(); i++)
        equationPNM.connCoeff[i] += connCoeffDiff[i];

    for (int i = 0; i < equationPNM.centralCoeff.size(); i++)
        equationPNM.centralCoeff[i] += centralCoeffDiff[i];
}

void DiffusionPNM::calcCoupledFreeVector() {

    std::vector<double> porFlowDiff;
    std::vector<double> porFlowDiffDer;

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++)
            coeffSum += diffFlow[equationPNM.porConns[i][j]];
        porFlowDiff.emplace_back(-1 * coeffSum);
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++)
            coeffSum += flowDerivDiff[equationPNM.porConns[i][j]] *
                        throatAvPress[equationPNM.porConns[i][j]];
        porFlowDiffDer.emplace_back(coeffSum);
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++)
        equationPNM.freeVector[i] = porFlowDiff[i] + porFlowDiffDer[i];

    equationPNM.calculateFreeVector(equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);
}

void DiffusionPNM::cfdProcedureDiff() {

    equationPNM.networkData.findBoundaryPores(
            equationPNM.networkData.poreCoordX);
    equationPNM.calcPorConns();
    equationPNM.calcThroatConns();

    calcRockVolume();
    calcEffRadius();

    equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    calcThroatAvPress();

    // Enhance and rethink later

    equation.calcConcIni(conc_ini);

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        for (int j = 0; j < equation.props.gridBlockN; j++)
            matrixConc[i][j] = equation.conc[equation.iCurr][j];

    calcThroatConc(dP);
    calcDiffFlow(diffFlow);

    dP = 1000;
    calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowPlus);

    calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowMinus);

    calcDiffFlowDeriv();
//     +++++++++++++++++++++++++++++

    calcMatCoeffDiff();
    calcMatCoupledCoeff();

    equationPNM.calculateMatrix();

    calcCoupledFreeVector();

    std::cout << "freeVector" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << equationPNM.freeVector[i] << std::endl;

    equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    equationPNM.calculateGuessVector();


//    std::cout << equationPNM.matrix << std::endl;
//    std::cout << std::endl;

    equationPNM.calculatePress();
//    std::cout << std::endl;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << equationPNM.pressure[i] << std::endl;

    std::cout << std::endl;


//    std::cout << equationPNM.matrix << std::endl;

    updateConc();

//    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
//        std::cout << i << std::endl;
//        for (int j = 0; j < equation.props.gridBlockN; j++) {
//            std::cout << matrixConc[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }

}













