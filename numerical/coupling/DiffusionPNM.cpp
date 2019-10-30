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
        diffFlow(equationPNM.networkData.throatN, 0),
        diffFlowPlus(equationPNM.networkData.throatN, 0),
        diffFlowMinus(equationPNM.networkData.throatN, 0),
        flowDerivDiff(equationPNM.networkData.throatN, 0),
        connCoeffDiff(equationPNM.networkData.throatN, 0),
        centralCoeffDiff(equationPNM.networkData.poreN, 0),
        conc_ini(3),
        dP(0) {

    //    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
//        std::cout << i << std::endl;
//        for (int j = 0; j < equation.props.gridBlockN; j++) {
//            std::cout << matrixConc[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
    std::cout << "Pressure_INI" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << equationPNM.pressure[i] << std::endl;

    cfdProcedureDiff();


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
        effRadius[i] = 0.00001 * rockVolume / equationPNM.networkData.throatN;
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

void DiffusionPNM::calcDiffFlow(std::vector<double> &diffFlowVector) {

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

        diffFlowVector[i] = equation.flowRate / 1000;

        for (int j = 0; j < equation.props.gridBlockN; j++)
            matrixConc[i][j] = equation.conc[equation.iPrev][j];
    }
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
        flowDerivDiff[i] = ((diffFlowPlus[i] - diffFlowMinus[i]) / dP);
}

void DiffusionPNM::calcMatCoeffDiff() {

    for (int i = 0; i < flowDerivDiff.size(); i++)
        connCoeffDiff[i] = flowDerivDiff[i] / 2;


    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++)
//            if (equationPNM.porConnsIsOut[i][j])
            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] / 2;
//            else
//        coeffSum -= flowDerivDiff[equationPNM.porConns[i][j]] / 2;
        centralCoeffDiff[i] = coeffSum;
    }
}

// potential error
void DiffusionPNM::calcMatCoupledCoeff() {

//    equationPNM.calcMatCoeff();

//    for (int i = 0; i < equationPNM.connCoeff.size(); i++)
//        equationPNM.connCoeff[i] += connCoeffDiff[i];

    for (int i = 0; i < equationPNM.centralCoeff.size(); i++)
        equationPNM.centralCoeff[i] += centralCoeffDiff[i];
}

void DiffusionPNM::calcCoupledFreeVector() {

    std::vector<double> porFlowDiff(equationPNM.networkData.poreN, 0);
    std::vector<double> porFlowDiffDer(equationPNM.networkData.poreN, 0);

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {
//            if (equationPNM.porConnsIsOut[i][j])
            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        diffFlow[equationPNM.porConns[i][j]];
//            else
//                coeffSum -= diffFlow[equationPNM.porConns[i][j]];
        }
        porFlowDiff[i] = -1 * coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {
//            if (equationPNM.porConnsIsOut[i][j])
            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] *
                        throatAvPress[equationPNM.porConns[i][j]];
//            else
//                coeffSum -= flowDerivDiff[equationPNM.porConns[i][j]] *
//                            throatAvPress[equationPNM.porConns[i][j]];
        }
        porFlowDiffDer[i] = coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++)
        equationPNM.freeVector[i] = porFlowDiff[i] + porFlowDiffDer[i];


    equationPNM.calculateFreeVector(equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);
}

void DiffusionPNM::setInitialCond() {

//    equationPNM.networkData.findBoundaryPores(
//            equationPNM.networkData.poreCoordX);
//    equationPNM.calcPorConns();
//    equationPNM.calcThroatConns();

    calcRockVolume();
    calcEffRadius();

//    equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
//                                    equationPNM.propsPNM.pressOut);

    equation.calcConcIni(conc_ini);

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        for (int j = 0; j < equation.props.gridBlockN; j++)
            matrixConc[i][j] = equation.conc[equation.iCurr][j];

//    equationPNM.getPorConnsIsOut();
}

void DiffusionPNM::calcCoupledFlow() {

    calcThroatAvPress();

    // Enhance and rethink later

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

    equationPNM.calculateMatrix(connCoeffDiff);

    calcCoupledFreeVector();

    equationPNM.calculateGuessVector();

    equationPNM.calculatePress();

    equationPNM.calcPorFlowRate();

    updateConc();

}

void DiffusionPNM::cfdProcedureDiff() {

    setInitialCond();

//    std::cout << "OUT" << std::endl;
//    for (int i = 0; i < equationPNM.porConns.size(); i++) {
//        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {
//            std::cout << equationPNM.porConnsIsOut[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }


    for (double t = equation.props.timeStep;
         t <= 1000 * equation.props.timeStep; t += equation.props.timeStep) {

        calcCoupledFlow();

//        double sum = 0;
//        for (int i = 0; i < equationPNM.porFlowRate.size();i++)
//            sum+=equationPNM.porFlowRate[i];
//
//        std::cout << "totalFlow" << std::endl;
//        std::cout << sum << std::endl;
//        std::cout << std::endl;
//
//        std::cout << "freeVector" << std::endl;
//        std::cout << equationPNM.freeVector << std::endl;
//        std::cout << std::endl;
//
        for (int i = 0; i < equationPNM.networkData.poreN; i++)
            std::cout << equationPNM.pressure[i] << std::endl;
        std::cout << std::endl;
//        std::cout << std::endl;
//
        for (int i = 0; i < equationPNM.networkData.throatN; i++) {
            std::cout << i << std::endl;
            for (int j = 0; j < equation.props.gridBlockN; j++) {
                std::cout << matrixConc[i][j] << ' ';
            }
            std::cout << std::endl;
        }
    }


//
//    for (int i = 0; i < equationPNM.networkData.throatN; i++)
//        std::cout << diffFlow[i] << std::endl;

//    std::cout << std::endl;
//
}















