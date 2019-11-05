#include <DiffusionPNM.h>
#include <numeric>

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
        conc_ini(1.),
        dP(0) {

    //    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
//        std::cout << i << std::endl;
//        for (int j = 0; j < equation.props.gridBlockN; j++) {
//            std::cout << matrixConc[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "Pressure_INI" << std::endl;
//    for (int i = 0; i < equationPNM.networkData.poreN; i++)
//        std::cout << equationPNM.pressure[i] << std::endl;

    equationPNM.setInitialCond();

    equationPNM.cfdProcedure(1,
                             equationPNM.networkData.boundaryPores,
                             equationPNM.pIn,
                             equationPNM.pOut);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << equationPNM.pressure[i] << std::endl;

//    std::cout << std::endl;
//
//    std::cout << "inletFlow" << std::endl;
//
//    for (int i = 0; i < equationPNM.inletFlow.size(); i++)
//        std::cout << equationPNM.inletFlow[i] << std::endl;

    equationPNM.cfdProcedure(0, equationPNM.networkData.boundaryPoresOut,
                             equationPNM.pIn,
                             equationPNM.pOut);

    std::cout << std::endl;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << equationPNM.pressure[i] << std::endl;

    std::cout << std::endl;

    std::cout << equationPNM.matrix << std::endl;

    std::cout << std::endl;

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
        effRadius[i] = equationPNM.networkData.throatRadius[i] * 1.65;
//        effRadius[i] = 0.0005 * rockVolume / equationPNM.networkData.throatN;
}

double DiffusionPNM::calcLangmConc(double pressure) {

    std::vector<double> &coeff = equation.props.langmuirCoeff;

    langmConc = 0;
    for (int i = 0; i < coeff.size(); i++)
        langmConc += coeff[i] * pow(pressure, i);

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

        diffFlowVector[i] = equation.flowRate / 1.653;

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

            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] / 2;

        centralCoeffDiff[i] = -1 * coeffSum;
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

            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        diffFlow[equationPNM.porConns[i][j]];
        }
        porFlowDiff[i] = -1 * coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] *
                        throatAvPress[equationPNM.porConns[i][j]];
        }
        porFlowDiffDer[i] = coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++)
        equationPNM.freeVector[i] = -1 * (porFlowDiffDer[i] + porFlowDiff[i]);

//    for (int i = 0; i < equationPNM.porConns.size(); i++)
//        std::cout << porFlowDiff[i] + porFlowDiffDer[i] << std::endl;

//    std::cout << "PorFLowDIFF" << std::endl;
//    for (int i = 0; i < porFlowDiff.size(); i++)
//        std::cout << porFlowDiff[i] << std::endl;
//
//    std::cout << "PorFLowDIFF" << std::endl;
//    for (int i = 0; i < porFlowDiffDer.size(); i++)
//        std::cout << porFlowDiffDer[i] << std::endl;


    equationPNM.calculateFreeVector(0,
                                    equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    for (int i = 0; i < equationPNM.networkData.boundaryPoresIn.size(); i++)
        equationPNM.freeVector[i] += -1 * (porFlowDiffDer[i] + porFlowDiff[i]);
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
    dP = 0;

    calcThroatConc(dP);

    calcDiffFlow(diffFlow);

//    for (int i = 0; i < diffFlow.size(); i++)
//        std::cout << diffFlow[i] << std::endl;

    dP = 100;
    calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowPlus);

    calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowMinus);

    calcDiffFlowDeriv();
//     +++++++++++++++++++++++++++++

    equationPNM.calcMatCoeff();

    calcMatCoeffDiff();

    calcMatCoupledCoeff();

    equationPNM.calculateMatrix(0,
                                equationPNM.connCoeff,
                                equationPNM.centralCoeff,
                                equationPNM.networkData.boundaryPoresOut,
                                equationPNM.porConnsIsOutByPressure,
                                connCoeffDiff);

    calcCoupledFreeVector();

    equationPNM.calculateGuessVector();

    equationPNM.calculatePress(1);

    equationPNM.getPorConnsIsOutByPressure();

    equationPNM.calcThrFlowRate();

    equationPNM.calcPorFlowRate();

    equationPNM.calcTotFlow();

    updateConc();

}

void DiffusionPNM::cfdProcedureDiff() {

    setInitialCond();

    double sum = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        sum += equationPNM.pressure[i];

    pressureAverage.emplace_back(sum / equationPNM.networkData.poreN);


    concAverage.emplace_back(conc_ini);

    for (double t = equation.props.timeStep; t <= equation.props.time;
         t += equation.props.timeStep) {

        std::vector<double> cAV(equationPNM.networkData.throatN, 0);

        calcCoupledFlow();


//
//        std::cout << std::endl;

//        for (int i = 0; i < diffFlow.size(); i++)
//            std::cout << diffFlow[i] << ' ' << diffFlowPlus[i] << ' '
//                      << diffFlowMinus[i] << std::endl;
//        std::cout << std::endl;
        double sum = 0;
        for (int i = 0; i < equationPNM.networkData.poreN; i++) {
            sum += equationPNM.pressure[i];
        }
        pressureAverage.emplace_back(sum / equationPNM.networkData.poreN);


        sum = 0;
        for (int i = 0; i < equationPNM.networkData.throatN; i++) {
            sum = 0;
            for (int j = 0; j < equation.props.gridBlockN; j++) {
                sum += matrixConc[i][j];
            }
            sum = sum / equation.props.gridBlockN;

            cAV[i] = sum;
        }

        for (int i = 0; i < cAV.size(); i++)
            std::cout << cAV[i] << ' ' << std::endl;
        std::cout << std::endl;


        double concAv = accumulate(cAV.begin(), cAV.end(), 0.0) / cAV.size();
        concAverage.emplace_back(concAv);
    }
}

const std::vector<double> DiffusionPNM::getPressureAverage() const {
    return pressureAverage;
}

const std::vector<double> DiffusionPNM::getConcAverage() const {
    return concAverage;
}


//        std::cout << "pressure" << std::endl;
//        for (int i = 0; i < equationPNM.networkData.poreN; i++)
//            std::cout << equationPNM.pressure[i] << std::endl;
//
//        std::cout << std::endl;
//

//        for (int i = 0; i < equationPNM.networkData.throatN; i++) {
//            std::cout << i << std::endl;
//            for (int j = 0; j < equation.props.gridBlockN; j++) {
//                std::cout << matrixConc[i][j] << ' ';
//            }
//            std::cout << std::endl;
//        }












