#include <DiffusionPNM.h>
#include <numeric>

#include <iomanip>


DiffusionPNM::DiffusionPNM(const std::vector<double> &propsPNM,
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
                           const std::vector<double> &langmuirCoeff) :

        equationPNM(propsPNM, throatList, throatHeight, throatLength,
                    throatWidth, connIndIn, connIndOut, poreCoordX, poreCoordY,
                    poreCoordZ, poreRadius, poreList, poreConns, connNumber,
                    porePerRow, poreLeftX, poreRightX, hydraulicCond),

        equationDiffusion(propsDiffusion),
        langmuirCoeff(langmuirCoeff),
        effRadius(equationPNM.networkData.throatN, 0),
        matrixWidth(equationPNM.networkData.throatN, 0),
        throatAvPress(equationPNM.networkData.throatN, 0),
        throatConc(equationPNM.networkData.throatN, 0),
        matrixConc(equationPNM.networkData.throatN,
                   std::vector<double>(
                           equationDiffusion.propsDiffusion.gridBlockN, 0)),
        diffFlow(equationPNM.networkData.throatN, 0),
        diffFlowPlus(equationPNM.networkData.throatN, 0),
        diffFlowMinus(equationPNM.networkData.throatN, 0),
        flowDerivDiff(equationPNM.networkData.throatN, 0),
        connCoeffDiff(equationPNM.networkData.throatN, 0),
        centralCoeffDiff(equationPNM.networkData.poreN, 0),
        conc_ini(3.0),
        dP(0) {

//
//    //    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
//        std::cout << i << std::endl;
//        for (int j = 0; j < equation.propsDiffusion.gridBlockN; j++) {
//            std::cout << matrixConc[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "Pressure_INI" << std::endl;
//    for (int i = 0; i < equationPNM.networkData.poreN; i++)
//        std::cout << equationPNM.pressure[i] << std::endl;
//
    equationPNM.setInitialCond();

    std::vector<int> boundPoresInput;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i] or
            equationPNM.networkData.poreRightX[i])
            boundPoresInput.emplace_back(i);

    equationPNM.cfdProcedure(1,
                             boundPoresInput,
                             equationPNM.pIn,
                             equationPNM.pOut);

    std::cout << "Pressure Dirichlet" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(9) << equationPNM.pressure[i];
    std::cout << std::endl;


    // std::cout << equationPNM.matrix << std::endl;
    //
    // std::cout << std::endl;
    //
    // std::cout << "centralCoeff" << std::endl;
    // for (int i = 0; i < equationPNM.centralCoeff.size(); i++)
    //     std::cout << equationPNM.centralCoeff[i] << std::endl;
    // std::cout << std::endl;

    std::cout << "inletFlow" << std::endl;

    for (int i = 0; i < equationPNM.inletFlow.size(); i++)
        std::cout << equationPNM.inletFlow[i] << std::endl;

    std::cout << "inletFlow" << std::endl;

    std::vector<int> boundPoresRight;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreRightX[i])
            boundPoresRight.emplace_back(i);

    equationPNM.cfdProcedure(0,
                             boundPoresRight,
                             equationPNM.pIn,
                             equationPNM.pOut);

    std::cout << "Pressure Dirichlet+Neumann" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(9) << equationPNM.pressure[i];
    std::cout << std::endl;

    // std::cout << std::endl;
    //
    // std::cout << equationPNM.matrix << std::endl;
    //
    // std::cout << std::endl;
    //
    // for (int i = 0; i < equationPNM.networkData.poreN; i++)
    //     std::cout << equationPNM.pressure[i] << std::endl;
    //
    // std::cout << std::endl;
    //
    // std::cout << "centralCoeff" << std::endl;
    // for (int i = 0; i < equationPNM.centralCoeff.size(); i++)
    //     std::cout << equationPNM.centralCoeff[i] << std::endl;
    // std::cout << std::endl;


    cfdProcedureDiff();

    std::cout << "Pressure Final" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(9) << equationPNM.pressure[i];
    std::cout << std::endl;
}

//
double DiffusionPNM::calcSideLength(std::vector<double> &poreCoord) {

    auto min = std::min_element(std::begin(poreCoord),
                                std::end(poreCoord));
    auto max = std::max_element(std::begin(poreCoord),
                                std::end(poreCoord));

    return *max - *min;
}

double DiffusionPNM::calcDensConst() {

    auto aGasDens = equationPNM.propsPNM.aGasDens;
    auto bGasDens = equationPNM.propsPNM.bGasDens;

    auto pressIn = equationPNM.propsPNM.pressIn;
    auto pressOut = equationPNM.propsPNM.pressOut;

    auto pressureAv = (pressIn + pressOut) / 2;

    return aGasDens * pressureAv + bGasDens;
}

void DiffusionPNM::calcRockVolume() {

    auto lengthX = calcSideLength(equationPNM.networkData.poreCoordX);
    auto lengthY = calcSideLength(equationPNM.networkData.poreCoordY);
    auto lengthZ = calcSideLength(equationPNM.networkData.poreCoordZ);

    // TODO: to make the option of imporing real rock volume like is done below
    // rockVolume = (lengthX * lengthY * lengthZ);
    rockVolume = 0.018 * 0.018 * 0.018 * 0.0001;

    // std::cout << "rockVolume= " << rockVolume << std::endl;
}

void DiffusionPNM::calcEffRadius() {

    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < effRadius.size(); i++)
        effRadius[i] = rockVolume / throatN;
}

void DiffusionPNM::calcMatrixWidth() {

    auto throatN = equationPNM.networkData.throatN;

    calcEffRadius();


    for (int i = 0; i < throatN; i++) {

        auto fracHeight = equationPNM.networkData.throatRadius[i];
        auto fracLength = equationPNM.networkData.throatLength[i];
        auto fracWidth = equationPNM.networkData.throatWidth[i];

        matrixWidth[i] = effRadius[i] / fracLength / fracHeight + fracWidth;
    }
}

double DiffusionPNM::calcLangmConc(double pressure) {

    langmConc = 0;
    for (int i = 0; i < langmuirCoeff.size(); i++)
        langmConc += langmuirCoeff[i] * pow(pressure, i);

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

        auto gridBlockN = equationDiffusion.propsDiffusion.gridBlockN;

        auto fracHeight = equationPNM.networkData.throatRadius[i];
        auto fracLength = equationPNM.networkData.throatLength[i];
        auto fracWidth = equationPNM.networkData.throatWidth[i];

        equationDiffusion.convectiveDiffusion.calcOmegaCartes(fracHeight,
                                                              fracLength);

        equationDiffusion.localDiffusion.calcVolCartesian(fracHeight,
                                                          matrixWidth[i],
                                                          fracLength,
                                                          fracWidth);

        for (int j = 0; j < gridBlockN; j++) {
            equationDiffusion.conc[0][j] = matrixConc[i][j];
            equationDiffusion.conc[1][j] = matrixConc[i][j];
        }

        equationDiffusion.cfdProcedureOneStep(
                throatConc[i],
                equationPNM.networkData.throatRadius[i],
                matrixWidth[i],
                equationPNM.networkData.throatLength[i],
                equationDiffusion.localDiffusion.volCartes,
                equationDiffusion.convectiveDiffusion.omegaCartesian);

        for (int j = 0; j < gridBlockN; j++) {
            matrixConc[i][j] = equationDiffusion.conc[equationDiffusion.iCurr][j];
        }

        diffFlowVector[i] = equationDiffusion.flowRate / densityConst;

        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++)
            matrixConc[i][j] = equationDiffusion.conc[equationDiffusion.iPrev][j];
    }

    // for (int i = 0; i < equationPNM.networkData.throatN; i++) {
    //     // std::cout << i << std::endl;
    //     for (int j = 0;
    //          j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
    //         std::cout << matrixConc[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
}

void DiffusionPNM::updateConc() {

    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < throatN; i++) {

        auto gridBlockN = equationDiffusion.propsDiffusion.gridBlockN;

        for (int j = 0; j < gridBlockN; j++) {
            equationDiffusion.conc[0][j] = matrixConc[i][j];
            equationDiffusion.conc[1][j] = matrixConc[i][j];
        }


        auto fracHeight = equationPNM.networkData.throatRadius[i];
        auto fracLength = equationPNM.networkData.throatLength[i];
        auto fracWidth = equationPNM.networkData.throatWidth[i];

        equationDiffusion.convectiveDiffusion.calcOmegaCartes(fracHeight,
                                                              fracLength);

        equationDiffusion.localDiffusion.calcVolCartesian(fracHeight,
                                                          matrixWidth[i],
                                                          fracLength,
                                                          fracWidth);

        equationDiffusion.cfdProcedureOneStep(
                throatConc[i],
                equationPNM.networkData.throatRadius[i],
                matrixWidth[i],
                equationPNM.networkData.throatLength[i],
                equationDiffusion.localDiffusion.volCartes,
                equationDiffusion.convectiveDiffusion.omegaCartesian);

        // equationDiffusion.cfdProcedure(
        //         boundCond,
        //         throatConc[i],
        //         equationPNM.networkData.throatRadius[i],
        //         matrixWidth[i],
        //         equationPNM.networkData.throatLength[i],
        //         equationDiffusion.localDiffusion.volCartes,
        //         equationDiffusion.convectiveDiffusion.omegaCartesian);

        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            matrixConc[i][j] = equationDiffusion.conc[equationDiffusion.iCurr][j];
        }
    }

    // for (int i = 0; i < equationPNM.networkData.throatN; i++) {
    //     // std::cout << i << std::endl;
    //     for (int j = 0;
    //          j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
    //         std::cout << matrixConc[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
}

//
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

        // TODO: understand plus or minus sign and properly name as derivDiff
        centralCoeffDiff[i] = -1 * coeffSum;
        // centralCoeffDiff[i] = coeffSum;
    }
}

//
// potential error
void DiffusionPNM::calcMatCoupledCoeff() {

    equationPNM.calcMatCoeff();

    // for (int i = 0; i < equationPNM.connCoeff.size(); i++)
    //    equationPNM.connCoeff[i] += connCoeffDiff[i];

    for (int i = 0; i < equationPNM.centralCoeff.size(); i++)
        equationPNM.centralCoeff[i] += centralCoeffDiff[i];

}

void DiffusionPNM::calcCoupledFreeVector() {

    std::vector<double> porFlowDiff(equationPNM.networkData.poreN, 0);
    std::vector<double> porFlowDiffDer(equationPNM.networkData.poreN, 0);

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            // TODO: understand plus or minus sign and properly name as derivDiff
            // coeffSum -= equationPNM.porConnsIsOutByPressure[i][j] *
            //            diffFlow[equationPNM.porConns[i][j]];
            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        diffFlow[equationPNM.porConns[i][j]];
        }
        porFlowDiff[i] = coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            // TODO: understand plus or minus sign and properly name as derivDiff
            // coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
            //            flowDerivDiff[equationPNM.porConns[i][j]] *
            //            throatAvPress[equationPNM.porConns[i][j]];

            coeffSum -= equationPNM.porConnsIsOutByPressure[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] *
                        throatAvPress[equationPNM.porConns[i][j]];
        }
        porFlowDiffDer[i] = coeffSum;
    }

//    for (int i = 0; i < equationPNM.porConns.size(); i++)
//        std::cout << porFlowDiff[i] + porFlowDiffDer[i] << std::endl;

    // std::cout << "PorFLowDIFF" << std::endl;
    // for (int i = 0; i < porFlowDiff.size(); i++)
    //     std::cout << porFlowDiff[i] << std::endl;
    //
    // std::cout << std::endl;
    //
    // std::cout << std::endl;
    // std::cout << "PorFLowDIFF_Der" << std::endl;
    // for (int i = 0; i < porFlowDiffDer.size(); i++)
    //     std::cout << porFlowDiffDer[i] << std::endl;
    // std::cout << std::endl;

    equationPNM.calculateFreeVector(0,
                                    equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    for (int i = 0; i < equationPNM.networkData.boundaryPoresIn.size(); i++)
        // TODO: understand plus or minus sign
        equationPNM.freeVector[i] += porFlowDiff[i] + porFlowDiffDer[i];
    // equationPNM.freeVector[i] += -1 * (porFlowDiff[i] - porFlowDiffDer[i]);
}

//
void DiffusionPNM::setInitialCond() {

    auto throatN = equationPNM.networkData.throatN;

    densityConst = calcDensConst();

    equationPNM.calcPorConns();
    equationPNM.calcThroatConns();

    calcRockVolume();
    calcMatrixWidth();

    equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    equationDiffusion.calcConcIni(conc_ini);

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++)
            matrixConc[i][j] = conc_ini;

    equationPNM.getPorConnsIsOut();
}

void DiffusionPNM::calcCoupledFlow() {

    calcThroatAvPress();

    // Enhance and rethink later
    dP = 0;

    calcThroatConc(dP);

    calcDiffFlow(diffFlow);

    for (int i = 0; i < diffFlow.size(); i++)
        std::cout << diffFlow[i] << std::endl;
    std::cout << std::endl;

    dP = 0.0001;
    calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowPlus);

    calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowMinus);

    calcDiffFlowDeriv();

    calcMatCoeffDiff();

    calcMatCoupledCoeff();

    std::vector<int> boundPoresRight;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreRightX[i])
            boundPoresRight.emplace_back(i);

    std::vector<int> boundPoresLeft;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i])
            boundPoresLeft.emplace_back(i);

    equationPNM.calculateMatrix(0,
                                equationPNM.connCoeff,
                                equationPNM.centralCoeff,
                                boundPoresRight,
                                equationPNM.porConnsIsOutByPressure,
                                connCoeffDiff);

    calcCoupledFreeVector();

    equationPNM.calculateGuessVector();

    equationPNM.calculatePress(1);

    // equationPNM.getPorConnsIsOutByPressure();

    equationPNM.calcThrFlowRate();

    equationPNM.calcPorFlowRate();

    equationPNM.calcTotFlow(boundPoresRight);
    totalFlowPoresOut.emplace_back(-1 * equationPNM.totFlowRate);

    equationPNM.calcTotFlow(boundPoresLeft);

    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate);

    updateConc();

}

//
void DiffusionPNM::cfdProcedureDiff() {

    setInitialCond();

    double sum = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        sum += equationPNM.pressure[i];

    pressureAverage.emplace_back(sum / equationPNM.networkData.poreN);

    concAverage.emplace_back(conc_ini);
    totalFlowDiff.emplace_back(0);

    totalFlowPoresOut.emplace_back(equationPNM.totFlowRate);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate);

    // TODO: inlet pressure to be rewritten later
    std::vector<int> boundPoresLeft;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i])
            boundPoresLeft.emplace_back(i);

    double pressInlet = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i])
            pressInlet += equationPNM.pressure[i];

    pressIn.emplace_back(pressInlet / boundPoresLeft.size());

    // To revrite latter


    for (double t = equationDiffusion.propsDiffusion.timeStep;
         t <= equationDiffusion.propsDiffusion.time;
         t += equationDiffusion.propsDiffusion.timeStep) {

        //std::cout << "freeVector" << std::endl;
        //std::cout << equationPNM.freeVector << std::endl;;
        //std::cout << std::endl;

        std::vector<double> cAV(equationPNM.networkData.throatN, 0);

        // for (int i = 0; i < equationPNM.networkData.throatN; i++) {
        //     // std::cout << i << std::endl;
        //     for (int j = 0;
        //          j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
        //         std::cout << matrixConc[i][j] << ' ';
        //     }
        //     std::cout << std::endl;
        // }

        calcCoupledFlow();

        double pressInlet = 0;
        for (int i = 0; i < equationPNM.networkData.poreN; i++)
            if (equationPNM.networkData.poreLeftX[i])
                pressInlet += equationPNM.pressure[i];

        pressIn.emplace_back(pressInlet / boundPoresLeft.size());


        double sum = 0;
        for (int i = 0; i < equationPNM.networkData.poreN; i++) {
            sum += equationPNM.pressure[i];
        }
        pressureAverage.emplace_back(sum / equationPNM.networkData.poreN);


        sum = 0;
        for (int i = 0; i < equationPNM.networkData.throatN; i++) {
            sum = 0;
            for (int j = 0;
                 j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
                sum += matrixConc[i][j];
            }
            sum = sum / equationDiffusion.propsDiffusion.gridBlockN;

            cAV[i] = sum;
        }


        double concAv = accumulate(cAV.begin(), cAV.end(), 0.0) / cAV.size();
        concAverage.emplace_back(concAv);


        double diffFlowAv = accumulate(diffFlow.begin(), diffFlow.end(), 0.0) /
                            diffFlow.size();

        totalFlowDiff.emplace_back(diffFlowAv);

        std::cout << std::endl;

        std::cout << "Pressure" << std::endl;
        for (int i = 0; i < equationPNM.networkData.poreN; i++)
            std::cout << std::setw(8) << equationPNM.pressure[i];
        std::cout << std::endl;

        // std::cout << equationPNM.pressure[i] << std::endl;
        // std::cout << std::endl;
    }
}
//
// Getters for Python

const std::vector<double> DiffusionPNM::getPressureAverage() const {
    return pressureAverage;
}

const std::vector<double> DiffusionPNM::getConcAverage() const {
    return concAverage;
}

const std::vector<double> DiffusionPNM::getTotalFlowPoresOut() const {
    return totalFlowPoresOut;
}

const std::vector<double> DiffusionPNM::getTotalFlowPoresIn() const {
    return totalFlowPoresIn;
}

const std::vector<double> DiffusionPNM::getTotalFlowDiff() const {
    return totalFlowDiff;
}

const std::vector<double> DiffusionPNM::getInletPressure() const {
    return pressIn;
}

const std::vector<double> DiffusionPNM::getPorePressure() const {
    return equationPNM.pressure;
}






//        std::cout << "pressure" << std::endl;
//        for (int i = 0; i < equationPNM.networkData.poreN; i++)
//            std::cout << equationPNM.pressure[i] << std::endl;
//
//        std::cout << std::endl;
//

//        for (int i = 0; i < equationPNM.networkData.throatN; i++) {
//            std::cout << i << std::endl;
//            for (int j = 0; j < equation.propsDiffusion.gridBlockN; j++) {
//                std::cout << matrixConc[i][j] << ' ';
//            }
//            std::cout << std::endl;
//        }












