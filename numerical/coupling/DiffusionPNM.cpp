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
        diffFlowInst(equationPNM.networkData.throatN, 0),
        diffFlowInstPlus(equationPNM.networkData.throatN, 0),
        diffFlowInstMinus(equationPNM.networkData.throatN, 0),
        flowDerivDiff(equationPNM.networkData.throatN, 0),
        connCoeffDiff(equationPNM.networkData.throatN, 0),
        centralCoeffDiff(equationPNM.networkData.poreN, 0),
        conc_ini(0.77),
        dP(0) {

    equationPNM.setInitialCond();

    // calculate PN no diffusion with Direchlet-Newman for finding Gamma
    std::vector<bool> boundPoresInputForGamma;

    for (int i = 0; i < equationPNM.networkData.poreN; i++) {
        if (equationPNM.networkData.poreRightX[i])
            boundPoresInputForGamma.emplace_back(false);
        else boundPoresInputForGamma.emplace_back(true);
    }

    std::vector<bool> poreLeftXSaved;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        poreLeftXSaved.push_back(equationPNM.networkData.poreLeftX[i]);

    std::cout << "Boundary pores LeftX" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(12) << equationPNM.networkData.poreLeftX[i];
    std::cout << std::endl;


    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.networkData.poreLeftX[i] = boundPoresInputForGamma[i];


    std::cout << "Boundary pores In" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(12) << boundPoresInputForGamma[i];
    std::cout << std::endl;

    std::vector<int> boundPoresRight;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreRightX[i])
            boundPoresRight.emplace_back(i);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.porFlowRate[i] = 1.e-12;

    equationPNM.cfdProcedure(0,
                             boundPoresRight,
                             equationPNM.pIn,
                             equationPNM.pOut);

    std::cout << "Pressure for Gamma" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(12) << equationPNM.pressure[i];
    std::cout << std::endl;


    std::cout << "Flowrate for Gamma" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(12) << equationPNM.porFlowRate[i];
    std::cout << std::endl;


    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.networkData.poreLeftX[i] = poreLeftXSaved[i];

    std::vector<std::vector<int>> poreConnsIsOutByPressureSaved(
            equationPNM.networkData.poreN);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        for (int j = 0; j < equationPNM.porConnsIsOutByPressure[i].size(); j++)
            poreConnsIsOutByPressureSaved[i].push_back(
                    equationPNM.porConnsIsOutByPressure[i][j]);

    // calculate PN no diffusion with Direchlet
    std::vector<int> boundPoresInput;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i] or
            equationPNM.networkData.poreRightX[i])
            boundPoresInput.emplace_back(i);

    std::cout << "boundPoresInput" << std::endl;
    for (int i = 0; i < boundPoresInput.size(); i++)
        std::cout << std::setw(12) << boundPoresInput[i];
    std::cout << std::endl;


    equationPNM.cfdProcedure(1,
                             boundPoresInput,
                             equationPNM.pIn,
                             equationPNM.pOut);

    std::cout << "Pressure Dirichlet" << std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(9) << equationPNM.pressure[i];
    std::cout << std::endl;


    // calculate PN no diffusion with mixed Dirichlet-Newman
    std::cout << "inletFlow" << std::endl;

    for (int i = 0; i < equationPNM.inletFlow.size(); i++)
        std::cout << equationPNM.inletFlow[i] << std::endl;

    equationPNM.cfdProcedure(0,
                             boundPoresRight,
                             equationPNM.pIn,
                             equationPNM.pOut);

    std::cout << "Pressure Dirichlet+Neumann" <<
              std::endl;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        std::cout << std::setw(9) << equationPNM.pressure[i];
    std::cout << std::endl;


    // calculate PN with diffusion with mixed Dirichlet-Newman


    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        for (int j = 0; j < equationPNM.porConnsIsOutByPressure[i].size(); j++)
            equationPNM.porConnsIsOutByPressure[i][j] = poreConnsIsOutByPressureSaved[i][j];

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

    return 1.986;
    // return aGasDens * pressureAv + bGasDens;
}

void DiffusionPNM::calcRockVolume() {

    auto lengthX = calcSideLength(equationPNM.networkData.poreCoordX);
    auto lengthY = calcSideLength(equationPNM.networkData.poreCoordY);
    auto lengthZ = calcSideLength(equationPNM.networkData.poreCoordZ);

    // TODO: to make the option of imporing real rock volume like is done below
    // rockVolume = (lengthX * lengthY * lengthZ);
    rockVolume = 0.018 * 0.018 * 0.018 * 0.1;

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

        // diffFlowVector[i] = equationDiffusion.flowRate / densityConst;

        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++)
            matrixConc[i][j] = equationDiffusion.conc[equationDiffusion.iPrev][j];

        double flowSum = 0;
        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++)
            flowSum +=
                    -1 * (equationDiffusion.conc[equationDiffusion.iCurr][j] -
                          equationDiffusion.conc[equationDiffusion.iPrev][j]) *
                    equationDiffusion.localDiffusion.volCartes[j] /
                    equationDiffusion.propsDiffusion.timeStep;

        diffFlowVector[i] = flowSum / densityConst;
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
        flowDerivDiff[i] = ((diffFlowInstPlus[i] - diffFlowInstMinus[i]) / dP);
}

void DiffusionPNM::calcMatCoeffDiff() {

    for (int i = 0; i < flowDerivDiff.size(); i++)
        connCoeffDiff[i] = 0;
    // connCoeffDiff[i] = flowDerivDiff[i] / 2;


    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++)

            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] / 2;

        // TODO: understand plus or minus sign and properly name as derivDiff
        centralCoeffDiff[i] = 0;
        // centralCoeffDiff[i] = -1 * coeffSum;
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
        std::cout << "i: " << i << std::endl;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            // TODO: understand plus or minus sign and properly name as derivDiff
            // coeffSum -= equationPNM.porConnsIsOutByPressure[i][j] *
            //            diffFlowInst[equationPNM.porConns[i][j]];
            coeffSum += equationPNM.porConnsIsOutByPressure[i][j] *
                        diffFlowInst[equationPNM.porConns[i][j]];
            std::cout << std::endl << "gamma: "
                      << equationPNM.porConnsIsOutByPressure[i][j]
                      << " " << equationPNM.porConns[i][j] << std::endl;
            if (equationPNM.porConnsIsOutByPressure[i][j] == 1)
                std::cout << equationPNM.porConns[i][j] << std::endl;
        }
        porFlowDiff[i] = coeffSum;
    }

    std::cout << "porFlowDiffAcum: "
              << accumulate(porFlowDiff.begin(), porFlowDiff.end(),
                            0.0) << std::endl;

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
        porFlowDiffDer[i] = 0;
        // porFlowDiffDer[i] = coeffSum;
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

    // for (int i = 0; i < equationPNM.networkData.boundaryPoresIn.size(); i++)
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        // TODO: understand plus or minus sign
        if (!equationPNM.networkData.poreRightX[i])
            equationPNM.freeVector[i] += porFlowDiff[i];
            // equationPNM.freeVector[i] += porFlowDiff[i] + porFlowDiffDer[i];
    // equationPNM.freeVector[i] += -1 * (porFlowDiff[i] - porFlowDiffDer[i]);
}

//
void DiffusionPNM::setInitialCond() {

    auto throatN = equationPNM.networkData.throatN;

    densityConst = calcDensConst();

    // equationPNM.calcPorConns();
    equationPNM.calcThroatConns();

    calcRockVolume();
    calcMatrixWidth();

    // equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
    //                                equationPNM.propsPNM.pressOut);

    equationDiffusion.calcConcIni(conc_ini);

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++)
            matrixConc[i][j] = conc_ini;

    // equationPNM.getPorConnsIsOut();
}

void DiffusionPNM::calcCoupledFlow() {

    calcThroatAvPress();

    // Enhance and rethink later
    dP = 0;

    calcThroatConc(dP);
    calcDiffFlow(diffFlowInst);

    dP = (equationPNM.pIn - equationPNM.pOut) / 10.e+5;
    calcThroatConc(dP / 2);
    calcDiffFlow(diffFlowInstPlus);

    calcThroatConc(-1 * dP / 2);
    calcDiffFlow(diffFlowInstMinus);

    calcDiffFlowDeriv();

    for (int i = 0; i < diffFlowInst.size(); i++)
        std::cout << diffFlowInst[i] << std::endl;
    std::cout << std::endl;

    updateConc();

    // Made To calculate diff flow properly
    double diffFlowThroat = 0;

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
        diffFlowThroat += diffFlowInst[i];
        std::cout << i << std::endl;
    }
    totalFlowDiff.emplace_back(diffFlowThroat * densityConst);
    // diffFlowThroat += diffFlowInst[i] - flowDerivDiff[i] * throatAvPress[i];

    std::cout << "porFlowDiffSecond: "
              << accumulate(totalFlowDiff.begin(), totalFlowDiff.end(),
                            0.0) << std::endl;
    std::cout << "Throat9: " << diffFlowInst[9] << std::endl;

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

    // Made to calculate diff flow properly

    // for (int i = 0; i < equationPNM.networkData.throatN; i++)
    //    diffFlowThroat += diffFlowInst[i];
    //     diffFlowThroat += flowDerivDiff[i] * throatAvPress[i];

    // equationPNM.getPorConnsIsOutByPressure();

    equationPNM.calcThrFlowRate();
    equationPNM.calcPorFlowRate();


    equationPNM.calcTotFlow(boundPoresRight);


    // for calculating diffusive flow to outlet throat
    for (int i = 0; i < equationPNM.porFlowRate.size(); i++)
        if (equationPNM.networkData.poreRightX[i])
            for (int j = 0; j < equationPNM.porConns[i].size(); j++) {
                equationPNM.totFlowRate +=
                        -1. * diffFlowInst[equationPNM.porConns[i][j]];
            }
    // finished

    totalFlowPoresOut.emplace_back(-equationPNM.totFlowRate * densityConst);

    // totalFlowPoresOut.emplace_back(
    //         -1 * (equationPNM.totFlowRate-diffFlowInst[9]) * densityConst);

    equationPNM.calcTotFlow(boundPoresLeft);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate * densityConst);

    // updateConc();

}

//
void DiffusionPNM::cfdProcedureDiff() {

    setInitialCond();

    double sum = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        sum += equationPNM.pressure[i];

    pressureAverage.emplace_back(sum / equationPNM.networkData.poreN);


    // ZATYCHKA
    sum = 0;
    std::vector<double> matrixMass;
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
        sum = 0;
        for (int j = 0;
             j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            sum += conc_ini *
                   equationDiffusion.localDiffusion.volCartes[j];
        }
//        sum = sum / equationDiffusion.propsDiffusion.gridBlockN;
        // sum = sum;

        matrixMass.emplace_back(sum);
    }

    matrixMassTotal.emplace_back(
            accumulate(matrixMass.begin(), matrixMass.end(), 0.0));
    // ZATYCHKA


    // matrixMassTotal.emplace_back(conc_ini);
    totalFlowDiff.emplace_back(0);

    totalFlowPoresOut.emplace_back(equationPNM.totFlowRate * densityConst);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate * densityConst);

    // TODO: inlet pressure to be rewritten later
    std::vector<int> boundPoresLeft;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i])
            boundPoresLeft.emplace_back(i);

    double pressInlet = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++) {
        std::cout << equationPNM.pressure[i] << std::endl;
        if (equationPNM.networkData.poreLeftX[i])
            pressInlet += equationPNM.pressure[i];
    }

    std::cout << "pressInlet" << pressInlet << std::endl;

    pressIn.emplace_back(pressInlet / boundPoresLeft.size());

    // To rewrite latter
    //for (int i = 0; i < equationPNM.porConnsIsOut.size(); i++)
    //    for (int j = 0; j < equationPNM.porConns[i].size(); j++) {
    //        std::cout << equationPNM.porConnsIsOutByPressure[i][j] << std::endl;
    //    }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    //    }


    for (double t = equationDiffusion.propsDiffusion.timeStep;
         t < equationDiffusion.propsDiffusion.time * (1. + 1.e-3);
         t += equationDiffusion.propsDiffusion.timeStep) {

        std::cout << "Iteration: " << t << std::endl;

        //std::cout << "freeVector" << std::endl;
        //std::cout << equationPNM.freeVector << std::endl;;
        //std::cout << std::endl;

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


        for (int i = 0; i < equationPNM.networkData.throatN; i++)
            matrixMass[i] = 0;

        sum = 0;
        for (int i = 0; i < equationPNM.networkData.throatN; i++) {
            sum = 0;
            for (int j = 0;
                 j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
                sum += matrixConc[i][j] *
                       equationDiffusion.localDiffusion.volCartes[j];
            }
            sum = sum;

            matrixMass[i] = sum;
        }

        matrixMassTotal.emplace_back(
                accumulate(matrixMass.begin(), matrixMass.end(), 0.0));


        // double diffFlowAv =
        //         accumulate(diffFlowInst.begin(), diffFlowInst.end(), 0.0) /
        //         diffFlowInst.size();

        // totalFlowDiff.emplace_back(diffFlowAv);

        // std::cout << std::endl;
        //
        // std::cout << "Pressure" << std::endl;
        // for (int i = 0; i < equationPNM.networkData.poreN; i++)
        //     std::cout << std::setw(8) << equationPNM.pressure[i];
        // std::cout << std::endl;

        // std::cout << equationPNM.pressure[i] << std::endl;
        // std::cout << std::endl;
    }
}
//
// Getters for Python

const std::vector<double> DiffusionPNM::getPressureAverage() const {
    return pressureAverage;
}

const std::vector<double> DiffusionPNM::getMatrixMassTotal() const {
    return matrixMassTotal;
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












