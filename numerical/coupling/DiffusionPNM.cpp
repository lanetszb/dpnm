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
                           const std::vector<double> &langmuirCoeff,
                           const double &matrixVolume) :

        equationPNM(propsPNM, throatList, throatHeight, throatLength,
                    throatWidth, connIndIn, connIndOut, poreCoordX, poreCoordY,
                    poreCoordZ, poreRadius, poreList, poreConns, connNumber,
                    porePerRow, poreLeftX, poreRightX, hydraulicCond),

        equationDiffusion(propsDiffusion),
        langmuirCoeff(langmuirCoeff),
        effRadius(equationPNM.networkData.throatN, 0),
        // TODO: don't forget to remove 0.1
        matrixVolume(matrixVolume * 0.1),
        matrixWidth(equationPNM.networkData.throatN, 0),
        throatAvPress(equationPNM.networkData.throatN, 0),
        throatConc(equationPNM.networkData.throatN, 0),
        conc_ini(equationDiffusion.propsDiffusion.concIni),
        matrixConc(equationPNM.networkData.throatN,
                   std::vector<double>(
                           equationDiffusion.propsDiffusion.gridBlockN, 0)),
        matricesOmega(equationPNM.networkData.throatN,
                      std::vector<double>(
                              equationDiffusion.propsDiffusion.gridBlockN, 0)),
        matricesVolume(equationPNM.networkData.throatN,
                       std::vector<double>(
                               equationDiffusion.propsDiffusion.gridBlockN, 0)),
        diffFlowInst(equationPNM.networkData.throatN, 0),
        diffFlowInstPlus(equationPNM.networkData.throatN, 0),
        diffFlowInstMinus(equationPNM.networkData.throatN, 0),
        flowDerivDiff(equationPNM.networkData.throatN, 0),
        connCoeffDiff(equationPNM.networkData.throatN, 0),
        centralCoeffDiff(equationPNM.networkData.poreN, 0),
        dP(0) {}

std::vector<std::vector<int>> DiffusionPNM::getGamma() {

    // calculate PN no diffusion with Direchlet-Newman for finding Gamma
    std::vector<bool> boundPoresInputForGamma;

    for (int i = 0; i < equationPNM.networkData.poreN; i++) {
        boundPoresInputForGamma.emplace_back(
                !equationPNM.networkData.poreRightX[i]);
    }

    std::vector<bool> poreLeftXSaved;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        poreLeftXSaved.push_back(equationPNM.networkData.poreLeftX[i]);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.networkData.poreLeftX[i] = boundPoresInputForGamma[i];

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.porFlowRate[i] = 1.e-12;

    equationPNM.cfdProcedure("mixed",
                             equationPNM.networkData.poreRightX,
                             equationPNM.pIn,
                             equationPNM.pOut);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.networkData.poreLeftX[i] = poreLeftXSaved[i];

    std::vector<std::vector<int>> poreConnsIsOutByPressureSaved(
            equationPNM.networkData.poreN);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        for (int j = 0; j < equationPNM.gammaPnm[i].size(); j++)
            poreConnsIsOutByPressureSaved[i].push_back(
                    equationPNM.gammaPnm[i][j]);

    return poreConnsIsOutByPressureSaved;
}

void DiffusionPNM::getInletFlow() {

    // calculate PN no diffusion with Direchlet
    std::vector<bool> boundPoresInput;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        boundPoresInput.emplace_back(equationPNM.networkData.poreLeftX[i] or
                                     equationPNM.networkData.poreRightX[i]);

    equationPNM.cfdProcedure("dirichlet",
                             boundPoresInput,
                             equationPNM.pIn,
                             equationPNM.pOut);

    // calculate PN no diffusion with mixed Dirichlet-Newman

    equationPNM.cfdProcedure("mixed",
                             equationPNM.networkData.poreRightX,
                             equationPNM.pIn,
                             equationPNM.pOut);

}

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

    // TODO: use functional dependence
    return 1.986;
    // return aGasDens * pressureAv + bGasDens;
}

void DiffusionPNM::calcRockVolume() {

    if (matrixVolume <= 0.) {

        auto lengthX = calcSideLength(equationPNM.networkData.poreCoordX);
        auto lengthY = calcSideLength(equationPNM.networkData.poreCoordY);
        auto lengthZ = calcSideLength(equationPNM.networkData.poreCoordZ);

        matrixVolume = lengthX * lengthY * lengthZ;
    }
}

void DiffusionPNM::calcEffRadius() {

    // TODO: connect effRadii to fracture area
    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < effRadius.size(); i++)
        effRadius[i] = matrixVolume / throatN;
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

void DiffusionPNM::calcMatricesOmega() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {

        auto gridBlockN = equationDiffusion.propsDiffusion.gridBlockN;
        auto fracHeight = equationPNM.networkData.throatRadius[i];
        auto fracLength = equationPNM.networkData.throatLength[i];

        equationDiffusion.convectiveDiffusion.calcOmegaCartes(fracHeight,
                                                              fracLength);

        for (int j = 0; j < gridBlockN; j++) {
            auto omega = equationDiffusion.convectiveDiffusion.omegaCartesian[j];
            matricesOmega[i][j] = omega;
        }
    }
}

void DiffusionPNM::calcMatricesVolume() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {

        auto gridBlockN = equationDiffusion.propsDiffusion.gridBlockN;
        auto fracHeight = equationPNM.networkData.throatRadius[i];
        auto fracLength = equationPNM.networkData.throatLength[i];
        auto fracWidth = equationPNM.networkData.throatWidth[i];


        equationDiffusion.localDiffusion.calcVolCartesian(fracHeight,
                                                          matrixWidth[i],
                                                          fracLength,
                                                          fracWidth);

        for (int j = 0; j < gridBlockN; j++) {
            auto volume = equationDiffusion.localDiffusion.volCartes[j];
            matricesVolume[i][j] = volume;
        }
    }
}

void DiffusionPNM::calcDiffFlow(std::vector<double> &diffFlowVector) {

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {

        auto gridBlockN = equationDiffusion.propsDiffusion.gridBlockN;

        for (int j = 0; j < gridBlockN; j++) {
            equationDiffusion.conc[0][j] = matrixConc[i][j];
            equationDiffusion.conc[1][j] = matrixConc[i][j];
        }

        equationDiffusion.cfdProcedureOneStep(
                throatConc[i],
                equationPNM.networkData.throatRadius[i],
                matrixWidth[i],
                equationPNM.networkData.throatLength[i],
                matricesVolume[i],
                matricesOmega[i]);

        // diffFlowVector[i] = equationDiffusion.flowRate / densityConst;

        double flowSum = 0;
        for (int j = 0; j < gridBlockN; j++) {
            auto conc_curr = equationDiffusion.conc[equationDiffusion.iCurr][j];
            auto conc_prev = equationDiffusion.conc[equationDiffusion.iPrev][j];

            flowSum += -1 * (conc_curr - conc_prev) *
                       equationDiffusion.localDiffusion.volCartes[j] /
                       equationDiffusion.propsDiffusion.timeStep;
        }

        diffFlowVector[i] = flowSum / densityConst;

        for (int j = 0; j < gridBlockN; j++) {
            auto conc_prev = equationDiffusion.conc[equationDiffusion.iPrev][j];
            matrixConc[i][j] = conc_prev;
        }
    }
}

void DiffusionPNM::updateConc() {

    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < throatN; i++) {

        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            auto conc_curr = equationDiffusion.conc[equationDiffusion.iCurr][j];
            matrixConc[i][j] = conc_curr;
        }
    }
}

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

            coeffSum += equationPNM.gammaPnm[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] / 2;

        // TODO: understand plus or minus sign and properly name as derivDiff
        centralCoeffDiff[i] = 0;
        // centralCoeffDiff[i] = -1 * coeffSum;
        // centralCoeffDiff[i] = coeffSum;
    }
}

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

            coeffSum += equationPNM.gammaPnm[i][j] *
                        diffFlowInst[equationPNM.porConns[i][j]];
        }
        porFlowDiff[i] = coeffSum;
    }

    for (int i = 0; i < equationPNM.porConns.size(); i++) {
        double coeffSum = 0;
        for (int j = 0; j < equationPNM.porConns[i].size(); j++) {

            coeffSum -= equationPNM.gammaPnm[i][j] *
                        flowDerivDiff[equationPNM.porConns[i][j]] *
                        throatAvPress[equationPNM.porConns[i][j]];
        }
        porFlowDiffDer[i] = 0;
        // porFlowDiffDer[i] = coeffSum;
    }

    equationPNM.calculateFreeVector("mixed",
                                    equationPNM.propsPNM.pressIn,
                                    equationPNM.propsPNM.pressOut);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        // TODO: understand plus or minus sign
        if (!equationPNM.networkData.poreRightX[i])
            equationPNM.freeVector[i] += porFlowDiff[i];
    // equationPNM.freeVector[i] += porFlowDiff[i] + porFlowDiffDer[i];
    // equationPNM.freeVector[i] += -1 * (porFlowDiff[i] - porFlowDiffDer[i]);
}

void DiffusionPNM::setInitialCondCoupledMod() {

    densityConst = calcDensConst();

    // equationPNM.calcPorConns();
    equationPNM.calcThroatConns();

    calcRockVolume();
    calcMatrixWidth();
    calcMatricesOmega();
    calcMatricesVolume();

    // equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
    //                                equationPNM.propsPNM.pressOut);

    equationDiffusion.calcConcIni(conc_ini);

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++)
            matrixConc[i][j] = conc_ini;
}

void DiffusionPNM::calcMatrixMassTot() {

    double sum = 0;
    std::vector<double> matrixMass;
    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
        sum = 0;
        for (int j = 0;
             j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            sum += conc_ini *
                   equationDiffusion.localDiffusion.volCartes[j];
        }
        matrixMass.emplace_back(sum);
    }

    matrixMassTotal.emplace_back(
            accumulate(matrixMass.begin(), matrixMass.end(), 0.0));
}

void DiffusionPNM::calcVecSum(const int &iter,
                              const std::vector<double> &vectorToSum,
                              std::vector<double> &vectorSum,
                              const double &mult) {

    double sum = 0;
    for (int i = 0; i < iter; i++)
        sum += vectorToSum[i];

    vectorSum.emplace_back(sum * mult);
}

void DiffusionPNM::calcPressInlet() {

    // TODO: Should be already known, rather than calculated every time
    double boundPoresLeftSize = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        if (equationPNM.networkData.poreLeftX[i])
            boundPoresLeftSize += 1;

    double pressInlet = 0;
    for (int i = 0; i < equationPNM.networkData.poreN; i++) {
        if (equationPNM.networkData.poreLeftX[i])
            pressInlet += equationPNM.pressure[i];
    }

    pressIn.emplace_back(pressInlet / boundPoresLeftSize);
}

void DiffusionPNM::calcDiffPart() {

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
    updateConc();
}

void DiffusionPNM::solveCoupledMatrix() {

    calcMatCoeffDiff();
    calcMatCoupledCoeff();

    equationPNM.calculateMatrix(
            equationPNM.connCoeff,
            equationPNM.centralCoeff,
            equationPNM.networkData.poreRightX,
            equationPNM.gammaPnm,
            connCoeffDiff);

    calcCoupledFreeVector();
    equationPNM.calculateGuessVector();
    // TODO: the solution with guess vector should be preffered later
    equationPNM.calculatePress(1);
}

void DiffusionPNM::calcCoupledFlowParams() {

    // Total flow from diffusion
    calcVecSum(equationPNM.networkData.throatN, diffFlowInst,
               totalFlowDiff, densityConst);
    // diffFlowThroat += diffFlowInst[i] - flowDerivDiff[i] * throatAvPress[i];

    // equationPNM.getGammaByPressure();

    equationPNM.calcThrFlowRate();
    equationPNM.calcPorFlowRate();
    equationPNM.calcTotFlow(equationPNM.networkData.poreLeftX);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate * densityConst);
}

void DiffusionPNM::calcCoupledFlow() {

    calcDiffPart();
    solveCoupledMatrix();
    calcCoupledFlowParams();
}

void DiffusionPNM::cfdProcedurePnmDiff() {

    equationPNM.setInitialCondPurePnm();
    std::vector<std::vector<int>> gammaByPressureSaved = getGamma();
    getInletFlow();

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        for (int j = 0; j < equationPNM.gammaPnm[i].size(); j++)
            equationPNM.gammaPnm[i][j] = gammaByPressureSaved[i][j];

    setInitialCondCoupledMod();
    calcVecSum(equationPNM.networkData.poreN, equationPNM.pressure,
               pressureAv, 1.0 / equationPNM.networkData.poreN);
    calcMatrixMassTot();

    totalFlowDiff.emplace_back(0);
    totalFlowPoresOut.emplace_back(equationPNM.totFlowRate * densityConst);
    totalFlowPoresIn.emplace_back(equationPNM.totFlowRate * densityConst);

    calcPressInlet();

    // TODO: Remove castyl
    for (double t = equationDiffusion.propsDiffusion.timeStep;
         t < equationDiffusion.propsDiffusion.time * (1. + 1.e-3);
         t += equationDiffusion.propsDiffusion.timeStep) {

        calcCoupledFlow();
        calcPressInlet();
        calcVecSum(equationPNM.networkData.poreN, equationPNM.pressure,
                   pressureAv, 1.0 / equationPNM.networkData.poreN);
        calcMatrixMassTot();
    }
}

// Getters for Python
const std::vector<double> DiffusionPNM::getPressureAverage() const {
    return pressureAv;
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












