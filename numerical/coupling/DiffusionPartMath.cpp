#include <DiffusionPartMath.h>
#include <numeric>
#include <iomanip>


DiffusionPartMath::DiffusionPartMath(const std::vector<double> &propsPNM,
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
        // TODO: don't forget to remove 0.1 castyl
        matrixVolume(matrixVolume * 0.1),
        matrixWidth(equationPNM.networkData.throatN, 0),
        throatAvPress(equationPNM.networkData.throatN, 0),
        throatConc(equationPNM.networkData.throatN, 0),
        matricesOmega(equationPNM.networkData.throatN,
                      std::vector<double>(
                              equationDiffusion.propsDiffusion.gridBlockN, 0)),
        matricesVolume(equationPNM.networkData.throatN,
                       std::vector<double>(
                               equationDiffusion.propsDiffusion.gridBlockN,
                               0)) {}

double DiffusionPartMath::calcSideLength(std::vector<double> &poreCoord) {

    auto min = std::min_element(std::begin(poreCoord),
                                std::end(poreCoord));
    auto max = std::max_element(std::begin(poreCoord),
                                std::end(poreCoord));

    return *max - *min;
}

double DiffusionPartMath::calcDensConst() {

    auto aGasDens = equationPNM.propsPNM.aGasDens;
    auto bGasDens = equationPNM.propsPNM.bGasDens;

    auto pressIn = equationPNM.propsPNM.pressIn;
    auto pressOut = equationPNM.propsPNM.pressOut;

    auto pressureAv = (pressIn + pressOut) / 2;

    // TODO: use functional dependence
    return 1.986;
    // return aGasDens * pressureAv + bGasDens;
}

void DiffusionPartMath::calcRockVolume() {

    if (matrixVolume <= 0.) {

        auto lengthX = calcSideLength(equationPNM.networkData.poreCoordX);
        auto lengthY = calcSideLength(equationPNM.networkData.poreCoordY);
        auto lengthZ = calcSideLength(equationPNM.networkData.poreCoordZ);

        matrixVolume = lengthX * lengthY * lengthZ;
    }
}

void DiffusionPartMath::calcEffRadius() {

    // TODO: connect effRadii to fracture area
    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < effRadius.size(); i++)
        effRadius[i] = matrixVolume / throatN;
}

void DiffusionPartMath::calcMatrixWidth() {

    auto throatN = equationPNM.networkData.throatN;

    calcEffRadius();

    for (int i = 0; i < throatN; i++) {

        auto fracHeight = equationPNM.networkData.throatRadius[i];
        auto fracLength = equationPNM.networkData.throatLength[i];
        auto fracWidth = equationPNM.networkData.throatWidth[i];

        matrixWidth[i] = effRadius[i] / fracLength / fracHeight + fracWidth;
    }
}

double DiffusionPartMath::calcLangmConc(double pressure) {

    langmConc = 0;
    for (int i = 0; i < langmuirCoeff.size(); i++)
        langmConc += langmuirCoeff[i] * pow(pressure, i);

    return langmConc;
}

void DiffusionPartMath::calcThroatAvPress() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        throatAvPress[i] =
                (equationPNM.pressure[equationPNM.throatConns[i].first]
                 + equationPNM.pressure[equationPNM.throatConns[i].second]) / 2;
}

void DiffusionPartMath::calcThroatConc(const double &dP) {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        throatConc[i] = calcLangmConc(throatAvPress[i] + dP);
}

void DiffusionPartMath::calcMatricesOmega() {

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

void DiffusionPartMath::calcMatricesVolume() {

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