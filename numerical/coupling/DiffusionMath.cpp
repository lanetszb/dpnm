#include <DiffusionMath.h>
#include <numeric>
#include <iomanip>


DiffusionMath::DiffusionMath(EquationPNM &equationPNM,
                             EquationDiffusion &equationDiffusion,
                             const std::vector<double> &langmuirCoeff,
                             const double &matrixVolume) :

        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        langmuirCoeff(langmuirCoeff),
        matrixVolume(matrixVolume * 0.1),
        effRadius(equationPNM.networkData.throatN, 0),
        // TODO: don't forget to remove 0.1 castyl
        matrixWidth(equationPNM.networkData.throatN, 0),
        throatAvPress(equationPNM.networkData.throatN, 0),
        throatConc(equationPNM.networkData.throatN, 0),
        conc_ini(equationDiffusion.propsDiffusion.concIni),
        matricesOmega(equationPNM.networkData.throatN,
                      std::vector<double>(
                              equationDiffusion.propsDiffusion.gridBlockN, 0)),
        matricesVolume(equationPNM.networkData.throatN,
                       std::vector<double>(
                               equationDiffusion.propsDiffusion.gridBlockN,
                               0)) {}

double DiffusionMath::calcSideLength(std::vector<double> &poreCoord) {

    auto min = std::min_element(std::begin(poreCoord),
                                std::end(poreCoord));
    auto max = std::max_element(std::begin(poreCoord),
                                std::end(poreCoord));

    return *max - *min;
}

double DiffusionMath::calcDensConst() {

    auto aGasDens = equationPNM.propsPNM.aGasDens;
    auto bGasDens = equationPNM.propsPNM.bGasDens;

    auto pressIn = equationPNM.propsPNM.pressIn;
    auto pressOut = equationPNM.propsPNM.pressOut;

    auto pressureAv = (pressIn + pressOut) / 2;

    // TODO: use functional dependence
    return 1.986;
    // return aGasDens * pressureAv + bGasDens;
}

void DiffusionMath::calcRockVolume() {

    if (matrixVolume <= 0.) {

        auto lengthX = calcSideLength(equationPNM.networkData.poreCoordX);
        auto lengthY = calcSideLength(equationPNM.networkData.poreCoordY);
        auto lengthZ = calcSideLength(equationPNM.networkData.poreCoordZ);

        matrixVolume = lengthX * lengthY * lengthZ;
    }
}

void DiffusionMath::calcEffRadius() {

    // TODO: connect effRadii to fracture area
    auto throatN = equationPNM.networkData.throatN;

    for (int i = 0; i < effRadius.size(); i++)
        effRadius[i] = matrixVolume / throatN;
}

void DiffusionMath::calcMatrixWidth() {

    auto throatN = equationPNM.networkData.throatN;

    calcEffRadius();

    for (int i = 0; i < throatN; i++) {

        auto fracHeight = equationPNM.networkData.throatRadius[i];
        auto fracLength = equationPNM.networkData.throatLength[i];
        auto fracWidth = equationPNM.networkData.throatWidth[i];

        matrixWidth[i] = effRadius[i] / fracLength / fracHeight + fracWidth;
    }
}

double DiffusionMath::calcLangmConc(double pressure) {

    langmConc = 0;
    for (int i = 0; i < langmuirCoeff.size(); i++)
        langmConc += langmuirCoeff[i] * pow(pressure, i);

    return langmConc;
}

void DiffusionMath::calcThroatAvPress() {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        throatAvPress[i] =
                (equationPNM.pressure[equationPNM.throatConns[i].first]
                 + equationPNM.pressure[equationPNM.throatConns[i].second]) / 2;
}

void DiffusionMath::calcThroatConc(const double &dP) {

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        throatConc[i] = calcLangmConc(throatAvPress[i] + dP);
}

void DiffusionMath::calcMatricesOmega() {

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

void DiffusionMath::calcMatricesVolume() {

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
