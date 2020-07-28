#include <DiffusionMath.h>
#include <numeric>
#include <iomanip>


DiffusionMath::DiffusionMath(PropsPnm &propsPnm, NetworkData &networkData,
                             EquationPNM &equationPNM,
                             EquationDiffusion &equationDiffusion,
                             const std::vector<double> &langmuirCoeff,
                             const double &matrixVolume) :

        propsPnm(propsPnm),
        networkData(networkData),
        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        langmuirCoeff(langmuirCoeff),
        matrixVolume(matrixVolume * 0.1),
        effRadius(networkData.fracturesN, 0),
        gridBlockN(equationDiffusion.propsDiffusion.gridBlockN),
        // TODO: don't forget to remove 0.1 castyl
        matrixWidth(networkData.fracturesN, 0),
        throatAvPress(networkData.fracturesN, 0),
        throatConc(networkData.fracturesN, 0),
        conc_ini(equationDiffusion.propsDiffusion.concIni),
        matricesOmega(networkData.fracturesN,
                      std::vector<double>(gridBlockN, 0)),
        matricesVolume(networkData.fracturesN,
                       std::vector<double>(gridBlockN, 0)) {}

double DiffusionMath::calcSideLength(std::vector<double> &poreCoord) {

    auto min = std::min_element(std::begin(poreCoord),
                                std::end(poreCoord));
    auto max = std::max_element(std::begin(poreCoord),
                                std::end(poreCoord));

    return *max - *min;
}

double DiffusionMath::calcDensConst() {

    auto pressureAv = (propsPnm.pressIn + propsPnm.pressOut) / 2;

    // TODO: use functional dependence
    return 1.986;
    // return propsPnm.aGasDens * pressureAv + propsPnm.bGasDens;
}

void DiffusionMath::calcRockVolume() {
    // ToDo hint: <= 0.0...01
    if (matrixVolume < 0) {

        auto lengthX = calcSideLength(networkData.poresCoordsX);
        auto lengthY = calcSideLength(networkData.poresCoordsY);
        auto lengthZ = calcSideLength(networkData.poresCoordsZ);

        matrixVolume = lengthX * lengthY * lengthZ;
    }
}

void DiffusionMath::calcEffRadius() {

    // TODO: connect effRadii to fracture area
    auto throatN = networkData.fracturesN;

    for (int i = 0; i < effRadius.size(); i++)
        effRadius[i] = matrixVolume / throatN;
}

void DiffusionMath::calcMatrixWidth() {

    auto throatN = networkData.fracturesN;

    calcEffRadius();

    for (int i = 0; i < throatN; i++) {

        auto fracHeight = networkData.fracturesHeights[i];
        auto fracLength = networkData.fracturesLengths[i];
        auto fracWidth = networkData.fracturesWidths[i];

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

    for (int i = 0; i < networkData.fracturesN; i++)
        throatAvPress[i] =
                (equationPNM.pressure[networkData.throatConns[i].first]
                 + equationPNM.pressure[networkData.throatConns[i].second]) / 2;
}

void DiffusionMath::calcThroatConc(const double &dP) {

    for (int i = 0; i < networkData.fracturesN; i++)
        throatConc[i] = calcLangmConc(throatAvPress[i] + dP);
}

void DiffusionMath::calcMatricesOmega() {

    for (int i = 0; i < networkData.fracturesN; i++) {

        auto fracHeight = networkData.fracturesHeights[i];
        auto fracLength = networkData.fracturesLengths[i];

        equationDiffusion.convectiveDiffusion.calcOmegaCartes(fracHeight,
                                                              fracLength);

        matricesOmega[i] = equationDiffusion.convectiveDiffusion.omegaCartesian;
    }
}

void DiffusionMath::calcMatricesVolume() {

    for (int i = 0; i < networkData.fracturesN; i++) {

        auto fracHeight = networkData.fracturesHeights[i];
        auto fracLength = networkData.fracturesLengths[i];
        auto fracWidth = networkData.fracturesWidths[i];


        equationDiffusion.localDiffusion.calcVolCartesian(fracHeight,
                                                          matrixWidth[i],
                                                          fracLength,
                                                          fracWidth);

        matricesVolume[i] = equationDiffusion.localDiffusion.volCartes;

    }
}
