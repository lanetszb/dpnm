#include <DiffusionMath.h>
#include <numeric>
#include <iomanip>


DiffusionMath::DiffusionMath(PropsPNM &propsPnm, NetworkData &networkData,
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
        effRadius(networkData.throatN, 0),
        gridBlockN(equationDiffusion.propsDiffusion.gridBlockN),
        // TODO: don't forget to remove 0.1 castyl
        matrixWidth(networkData.throatN, 0),
        throatAvPress(networkData.throatN, 0),
        throatConc(networkData.throatN, 0),
        conc_ini(equationDiffusion.propsDiffusion.concIni),
        matricesOmega(networkData.throatN,
                      std::vector<double>(gridBlockN, 0)),
        matricesVolume(networkData.throatN,
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
    if (matrixVolume <= 0.00000001) {

        auto lengthX = calcSideLength(networkData.poreCoordX);
        auto lengthY = calcSideLength(networkData.poreCoordY);
        auto lengthZ = calcSideLength(networkData.poreCoordZ);

        matrixVolume = lengthX * lengthY * lengthZ;
    }
}

void DiffusionMath::calcEffRadius() {

    // TODO: connect effRadii to fracture area
    auto throatN = networkData.throatN;

    for (int i = 0; i < effRadius.size(); i++)
        effRadius[i] = matrixVolume / throatN;
}

void DiffusionMath::calcMatrixWidth() {

    auto throatN = networkData.throatN;

    calcEffRadius();

    for (int i = 0; i < throatN; i++) {

        auto fracHeight = networkData.throatRadius[i];
        auto fracLength = networkData.throatLength[i];
        auto fracWidth = networkData.throatWidth[i];

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
    // ToDo provide pressure in this class without EquationPNM
    // ToDo transfer throatConns from equationPNM to NetworkData
    for (int i = 0; i < networkData.throatN; i++)
        throatAvPress[i] =
                (equationPNM.pressure[equationPNM.throatConns[i].first]
                 + equationPNM.pressure[equationPNM.throatConns[i].second]) / 2;
}

void DiffusionMath::calcThroatConc(const double &dP) {

    for (int i = 0; i < networkData.throatN; i++)
        throatConc[i] = calcLangmConc(throatAvPress[i] + dP);
}

void DiffusionMath::calcMatricesOmega() {

    for (int i = 0; i < networkData.throatN; i++) {

        auto fracHeight = networkData.throatRadius[i];
        auto fracLength = networkData.throatLength[i];

        equationDiffusion.convectiveDiffusion.calcOmegaCartes(fracHeight,
                                                              fracLength);

        for (int j = 0; j < gridBlockN; j++) {
            auto omega = equationDiffusion.convectiveDiffusion.omegaCartesian[j];
            matricesOmega[i][j] = omega;
        }
    }
}

void DiffusionMath::calcMatricesVolume() {

    for (int i = 0; i < networkData.throatN; i++) {

        auto fracHeight = networkData.throatRadius[i];
        auto fracLength = networkData.throatLength[i];
        auto fracWidth = networkData.throatWidth[i];


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
