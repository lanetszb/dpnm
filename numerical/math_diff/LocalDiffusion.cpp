#include <LocalDiffusion.h>
#include <cmath>

LocalDiffusion::LocalDiffusion(const std::vector<double> &propsVector,
                               const std::vector<double> &langmuirCoeff) :
        propsDiffusion(propsVector, langmuirCoeff),
        dRadius(0),
        alpha(propsDiffusion.gridBlockN, 0),
        radiusCurr(propsDiffusion.gridBlockN + 1, 0),
        volCylindr(propsDiffusion.gridBlockN, 0),
        volCartes(propsDiffusion.gridBlockN, 0) {}

int LocalDiffusion::left(const int &index) {
    return index;
}

int LocalDiffusion::right(const int &index) {
    return index + 1;
}

double LocalDiffusion::calcDelRadius(const double &radius,
                                     const double &effRadius,
                                     const int &gridBlockN) {

    return (effRadius - radius) / gridBlockN;
}

void LocalDiffusion::calcRadiusCurr(const double &radius,
                                    const double &effRadius,
                                    const int &gridBlockN) {

    dRadius = calcDelRadius(radius, effRadius, gridBlockN);

    for (int i = 0; i < gridBlockN + 1; i++)
        radiusCurr[i] = radius + i * dRadius;
}

void LocalDiffusion::calcVolCartesian(const double &frac_height,
                                      const double &matrix_width,
                                      const double &frac_length,
                                      const double &frac_width) {

    for (int i = 0; i < alpha.size(); i++)
        volCartes[i] = frac_length * frac_height * (matrix_width - frac_width) /
                       propsDiffusion.gridBlockN;
}

void LocalDiffusion::calcVolCylindr(const double &radius,
                                    const double &effRadius,
                                    const int &gridBlockN,
                                    const double &thrLength) {

    calcRadiusCurr(radius, effRadius, gridBlockN);

    for (int i = 0; i < alpha.size(); i++)
        volCylindr[i] = (M_PI * thrLength *
                         (radiusCurr[right(i)] * radiusCurr[right(i)] -
                          radiusCurr[left(i)] * radiusCurr[left(i)]));
}

void LocalDiffusion::calculateAlpha(const double &dt,
                                    const std::vector<double> &vol) {

    for (int i = 0; i < alpha.size(); i++)
        alpha[i] = vol[i] / dt;
}

const std::vector<double> LocalDiffusion::getRadCurr() const {
    return radiusCurr;
}

const std::vector<double> LocalDiffusion::getVolCylindr() const {
    return volCylindr;
}

const std::vector<double> LocalDiffusion::getVolCartes() const {
    return volCartes;
}

const std::vector<double> LocalDiffusion::getAlpha() const {
    return alpha;
}













