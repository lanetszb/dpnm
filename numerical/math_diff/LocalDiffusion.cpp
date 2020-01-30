#include <LocalDiffusion.h>
#include <cmath>

LocalDiffusion::LocalDiffusion(const std::vector<double> &propsVector,
                               const std::vector<double> &langmuirCoeff) :
        propsDiffusion(propsVector, langmuirCoeff),
        dRadius(0),
        alpha(propsDiffusion.gridBlockN, 0),
        radiusCurr(propsDiffusion.gridBlockN + 1, 0) {}

int LocalDiffusion::left(const int &index) {
    return index;
}

int LocalDiffusion::right(const int &index) {
    return index + 1;
}

double LocalDiffusion::calcDelRadius(const double &radius, const double &effRadius,
                                     const int &gridBlockN) {

    return (effRadius - radius) / gridBlockN;
}

void LocalDiffusion::calcRadiusCurr(const double &radius, const double &effRadius,
                                    const int &gridBlockN) {

    dRadius = calcDelRadius(radius, effRadius, gridBlockN);

    for (int i = 0; i < gridBlockN + 1; i++)
        radiusCurr[i] = radius + i * dRadius;
}

void LocalDiffusion::calculateAlpha(const double &dt,
                                    const double &radius,
                                    const double &effRadius,
                                    const double &thrLength) {

    calcRadiusCurr(radius, effRadius, propsDiffusion.gridBlockN);

    for (int i = 0; i < alpha.size(); i++)
        alpha[i] = (M_PI * thrLength *
                    (radiusCurr[right(i)] * radiusCurr[right(i)] -
                     radiusCurr[left(i)] * radiusCurr[left(i)])) / dt;
}
const std::vector<double> LocalDiffusion::getAlpha() const {
    return alpha;
}

const std::vector<double> LocalDiffusion::getRadCurr() const {
    return radiusCurr;
}











