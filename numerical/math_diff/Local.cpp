#include <Local.h>
#include <cmath>

Local::Local(const std::vector<double> &propsVector,
             const std::vector<double> &langmuirCoeff) :
        props(propsVector, langmuirCoeff),
        dRadius(0),
        alpha(props.gridBlockN, 0),
        radiusCurr(props.gridBlockN + 1, 0) {}

int Local::left(const int &index) {
    return index;
}

int Local::right(const int &index) {
    return index + 1;
}

double Local::calcDelRadius(const double &radius, const double &effRadius,
                            const int &gridBlockN) {

    return (effRadius - radius) / gridBlockN;
}

void Local::calcRadiusCurr(const double &radius, const double &effRadius,
                           const int &gridBlockN) {

    dRadius = calcDelRadius(radius, effRadius, gridBlockN);

    for (int i = 0; i < gridBlockN + 1; i++)
        radiusCurr[i] = radius + i * dRadius;
}

void Local::calculateAlpha(const double &dt,
                           const double &radius,
                           const double &effRadius,
                           const double &thrLength) {

    calcRadiusCurr(radius, effRadius, props.gridBlockN);

    for (int i = 0; i < alpha.size(); i++)
        alpha[i] = (M_PI * thrLength *
                    (radiusCurr[right(i)] * radiusCurr[right(i)] -
                     radiusCurr[left(i)] * radiusCurr[left(i)])) / dt;
}
const std::vector<double> Local::getAlpha() const {
    return alpha;
}

const std::vector<double> Local::getRadCurr() const {
    return radiusCurr;
}











