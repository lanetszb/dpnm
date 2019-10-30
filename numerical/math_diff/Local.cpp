#include <Local.h>
#include <cmath>

Local::Local(const std::vector<double> &propsVector,
             const std::vector<double> &langmuirCoeff) :
        props(propsVector, langmuirCoeff),
        dRadius(0),
        alpha(std::vector<double>(props.gridBlockN, 0)) {}

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

void Local::calculateAlpha(const double &dt,
                           const double &radius,
                           const double &effRadius,
                           const double &thrLength) {

    dRadius = calcDelRadius(radius, effRadius, alpha.size());

    for (int i = 0; i < alpha.size(); i++) {

        auto r_out = radius + (i + 1) * dRadius;
        auto r_in = radius + i * dRadius;

        alpha[i] = (M_PI * thrLength * (r_out * r_out - r_in * r_in)) / dt;
    }
}

const std::vector<double> Local::getAlpha() const {
    return alpha;
}











