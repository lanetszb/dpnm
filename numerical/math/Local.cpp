#include <Local.h>

#include <fstream>
#include <cmath>

Local::Local(const std::vector<double> &propsVector) :
        props(propsVector),
        dRadius(0),
        alpha(std::vector<double>(props.gridBlockN, 0)) {}

int Local::left(const int &index) {
    return index;
}

int Local::right(const int &index) {
    return index + 1;
}


std::vector<double> Local::calc_concListIni(const int &gridBlockN,
                                            const double &concentration) {

    return std::vector<double>(gridBlockN, concentration);
}

double Local::calcDelRadius(const double &radius, const double &effRadius,
                            const int &gridBlockN) {

    return (effRadius - radius) / gridBlockN;
}


void Local::calculateAlpha(const double &dt, const double &radius,
                           const double &effRadius, const int &gridBlockN) {

    dRadius = calcDelRadius(radius, effRadius, gridBlockN);

    for (int i = 0; i < alpha.size(); i++) {

        auto r_out = effRadius - i * dRadius;
        auto r_in = effRadius - (i + 1) * dRadius;

        alpha[i] = (M_PI * props.length * (r_out * r_out - r_in * r_in)) / dt;
    }


}

const std::vector<double> Local::getAlpha() const {
    return alpha;
}











