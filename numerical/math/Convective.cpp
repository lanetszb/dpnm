#include <Convective.h>
#include <cmath>

Convective::Convective(const std::vector<double> &propsVector) :
        Local(propsVector),
        props(propsVector),
        beta(std::vector<double>(props.gridBlockN + 1, 0)) {}

int Convective::left(const int &index) {
    return index - 1;
}

int Convective::right(const int &index) {
    return index;
}

double Convective::omegaCylindric(const double &radius, const double &length) {
    return 2 * M_PI * length * radius;
}

double Convective::omegaSpheric(const double &radius) {
    return 4 * M_PI * radius * radius;
}

std::vector<double> Convective::calc_diffusivityList(const int &gridBlockN,
                                                     const double &diffusivity) {

    return std::vector<double>(gridBlockN + 1, diffusivity);
}

void Convective::calculateBeta(const double &radius,
                               const double &effRadius,
                               const double &length,
                               const double &diffusivity,
                               const int &gridBlockN) {

    dRadius = calcDelRadius(radius, effRadius, gridBlockN);
    auto diffusivityList = calc_diffusivityList(gridBlockN, diffusivity);

    for (int i = 0; i < beta.size(); i++) {

        auto radius = effRadius - i * dRadius;
        auto omega = omegaCylindric(radius, length);

        beta[i] = diffusivityList[i] * omega / dRadius;
    }
}

const std::vector<double> Convective::getBeta() const {
    return beta;
}


