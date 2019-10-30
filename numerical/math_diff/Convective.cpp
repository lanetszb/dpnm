#include <Convective.h>
#include <cmath>

Convective::Convective(const std::vector<double> &propsVector,
                       const std::vector<double> &langmuirCoeff) :
        props(propsVector, langmuirCoeff),
        local(propsVector, langmuirCoeff),
        beta(props.gridBlockN + 1, 0) {}

int Convective::left(const int &index) {
    return index - 1;
}

int Convective::right(const int &index) {
    return index;
}

double Convective::omegaCylindric(const double &radius, const double &length) {
    return 2 * M_PI * radius * length;
}

double Convective::omegaSpheric(const double &radius) {
    return 4 * M_PI * radius * radius;
}

std::vector<double> Convective::calc_diffusivityList(const int &gridBlockN,
                                                     const double &diffusivity) {

    return std::vector<double>(gridBlockN, diffusivity);
}

void Convective::calculateBeta(const double &radius,
                               const double &effRadius,
                               const double &length,
                               const double &diffusivity,
                               const double &gridBlockN) {

    local.calcRadiusCurr(radius, effRadius, gridBlockN);
    auto diffusivityList = calc_diffusivityList(beta.size(), diffusivity);

    for (int i = 0; i < beta.size(); i++) {
        auto omega = omegaCylindric(local.radiusCurr[i], length);
        beta[i] = diffusivityList[i] * omega / local.dRadius;
    }
}

const std::vector<double> Convective::getBeta() const {
    return beta;
}


