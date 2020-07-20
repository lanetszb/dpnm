#include <ConvectiveDiffusion.h>
#include <cmath>

ConvectiveDiffusion::ConvectiveDiffusion(const std::vector<double> &propsVector)
        :
        propsDiffusion(propsVector),
        localDiffusion(propsVector),
        beta(propsDiffusion.gridBlockN + 1, 0),
        omegaCylindric(propsDiffusion.gridBlockN + 1, 0),
        omegaCartesian(propsDiffusion.gridBlockN + 1, 0) {}

void ConvectiveDiffusion::calcOmegaCylindr(const double &length,
                                           const double &radius,
                                           const double &effRadius) {

    localDiffusion.calcMatrCoordCurr(radius, effRadius);

    for (int i = 0; i < beta.size(); i++)
        omegaCylindric[i] = 2 * M_PI * localDiffusion.radiusCurr[i] * length;
}

void ConvectiveDiffusion::calcOmegaCartes(const double &frac_height,
                                          const double &frac_length) {

    for (int i = 0; i < beta.size(); i++)
        omegaCartesian[i] = frac_height * frac_length;
}

double ConvectiveDiffusion::omegaSpheric(const double &radius) {
    return 4 * M_PI * radius * radius;
}

std::vector<double>
ConvectiveDiffusion::calc_diffusivityList(const int &gridBlockN,
                                          const double &diffusivity) {

    return std::vector<double>(gridBlockN, diffusivity);
}

void ConvectiveDiffusion::calculateBeta(const double &radius,
                                        const double &effRadius,
                                        const double &length,
                                        const double &diffusivity,
                                        const int &gridBlockN,
                                        const std::vector<double> &omega) {

    auto dRadius = localDiffusion.calcDelRadius(radius, effRadius);
    auto diffusivityList = calc_diffusivityList(beta.size(),
                                                diffusivity);

    for (int i = 0; i < beta.size(); i++)
        beta[i] = diffusivityList[i] * omega[i] / dRadius;
}

const std::vector<double> ConvectiveDiffusion::getOmegaCylindr() const {
    return omegaCylindric;
}

const std::vector<double> ConvectiveDiffusion::getOmegaCartes() const {
    return omegaCartesian;
}


