#ifndef PNFLOW_CONVECTIVEDIFFUSION_H
#define PNFLOW_CONVECTIVEDIFFUSION_H

#include <PropsDiffusion.h>
#include <LocalDiffusion.h>

class ConvectiveDiffusion {

public:

    explicit ConvectiveDiffusion(const std::vector<double> &propsVector);

    virtual ~ConvectiveDiffusion() = default;

    void calculateBeta(const double &radius,
                       const double &effRadius,
                       const double &length,
                       const double &diffusivity,
                       const int &gridBlockN,
                       const std::vector<double> &omega);

    void calcOmegaCylindr(const double &length,
                          const double &radius,
                          const double &effRadius);

    void calcOmegaCartes(const double &frac_height,
                       const double &frac_length);

    double omegaSpheric(const double &radius);

    std::vector<double> calc_diffusivityList(const int &gridBlockN,
                                             const double &diffusivity);

    PropsDiffusion propsDiffusion;
    LocalDiffusion localDiffusion;

    std::vector<double> omegaCylindric;
    std::vector<double> omegaCartesian;

    std::vector<double> beta;

    const std::vector<double> getOmegaCylindr() const;

    const std::vector<double> getOmegaCartes() const;

};

#endif
