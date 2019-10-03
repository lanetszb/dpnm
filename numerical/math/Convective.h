#ifndef PNFLOW_CONVECTIVE_H
#define PNFLOW_CONVECTIVE_H

#include <Props.h>
#include <Local.h>

class Convective : public Local {

public:

    explicit Convective(const std::vector<double> &propsVector);

    virtual ~Convective() = default;


    void calculateBeta(const double &radius,
                       const double &effRadius,
                       const double &length,
                       const double &diffusivity,
                       const int &gridBlockN);

    double omegaCylindric(const double &radius, const double &length);

    double omegaSpheric(const double &radius);

    std::vector<double> calc_diffusivityList(const int &gridBlockN,
                                             const double &diffusivity);


    Props props;

    std::vector<double> beta;

    const std::vector<double> getBeta() const;

};


#endif
