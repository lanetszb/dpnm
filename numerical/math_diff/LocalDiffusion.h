#ifndef PNFLOW_LOCALDIFFUSION_H
#define PNFLOW_LOCALDIFFUSION_H

#include <PropsDiffusion.h>

class LocalDiffusion {

public:

    explicit LocalDiffusion(const std::vector<double> &propsVector,
                            const std::vector<double> &langmuirCoeff);

    virtual ~LocalDiffusion() = default;

    static int left(const int &index);

    static int right(const int &index);

    double calcDelRadius(const double &radius, const double &effRadius,
                         const int &gridBlockN);

    void calculateAlpha(const double &dt,
                        const double &radius,
                        const double &effRadius,
                        const double &thrLength);

    void calcRadiusCurr(const double &radius, const double &effRadius,
                        const int &gridBlockN);

    PropsDiffusion propsDiffusion;

    double dRadius;

    std::vector<double> radiusCurr;

    std::vector<double> alpha;

    const std::vector<double> getAlpha() const;

    const std::vector<double> getRadCurr() const;

};

#endif
