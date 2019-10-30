#ifndef PNFLOW_LOCAL_H
#define PNFLOW_LOCAL_H

#include <Props.h>

class Local {

public:

    explicit Local(const std::vector<double> &propsVector,
                   const std::vector<double> &langmuirCoeff);

    virtual ~Local() = default;

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

    const std::vector<double> getAlpha() const;

    const std::vector<double> getRadCurr() const;

    Props props;

    double dRadius;
    std::vector<double> radiusCurr;

    std::vector<double> alpha;


};


#endif
