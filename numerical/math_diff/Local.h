#ifndef PNFLOW_LOCAL_H
#define PNFLOW_LOCAL_H

#include <Props.h>

class Local {

public:

    explicit Local(const std::vector<double> &propsVector);

    virtual ~Local() = default;

    static int left(const int &index);

    static int right(const int &index);

    double calcDelRadius(const double &radius, const double &effRadius,
                         const int &gridBlockN);

    void calculateAlpha(const double &dt, const double &radius,
                        const double &effRadius);

    Props props;

    double dRadius;

    std::vector<double> alpha;

    const std::vector<double> getAlpha() const;

};


#endif
