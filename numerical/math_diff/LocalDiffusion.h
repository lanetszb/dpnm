#ifndef PNFLOW_LOCALDIFFUSION_H
#define PNFLOW_LOCALDIFFUSION_H

#include <PropsDiffusion.h>

class LocalDiffusion {

public:

    explicit LocalDiffusion(const std::vector<double> &propsVector);

    virtual ~LocalDiffusion() = default;

    static int left(const int &index);

    static int right(const int &index);

    double calcDelRadius(const double &radius, const double &effRadius);

    void calcVolCartesian(const double &frac_height, const double &matrix_width,
                          const double &frac_length, const double &frac_width);

    void calcVolCylindr(const double &radius, const double &effRadius,
                        const double &thrLength);

    void calcRadiusCurr(const double &radius, const double &effRadius);

    void calculateAlpha(const double &dt,
                        const std::vector<double> &vol);

    PropsDiffusion propsDiffusion;

    double dRadius;

    std::vector<double> radiusCurr;
    std::vector<double> volCylindr;
    std::vector<double> volCartes;
    std::vector<double> alpha;

    const std::vector<double> getVolCylindr() const;

    const std::vector<double> getVolCartes() const;

    const std::vector<double> getRadCurr() const;

    const std::vector<double> getAlpha() const;

    // ToDo: const int getGridBlockN() const;


};

#endif
