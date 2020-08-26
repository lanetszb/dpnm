#ifndef PNFLOW_LOCALDIFFUSION_H
#define PNFLOW_LOCALDIFFUSION_H

#include <PropsDiffusion.h>

class LocalDiffusion {

public:

    explicit LocalDiffusion(PropsDiffusion &propsDiffusion);

    virtual ~LocalDiffusion() = default;

    static int left(const int &index);

    static int right(const int &index);

    double calcDelRadius(const double &radius, const double &effRadius);

    void calcVol(const std::string &coordType);

    void calcVolCartesian(const double &frac_height, const double &matrix_width,
                          const double &frac_length, const double &frac_width);

    void calcVolCylindr(const double &radius, const double &effRadius,
                        const double &thrLength);

    void calcMatrCoordCurr(const double &radius, const double &effRadius);

    void calculateAlpha(const double &dt,
                        const std::vector<double> &vol);

    PropsDiffusion &propsDiffusion;

    double dRadius;

    std::vector<double> radiusCurr;
    std::vector<double> volCylindr;
    std::vector<double> volCartes;
    std::vector<double> alpha;



    // ToDo: const int getGridBlockN() const;


};

#endif
