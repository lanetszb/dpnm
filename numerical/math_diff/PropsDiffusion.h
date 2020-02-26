#ifndef PNFLOW_PROPSDIFFUSION_H
#define PNFLOW_PROPSDIFFUSION_H

#include <vector>
#include <iostream>

class PropsDiffusion {

public:

    explicit PropsDiffusion(const std::vector<double> &_propsVector);

    virtual ~PropsDiffusion() {}

    friend std::ostream &operator<<(std::ostream &stream, const PropsDiffusion &props);

    double getTime() const;

    void setTime(double time);

    double getTimeStep() const;

    void setTimeStep(double timeStep);

    double getLength() const;

    void setLength(double length);

    double getRadius() const;

    void setRadius(double radius);

    double getEffRadius() const;

    void setEffRadius(double effRadius);

    int getGridBlockN() const;

    void setGridBlockN(int gridBlockN);

    double getConcentration() const;

    void setConcentration(double concentration);

    double getDiffusivity() const;

    void setDiffusivity(double diffusivity);

    double getIterativeAccuracy() const;

    void setIterativeAccuracy(double iterativeAccuracy);

    double getGridBlockSize() const;

    std::vector<double> getLangmuirCoeff() const;

    std::vector<double> getPropsVector() const;

    void printPropsVector();

    void printLangmuirCoeff();

    double time;
    double timeStep;
    double length;
    double radius;
    double effRadius;
    int gridBlockN;
    double concIni;
    double diffusivity;
    double iterativeAccuracy;

private:

    std::vector<double> propsVector;
public:


};


#endif
