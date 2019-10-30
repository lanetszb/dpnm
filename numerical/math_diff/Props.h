#ifndef PNFLOW_PROPS_H
#define PNFLOW_PROPS_H

#include <vector>
#include <iostream>

class Props {

public:

    explicit Props(const std::vector<double> &_propsVector, const
    std::vector<double> _langmuirCoeff);

    virtual ~Props() {}

    friend std::ostream &operator<<(std::ostream &stream, const Props &props);

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

    std::vector<double> getLangmuirCoeff() const;

    std::vector<double> getPropsVector() const;

    void printPropsVector();

    void printLangmuirCoeff();

    void calcRadiusCurr();

    double time;
    double timeStep;
    double length;
    double radius;
    double effRadius;
    int gridBlockN;
    double concIni;
    double diffusivity;
    double iterativeAccuracy;

    std::vector<double> langmuirCoeff;


private:

    std::vector<double> propsVector;
public:


};


#endif
