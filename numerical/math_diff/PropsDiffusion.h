#ifndef PNFLOW_PROPSDIFFUSION_H
#define PNFLOW_PROPSDIFFUSION_H

#include <vector>
#include <iostream>
#include <map>
#include <variant>

class PropsDiffusion {

public:

    explicit PropsDiffusion(
            const std::map<std::string, std::variant<int, double>> &params);

    virtual ~PropsDiffusion() {}

    void printParams();

    std::map<std::string, std::variant<int, double>> _params;

    double time;
    double timeStep;
    int gridBlockN;
    double length;
    double radius;
    double effRadius;
    double concIni;
    double diffusivity;
    double iterativeAccuracy;

private:

};


#endif
