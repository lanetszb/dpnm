#ifndef PNFLOW_PROPSPNM_H
#define PNFLOW_PROPSPNM_H

#include <vector>
#include <iostream>

class PropsPnm {

public:

    explicit PropsPnm(const std::vector<double> &_propsVector);

    virtual ~PropsPnm() {}
    
    friend std::ostream & operator<<(std::ostream &stream,
        const PropsPnm &propsPNM);

    double getAGasDens() const;

    void setAGasDens(double aGasDens);

    double getBGasDens() const;

    void setBGasDens(double bGasDens);

    double getGasVisc() const;

    void setGasVisc(double gasVisc);

    double getLiqDens() const;

    void setLiqDens(double liqDens);

    double getLiqVisc() const;

    void setLiqVisc(double liqVisc);

    void printPropsVector();


    double aGasDens;
    double bGasDens;
    double gasVisc;
    double liqDens;
    double liqVisc;
    double pressIn;
    double pressOut;
    double itAccuracy;

    std::vector<double> getPropsVector() const;

private:
    // ToDo make it map and use it properly
    std::vector<double> propsVector;
public:


};


#endif