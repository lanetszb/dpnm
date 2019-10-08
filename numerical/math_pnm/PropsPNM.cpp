#include <PropsPNM.h>
#include <vector>

PropsPNM::PropsPNM(const std::vector<double> &_propsVector) :

        aGasDens(_propsVector[0]),
        bGasDens(_propsVector[1]),
        gasVisc(_propsVector[2]),
        liqDens(_propsVector[3]),
        liqVisc(_propsVector[4]),
        propsVector(_propsVector) {}

std::ostream &operator<<(std::ostream &stream, const PropsPNM &propsPNM) {
    stream << "a_gas_dens " << propsPNM.aGasDens << std::endl;
    stream << "b_gas_dens " << propsPNM.bGasDens << std::endl;
    stream << "gas_visc " << propsPNM.gasVisc << std::endl;
    stream << "liq_dens " << propsPNM.liqDens << std::endl;
    stream << "liq_visc " << propsPNM.liqVisc << std::endl;

    return stream;
}

std::vector<double> PropsPNM::getPropsVector() const {
    return propsVector;
}

void PropsPNM::printPropsVector() {
    for (auto &element : propsVector)
        std::cout << element << std::endl;
}

double PropsPNM::getAGasDens() const {
    return aGasDens;
}

void PropsPNM::setAGasDens(double aGasDens) {
    PropsPNM::aGasDens = aGasDens;
}

double PropsPNM::getBGasDens() const {
    return bGasDens;
}

void PropsPNM::setBGasDens(double bGasDens) {
    PropsPNM::bGasDens = bGasDens;
}

double PropsPNM::getGasVisc() const {
    return gasVisc;
}

void PropsPNM::setGasVisc(double gasVisc) {
    PropsPNM::gasVisc = gasVisc;
}

double PropsPNM::getLiqDens() const {
    return liqDens;
}

void PropsPNM::setLiqDens(double liqDens) {
    PropsPNM::liqDens = liqDens;
}

double PropsPNM::getLiqVisc() const {
    return liqVisc;
}

void PropsPNM::setLiqVisc(double liqVisc) {
    PropsPNM::liqVisc = liqVisc;
}







