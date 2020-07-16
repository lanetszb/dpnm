#include <PropsPnm.h>
#include <vector>

PropsPnm::PropsPnm(const std::vector<double> &_propsVector) :

        aGasDens(_propsVector[0]),
        bGasDens(_propsVector[1]),
        gasVisc(_propsVector[2]),
        liqDens(_propsVector[3]),
        liqVisc(_propsVector[4]),
        pressIn(_propsVector[5]),
        pressOut(_propsVector[6]),
        itAccuracy(_propsVector[7]),
        propsVector(_propsVector) {}

std::ostream &operator<<(std::ostream &stream, const PropsPnm &propsPNM) {
    stream << "aGasDens " << propsPNM.aGasDens << std::endl;
    stream << "bGasDens " << propsPNM.bGasDens << std::endl;
    stream << "gasVisc " << propsPNM.gasVisc << std::endl;
    stream << "liqDens " << propsPNM.liqDens << std::endl;
    stream << "liqVisc " << propsPNM.liqVisc << std::endl;
    stream << "pressIn " << propsPNM.pressIn << std::endl;
    stream << "pressOut " << propsPNM.pressOut << std::endl;
    stream << "itAccuracy " << propsPNM.itAccuracy << std::endl;

    return stream;
}

std::vector<double> PropsPnm::getPropsVector() const {
    return propsVector;
}

void PropsPnm::printPropsVector() {
    for (auto &element : propsVector)
        std::cout << element << std::endl;
}

double PropsPnm::getAGasDens() const {
    return aGasDens;
}

void PropsPnm::setAGasDens(double aGasDens) {
    PropsPnm::aGasDens = aGasDens;
}

double PropsPnm::getBGasDens() const {
    return bGasDens;
}

void PropsPnm::setBGasDens(double bGasDens) {
    PropsPnm::bGasDens = bGasDens;
}

double PropsPnm::getGasVisc() const {
    return gasVisc;
}

void PropsPnm::setGasVisc(double gasVisc) {
    PropsPnm::gasVisc = gasVisc;
}

double PropsPnm::getLiqDens() const {
    return liqDens;
}

void PropsPnm::setLiqDens(double liqDens) {
    PropsPnm::liqDens = liqDens;
}

double PropsPnm::getLiqVisc() const {
    return liqVisc;
}

void PropsPnm::setLiqVisc(double liqVisc) {
    PropsPnm::liqVisc = liqVisc;
}







