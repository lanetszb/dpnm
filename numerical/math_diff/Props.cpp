#include <PropsPNM.h>

#include <vector>


PropsPNM::PropsPNM(const std::vector<double> &_propsVector) :

        time(_propsVector[0]),
        timeStep(_propsVector[1]),
        length(_propsVector[2]),
        radius(_propsVector[3]),
        effRadius(_propsVector[3] * 5),
        gridBlockN(_propsVector[4]),
        concIni(_propsVector[5]),
        diffusivity(_propsVector[6]),
        iterativeAccuracy(_propsVector[7]),
        propsVector(_propsVector) {}


std::ostream &operator<<(std::ostream &stream, const PropsPNM &props) {
    stream << "time " << props.time << std::endl;
    stream << "timeStep " << props.timeStep << std::endl;
    stream << "length " << props.length << std::endl;
    stream << "radius " << props.radius << std::endl;
    stream << "effRadius " << props.effRadius << std::endl;
    stream << "gridBlockN " << props.gridBlockN << std::endl;
    stream << "concIni " << props.length << std::endl;
    stream << "diffusivity " << props.diffusivity << std::endl;
    stream << "iterativeAccuracy " << props.iterativeAccuracy << std::endl;

    return stream;
}

std::vector<double> PropsPNM::getPropsVector() const {
    return propsVector;
}

double PropsPNM::getLength() const {
    return length;
}

void PropsPNM::setLength(double length) {
    PropsPNM::length = length;
}

double PropsPNM::getRadius() const {
    return radius;
}

void PropsPNM::setRadius(double radius) {
    PropsPNM::radius = radius;
}

int PropsPNM::getGridBlockN() const {
    return gridBlockN;
}

void PropsPNM::setGridBlockN(int gridBlockN) {
    gridBlockN = gridBlockN;
}

double PropsPNM::getConcentration() const {
    return concIni;
}

void PropsPNM::setConcentration(double concentration) {
    PropsPNM::concIni = concentration;
}

double PropsPNM::getDiffusivity() const {
    return diffusivity;
}

void PropsPNM::setDiffusivity(double diffusivity) {
    PropsPNM::diffusivity = diffusivity;
}

double PropsPNM::getIterativeAccuracy() const {
    return iterativeAccuracy;
}

void PropsPNM::setIterativeAccuracy(double iterativeAccuracy) {
    PropsPNM::iterativeAccuracy = iterativeAccuracy;
}

void PropsPNM::printPropsVector() {
    for (auto &element : propsVector)
        std::cout << element << std::endl;
}

double PropsPNM::getEffRadius() const {
    return effRadius;
}

void PropsPNM::setEffRadius(double effRadius) {
    PropsPNM::effRadius = effRadius;
}

double PropsPNM::getTime() const {
    return time;
}

void PropsPNM::setTime(double time) {
    PropsPNM::time = time;
}

double PropsPNM::getTimeStep() const {
    return timeStep;
}

void PropsPNM::setTimeStep(double timeStep) {
    PropsPNM::timeStep = timeStep;
}
