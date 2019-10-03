#include <Props.h>

#include <vector>


Props::Props(const std::vector<double> &_propsVector) :

        time(_propsVector[0]),
        timeStep(_propsVector[1]),
        length(_propsVector[2]),
        radius(_propsVector[3]),
        effRadius(_propsVector[3] * 5),
        gridBlockN(_propsVector[4]),
        concentration(_propsVector[5]),
        diffusivity(_propsVector[6]),
        iterativeAccuracy(_propsVector[7]),
        propsVector(_propsVector) {}


std::ostream &operator<<(std::ostream &stream, const Props &props) {
    stream << "time " << props.time << std::endl;
    stream << "timeStep " << props.timeStep << std::endl;
    stream << "length " << props.length << std::endl;
    stream << "radius " << props.radius << std::endl;
    stream << "effRadius " << props.effRadius << std::endl;
    stream << "gridBlockN " << props.gridBlockN << std::endl;
    stream << "concentration " << props.length << std::endl;
    stream << "diffusivity " << props.diffusivity << std::endl;
    stream << "iterativeAccuracy " << props.iterativeAccuracy << std::endl;

    return stream;
}

std::vector<double> Props::getPropsVector() const {
    return propsVector;
}

double Props::getLength() const {
    return length;
}

void Props::setLength(double length) {
    Props::length = length;
}

double Props::getRadius() const {
    return radius;
}

void Props::setRadius(double radius) {
    Props::radius = radius;
}

int Props::getGridBlockN() const {
    return gridBlockN;
}

void Props::setGridBlockN(int gridBlockN) {
    gridBlockN = gridBlockN;
}

double Props::getConcentration() const {
    return concentration;
}

void Props::setConcentration(double concentration) {
    Props::concentration = concentration;
}

double Props::getDiffusivity() const {
    return diffusivity;
}

void Props::setDiffusivity(double diffusivity) {
    Props::diffusivity = diffusivity;
}

double Props::getIterativeAccuracy() const {
    return iterativeAccuracy;
}

void Props::setIterativeAccuracy(double iterativeAccuracy) {
    Props::iterativeAccuracy = iterativeAccuracy;
}

void Props::printPropsVector() {
    for (auto &element : propsVector)
        std::cout << element << std::endl;
}

double Props::getEffRadius() const {
    return effRadius;
}

void Props::setEffRadius(double effRadius) {
    Props::effRadius = effRadius;
}

double Props::getTime() const {
    return time;
}

void Props::setTime(double time) {
    Props::time = time;
}

double Props::getTimeStep() const {
    return timeStep;
}

void Props::setTimeStep(double timeStep) {
    Props::timeStep = timeStep;
}
