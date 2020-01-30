#include <PropsDiffusion.h>
#include <vector>


PropsDiffusion::PropsDiffusion(const std::vector<double> &_propsVector, const
std::vector<double> _langmuirCoeff) :

        time(_propsVector[0]),
        timeStep(_propsVector[1]),
        length(_propsVector[2]),
        radius(_propsVector[3]),
        effRadius(_propsVector[3] * 5),
        gridBlockN(_propsVector[4]),
        concIni(_propsVector[5]),
        diffusivity(_propsVector[6]),
        iterativeAccuracy(_propsVector[7]),
        propsVector(_propsVector),
        langmuirCoeff(_langmuirCoeff) {}


std::ostream &
operator<<(std::ostream &stream, const PropsDiffusion &propsDiffusion) {
    stream << "time " << propsDiffusion.time << std::endl;
    stream << "timeStep " << propsDiffusion.timeStep << std::endl;
    stream << "length " << propsDiffusion.length << std::endl;
    stream << "radius " << propsDiffusion.radius << std::endl;
    stream << "effRadius " << propsDiffusion.effRadius << std::endl;
    stream << "gridBlockN " << propsDiffusion.gridBlockN << std::endl;
    stream << "concIni " << propsDiffusion.length << std::endl;
    stream << "diffusivity " << propsDiffusion.diffusivity << std::endl;
    stream << "iterativeAccuracy " << propsDiffusion.iterativeAccuracy
           << std::endl;

    return stream;
}

std::vector<double> PropsDiffusion::getPropsVector() const {
    return propsVector;
}

std::vector<double> PropsDiffusion::getLangmuirCoeff() const {
    return langmuirCoeff;
}

double PropsDiffusion::getLength() const {
    return length;
}

void PropsDiffusion::setLength(double length) {
    PropsDiffusion::length = length;
}

double PropsDiffusion::getRadius() const {
    return radius;
}

void PropsDiffusion::setRadius(double radius) {
    PropsDiffusion::radius = radius;
}

int PropsDiffusion::getGridBlockN() const {
    return gridBlockN;
}

void PropsDiffusion::setGridBlockN(int gridBlockN) {
    PropsDiffusion::gridBlockN = gridBlockN;
}

double PropsDiffusion::getConcentration() const {
    return concIni;
}

void PropsDiffusion::setConcentration(double concentration) {
    PropsDiffusion::concIni = concentration;
}

double PropsDiffusion::getDiffusivity() const {
    return diffusivity;
}

void PropsDiffusion::setDiffusivity(double diffusivity) {
    PropsDiffusion::diffusivity = diffusivity;
}

double PropsDiffusion::getIterativeAccuracy() const {
    return iterativeAccuracy;
}

void PropsDiffusion::setIterativeAccuracy(double iterativeAccuracy) {
    PropsDiffusion::iterativeAccuracy = iterativeAccuracy;
}

void PropsDiffusion::printPropsVector() {
    for (auto &element : propsVector)
        std::cout << element << std::endl;
}

void PropsDiffusion::printLangmuirCoeff() {
    for (auto &element : langmuirCoeff)
        std::cout << element << std::endl;
}

double PropsDiffusion::getEffRadius() const {
    return effRadius;
}

void PropsDiffusion::setEffRadius(double effRadius) {
    PropsDiffusion::effRadius = effRadius;
}

double PropsDiffusion::getTime() const {
    return time;
}

void PropsDiffusion::setTime(double time) {
    PropsDiffusion::time = time;
}

double PropsDiffusion::getTimeStep() const {
    return timeStep;
}

void PropsDiffusion::setTimeStep(double timeStep) {
    PropsDiffusion::timeStep = timeStep;
}
