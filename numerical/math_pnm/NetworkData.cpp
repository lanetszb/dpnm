#include <PropsPnm.h>
#include <Eigen/Sparse>
#include <cmath>
#include <string>
#include <vector>
#include <EquationPNM.h>
#include <NetworkData.h>

NetworkData::NetworkData(const std::vector<int> &_throatList,
                         const std::vector<double> &_throatHeight,
                         const std::vector<double> &_throatLength,
                         const std::vector<double> &_throatWidth,
                         const std::vector<double> &_connIndIn,
                         const std::vector<double> &_connIndOut,
                         const std::vector<double> &_poreCoordX,
                         const std::vector<double> &_poreCoordY,
                         const std::vector<double> &_poreCoordZ,
                         const std::vector<double> &_poreRadius,
                         const std::vector<int> &_poreList,
                         const std::vector<int> &_pore2poreConns,
                         const std::vector<int> &_connNumber,
                         const std::vector<int> &_porePerRow,
                         const std::vector<bool> &_poreLeftX,
                         const std::vector<bool> &_poreRightX,
                         const std::vector<double> &_hydraulicCond) :

        throatList(_throatList),
        throatRadius(_throatHeight),
        throatLength(_throatLength),
        throatWidth(_throatWidth),
        connIndIn(_connIndIn),
        connIndOut(_connIndOut),
        throatN(throatList.size()),
        throatConns(throatN, std::make_pair(-1, -1)),

        poreCoordX(_poreCoordX),
        poreCoordY(_poreCoordY),
        poreCoordZ(_poreCoordZ),

        poreRadius(_poreRadius),
        poreList(_poreList),
        poreN(poreList.size()),

        pore2poreConns(_pore2poreConns),
        por2thrConns(poreN),
        connNumber(_connNumber),
        porPerRow(_porePerRow),

        poreLeftX(_poreLeftX),
        poreRightX(_poreRightX),
        boundPoresLeftSize(0),
        boundPoresRightSize(0),
        boundPoresSize(0),

        hydraulicCond(_hydraulicCond) {}

std::vector<double> NetworkData::getThroatRadius() const {
    return throatRadius;
}

void NetworkData::printThroatRadius() {
    for (auto &element : throatRadius)
        std::cout << element << std::endl;
}

std::vector<double> NetworkData::getThroatLength() const {
    return throatLength;
}

void NetworkData::printThroatLength() {
    for (auto &element : throatLength)
        std::cout << element << std::endl;
}

std::vector<int> NetworkData::getPoreList() const {
    return poreList;
}

void NetworkData::printPoreList() {
    for (auto &element : poreList)
        std::cout << element << std::endl;
}


void NetworkData::calcThroatConns() {

    for (int i = 0; i < throatN; i++) {
        throatConns[i].first = connIndIn[i];
        throatConns[i].second = connIndOut[i];
    }
}

void NetworkData::calcPor2ThrConns() {

    for (int i = 0; i < throatConns.size(); i++) {
        por2thrConns[throatConns[i].first].emplace_back(i);
        por2thrConns[throatConns[i].second].emplace_back(i);
    }
}

void NetworkData::findBoundaryPores(std::vector<double> &poreCoord) {

    auto min = std::min_element(std::begin(poreCoord), std::end(poreCoord));
    auto max = std::max_element(std::begin(poreCoord), std::end(poreCoord));

    for (int i = 0; i < poreCoord.size(); i++) {
        if (poreCoord[i] == *min)
            boundaryPoresIn.emplace_back(i);
        else if (poreCoord[i] == *max)
            boundaryPoresOut.emplace_back(i);
    }

    for (int i = 0; i < boundaryPoresIn.size(); i++)
        boundaryPores.emplace_back(boundaryPoresIn[i]);

    boundaryPores.insert(boundaryPores.end(),
                         boundaryPoresOut.begin(),
                         boundaryPoresOut.end());
}

void NetworkData::calcBoundPoresSizes() {

    for (int i = 0; i < poreN; i++)
        if (poreLeftX[i])
            boundPoresLeftSize += 1;

    std::cout << boundPoresLeftSize << std::endl;

    for (int i = 0; i < poreN; i++)
        if (poreRightX[i])
            boundPoresRightSize += 1;

    boundPoresSize = boundPoresLeftSize + boundPoresRightSize;
}