#include <PropsPnm.h>
#include <Eigen/Sparse>
#include <cmath>
#include <string>
#include <vector>
#include <EquationPNM.h>
#include <NetworkData.h>

NetworkData::NetworkData(const std::vector<int> &_fracturesList,
                         const std::vector<double> &_fracturesHeights,
                         const std::vector<double> &_fracturesLengths,
                         const std::vector<double> &_fracturesWidths,
                         const std::vector<double> &_fracsConnIndIn,
                         const std::vector<double> &_fracsConnIndOut,
                         const std::vector<double> &_poresCoordsX,
                         const std::vector<double> &_poresCoordsY,
                         const std::vector<double> &_poresCoordsZ,
                         const std::vector<double> &_poresRadii,
                         const std::vector<int> &_poresList,
                         const std::vector<bool> &_poresInlet,
                         const std::vector<bool> &_poresOutlet,
                         const std::vector<double> &_hydraulicCond) :

        fracturesList(_fracturesList),
        fracturesHeights(_fracturesHeights),
        fracturesLengths(_fracturesLengths),
        fracturesWidths(_fracturesWidths),
        fracsConnIndIn(_fracsConnIndIn),
        fracsConnIndOut(_fracsConnIndOut),
        fracturesN(fracturesList.size()),
        throatConns(fracturesN, std::make_pair(-1, -1)),

        poresCoordsX(_poresCoordsX),
        poresCoordsY(_poresCoordsY),
        poresCoordsZ(_poresCoordsZ),

        poresRadii(_poresRadii),
        poresList(_poresList),
        poreN(poresList.size()),
        por2thrConns(poreN),

        poreInlet(_poresInlet),
        poreOutlet(_poresOutlet),
        boundPoresLeftSize(0),
        boundPoresRightSize(0),
        boundPoresSize(0),

        hydraulicCond(_hydraulicCond) {}


void NetworkData::calcThroatConns() {

    for (int i = 0; i < fracturesN; i++) {
        throatConns[i].first = fracsConnIndIn[i];
        throatConns[i].second = fracsConnIndOut[i];
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
        if (poreInlet[i])
            boundPoresLeftSize += 1;

    std::cout << boundPoresLeftSize << std::endl;

    for (int i = 0; i < poreN; i++)
        if (poreOutlet[i])
            boundPoresRightSize += 1;

    boundPoresSize = boundPoresLeftSize + boundPoresRightSize;
}