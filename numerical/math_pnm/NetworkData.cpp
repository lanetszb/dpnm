#include <string>
#include <vector>
#include <map>
#include <NetworkData.h>

NetworkData::NetworkData(
        const std::map<std::string, std::variant<std::vector<bool>,
                std::vector<int>, std::vector<double>>> &paramsNetwork) :


        _paramsNetwork(paramsNetwork),
        fracturesList(
                std::get<std::vector<int>>(_paramsNetwork["fracList"])),
        fracturesHeights(
                std::get<std::vector<double>>(_paramsNetwork["fracHeights"])),
        fracturesLengths(
                std::get<std::vector<double>>(_paramsNetwork["fracLengths"])),
        fracturesWidths(
                std::get<std::vector<double>>(_paramsNetwork["fracWidths"])),
        fracsConnIndIn(
                std::get<std::vector<int>>(_paramsNetwork["fracConnIndIn"])),
        fracsConnIndOut(std::get<std::vector<int>>(
                _paramsNetwork["fracConnIndOut"])),
        fracturesN(fracturesList.size()),
        throatConns(fracturesN, std::make_pair(-1, -1)),
        poresCoordsX(
                std::get<std::vector<double>>(_paramsNetwork["poresCoordsX"])),
        poresCoordsY(
                std::get<std::vector<double>>(_paramsNetwork["poresCoordsY"])),
        poresCoordsZ(
                std::get<std::vector<double>>(_paramsNetwork["poresCoordsZ"])),
        poresRadii(std::get<std::vector<double>>(_paramsNetwork["poresRadii"])),
        poresList(std::get<std::vector<int>>(_paramsNetwork["poreList"])),
        poreN(poresList.size()),
        por2thrConns(poreN),
        poreInlet(std::get<std::vector<bool>>(_paramsNetwork["poreInlet"])),
        poreOutlet(std::get<std::vector<bool>>(_paramsNetwork["poreOutlet"])),
        boundPoresLeftSize(0),
        boundPoresRightSize(0),
        boundPoresSize(0),

        hydraulicCond(std::get<std::vector<double>>(
                _paramsNetwork["hydraulicCond"])) {}


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