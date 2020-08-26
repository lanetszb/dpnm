#ifndef PNFLOW_NETWORKDATA_H
#define PNFLOW_NETWORKDATA_H

#include <vector>
#include <iostream>
#include <map>
#include <variant>


class NetworkData {

public :

    explicit NetworkData(const std::map<std::string, std::variant<std::vector<bool>,
            std::vector<int>, std::vector<double>>> &_paramsNetwork);

    virtual ~NetworkData() {}

    void findBoundaryPores(std::vector<double> &poreCoord);

    void calcThroatConns();

    void calcPor2ThrConns();

    void calcBoundPoresSizes();

    std::map<std::string, std::variant<std::vector<bool>, std::vector<int>,
            std::vector<double>>> _paramsNetwork;

    std::vector<int> fracturesList;
    std::vector<double> fracturesHeights;
    std::vector<double> fracturesLengths;
    std::vector<double> fracturesWidths;
    int fracturesN;

    std::vector<int> fracsConnIndIn;
    std::vector<int> fracsConnIndOut;


    std::vector<double> poresCoordsX;
    std::vector<double> poresCoordsY;
    std::vector<double> poresCoordsZ;

    std::vector<double> poresRadii;

    std::vector<int> poresList;
    int poreN;

    std::vector<bool> poreInlet;
    std::vector<bool> poreOutlet;

    std::vector<int> boundaryPores;
    std::vector<int> boundaryPoresIn;
    std::vector<int> boundaryPoresOut;

    std::vector<double> hydraulicCond;

    std::vector<std::pair<int, int>> throatConns;
    std::vector<std::vector<int>> por2thrConns;

    int boundPoresLeftSize;
    int boundPoresRightSize;
    int boundPoresSize;


private:

};


#endif
