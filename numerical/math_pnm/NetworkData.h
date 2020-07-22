#ifndef PNFLOW_NETWORKDATA_H
#define PNFLOW_NETWORKDATA_H

#include <vector>
#include <iostream>


class NetworkData {

public :

    explicit NetworkData(const std::vector<int> &_fracturesList,
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
                         const std::vector<double> &_hydraulicCond);

    virtual ~NetworkData() {}

    void findBoundaryPores(std::vector<double> &poreCoord);

    void calcThroatConns();

    void calcPor2ThrConns();

    void calcBoundPoresSizes();

    std::vector<int> fracturesList;
    std::vector<double> fracturesHeights;
    std::vector<double> fracturesLengths;
    std::vector<double> fracturesWidths;
    int fracturesN;

    std::vector<double> fracsConnIndIn;
    std::vector<double> fracsConnIndOut;


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
