#ifndef PNFLOW_NETWORKDATA_H
#define PNFLOW_NETWORKDATA_H

#include <vector>
#include <iostream>


class NetworkData {

public :

    explicit NetworkData(const std::vector<int> &_throatList,
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
                         const std::vector<int> &_poreConns,
                         const std::vector<int> &_connNumber,
                         const std::vector<int> &_porePerRow,
                         const std::vector<bool> &_poreLeftX,
                         const std::vector<bool> &_poreRightX,
                         const std::vector<double> &_hydraulicCond);

    virtual ~NetworkData() {}

    void printThroatRadius();

    void printThroatLength();

    void printPoreList();

    void findBoundaryPores(std::vector<double> &poreCoord);

    void calcThroatConns();

    void calcPorConns();

    void calcBoundPoresSizes();

    std::vector<double> getThroatRadius() const;

    std::vector<double> getThroatLength() const;

    std::vector<int> getPoreList() const;

    std::vector<double> throatRadius;
    std::vector<double> throatLength;
    std::vector<double> throatWidth;
    std::vector<int> throatList;
    int throatN;

    std::vector<double> connIndIn;
    std::vector<double> connIndOut;


    std::vector<double> poreCoordX;
    std::vector<double> poreCoordY;
    std::vector<double> poreCoordZ;

    std::vector<double> poreRadius;

    std::vector<int> poreList;
    int poreN;

    std::vector<int> poreConns;
    std::vector<int> connNumber;
    std::vector<int> porPerRow;

    std::vector<bool> poreLeftX;
    std::vector<bool> poreRightX;

    std::vector<int> boundaryPores;
    std::vector<int> boundaryPoresIn;
    std::vector<int> boundaryPoresOut;

    std::vector<double> hydraulicCond;

    std::vector<std::pair<int, int>> throatConns;
    std::vector<std::vector<int>> porConns;

    int boundPoresLeftSize;
    int boundPoresRightSize;
    int boundPoresSize;


private:

};


#endif
