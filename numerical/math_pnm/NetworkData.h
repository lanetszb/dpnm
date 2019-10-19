#ifndef PNFLOW_NETWORKDATA_H
#define PNFLOW_NETWORKDATA_H

#include <vector>
#include <iostream>


class NetworkData {

public :

    explicit NetworkData(const std::vector<int> & _throat_list,
                         const std::vector<double> &_throat_radius,
                         const std::vector<double> &_throat_length,
                         const std::vector<double> &_conn_ind_in,
                         const std::vector<double> &_conn_ind_out,

                         const std::vector<double> &_pore_coord_x,
                         const std::vector<double> &_pore_coord_y,
                         const std::vector<double> &_pore_coord_z,
                         const std::vector<double> &_pore_radius,
                         const std::vector<int> &_pore_list,
                         const std::vector<int> &_pore_conns,
                         const std::vector<int> &_conn_number,
                         const std::vector<int> &_pore_per_row
    );

    virtual ~NetworkData() {}

    friend std::ostream &operator<<(
            std::ostream &stream, const NetworkData &networkData);

    void printThroatRadius();

    void printThroatLength();

    void printPoreList();

    void findBoundaryPores(std::vector<double> &poreCoord);

    std::vector<double> getThroatRadius() const;

    std::vector<double> getThroatLength() const;

    std::vector<int> getPoreList() const;


    std::vector<double> throatRadius;
    std::vector<double> throatLength;
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

    std::vector<int> boundaryPores;
    std::vector<int> boundaryPoresIn;
    std::vector<int> boundaryPoresOut;


private:

};


#endif
