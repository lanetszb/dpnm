#ifndef PNFLOW_NETWORKDATA_H
#define PNFLOW_NETWORKDATA_H

#include <vector>
#include <iostream>


class NetworkData {

public :

    explicit NetworkData(const std::vector<double> &_throat_radius,
                         const std::vector<double> &_throat_length,
                         const std::vector<double> &_conn_ind_in,
                         const std::vector<double> &_conn_ind_out,
                         const std::vector<double> &_pore_coord_x,
                         const std::vector<double> &_pore_coord_y,
                         const std::vector<double> &_pore_coord_z,
                         const std::vector<double> &_pore_radius,
                         const std::vector<int> &_pore_list);

    virtual ~NetworkData() {}

    friend std::ostream &operator<<(
            std::ostream &stream, const NetworkData &networkData);

    void printThroatRadius();
    void printThroatLength();
    void printPoreList();


    std::vector<double> getThroatRadius() const;
    std::vector<double> getThroatLength() const;
    std::vector<int> getPoreList() const;


private:

    std::vector<double> throatRadius;
    std::vector<double> throatLength;

    std::vector<double> connIndIn;
    std::vector<double> connIndOut;

    std::vector<double> poreCoordX;
    std::vector<double> poreCoordY;
    std::vector<double> poreCoordZ;

    std::vector<double> poreRadius;
    std::vector<int> poreList;

public:

};


#endif
