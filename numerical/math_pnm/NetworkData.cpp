#include <NetworkData.h>
#include <vector>

NetworkData::NetworkData(const std::vector<double> &_throat_radius,
                         const std::vector<double> &_throat_length,
                         const std::vector<double> &_conn_ind_in,
                         const std::vector<double> &_conn_ind_out,
                         const std::vector<double> &_pore_coord_x,
                         const std::vector<double> &_pore_coord_y,
                         const std::vector<double> &_pore_coord_z,
                         const std::vector<double> &_pore_radius,
                         const std::vector<int> &_pore_list) :

        throatRadius(_throat_radius),
        throatLength(_throat_length),
        connIndIn(_conn_ind_in),
        connIndOut(_conn_ind_out),
        poreCoordX(_pore_coord_x),
        poreCoordY(_pore_coord_y),
        poreCoordZ(_pore_coord_z),
        poreRadius(_pore_radius),
        poreList(_pore_list){}

std::ostream &operator<<(std::ostream &stream, const NetworkData &networkData) {

    return stream;
}


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


