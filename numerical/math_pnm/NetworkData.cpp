#include <NetworkData.h>

NetworkData::NetworkData(const std::vector<double> &_throat_radius,
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
                         const std::vector<int> &_pore_per_row) :

        throatRadius(_throat_radius),
        throatLength(_throat_length),
        connIndIn(_conn_ind_in),
        connIndOut(_conn_ind_out),
        poreCoordX(_pore_coord_x),
        poreCoordY(_pore_coord_y),
        poreCoordZ(_pore_coord_z),
        poreRadius(_pore_radius),
        poreList(_pore_list),
        poreConns(_pore_conns),
        connNumber(_conn_number),
        porPerRow(_pore_per_row) {}

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

void NetworkData::findBoundaryPores(std::vector<double> &poreCoord) {

    auto min = std::min_element(std::begin(poreCoord), std::end(poreCoord));
    auto max = std::max_element(std::begin(poreCoord), std::end(poreCoord));

    for (int i = 0; i < poreCoord.size(); i++) {
        if (poreCoord[i] == *min or poreCoord[i] == *max)
        boundaryPores.emplace_back(i);
    }

}


