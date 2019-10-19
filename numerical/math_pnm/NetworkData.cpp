#include <NetworkData.h>

NetworkData::NetworkData(const std::vector<int> &_throat_list,
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
                         const std::vector<int> &_pore_per_row) :

        throatList(_throat_list),
        throatRadius(_throat_radius),
        throatLength(_throat_length),
        connIndIn(_conn_ind_in),
        connIndOut(_conn_ind_out),
        throatN(throatList.size()),

        poreCoordX(_pore_coord_x),
        poreCoordY(_pore_coord_y),
        poreCoordZ(_pore_coord_z),

        poreRadius(_pore_radius),
        poreList(_pore_list),
        poreN(poreList.size()),

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
        if (poreCoord[i] == *min)
            boundaryPoresIn.emplace_back(i);
        else if (poreCoord[i] == *max)
            boundaryPoresOut.emplace_back(i);
    }

    for (int i = 0; i < boundaryPoresIn.size(); i++)
        boundaryPores.emplace_back(boundaryPoresIn[i]);

    boundaryPores.insert(boundaryPores.end(),
                         boundaryPoresOut.begin(), boundaryPoresOut.end());
}


