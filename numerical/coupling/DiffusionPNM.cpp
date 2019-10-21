#include <DiffusionPNM.h>

DiffusionPNM::DiffusionPNM(
        const std::vector<double> &propsPNM,
        const std::vector<int> &throat_list,
        const std::vector<double> &throat_radius,
        const std::vector<double> &throat_length,
        const std::vector<double> &conn_ind_in,
        const std::vector<double> &conn_ind_out,
        const std::vector<double> &pore_coord_x,
        const std::vector<double> &pore_coord_y,
        const std::vector<double> &pore_coord_z,
        const std::vector<double> &pore_radius,
        const std::vector<int> &pore_list,
        const std::vector<int> &pore_conns,
        const std::vector<int> &conn_number,
        const std::vector<int> &pore_per_row,
        const std::vector<double> &propsDiffusion) :


        equationPNM(propsPNM, throat_list, throat_radius, throat_length,
                    conn_ind_in, conn_ind_out, pore_coord_x, pore_coord_y,
                    pore_coord_z, pore_radius, pore_list, pore_conns,
                    conn_number, pore_per_row),
        equation(propsDiffusion) {

//    calcRockVolume();
    std::cout << "return 42" << std::endl;
//    std::cout << "rockVolume= " << rockVolume << std::endl;
}

double DiffusionPNM::calcSideLength(std::vector<double> &poreCoord) {

    return 0;
}

//    auto min = std::min_element(std::begin(poreCoord),
//                                std::end(poreCoord));
//    auto max = std::max_element(std::begin(poreCoord),
//                                std::end(poreCoord));
//
//    return *max - *min;



void DiffusionPNM::calcRockVolume() {

    int i = 0;
}

//    auto lengthX = calcSideLength(networkData.poreCoordX);
//    auto lengthY = calcSideLength(networkData.poreCoordY);
//    auto lengthZ = calcSideLength(networkData.poreCoordZ);
//
//    rockVolume = lengthX * lengthY * lengthZ;









