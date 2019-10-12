#include <EquationPNM.h>

EquationPNM::EquationPNM(const std::vector<double> &_propsVector,
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
        propsPnm(_propsVector),
        networkData(_throat_radius, _throat_length, _conn_ind_in, _conn_ind_out,
                    _pore_coord_x, _pore_coord_y, _pore_coord_z, _pore_radius,
                    _pore_list, _pore_conns, _conn_number, _pore_per_row),
        dim(networkData.connNumber.size()),
        matrix(dim, dim),
        freeVector(dim),
        guessVector(dim),
        variable(dim) {

    networkData.findBoundaryPores(networkData.poreCoordX);

    for (int i = 0; i < networkData.boundaryPores.size(); i++)
        std::cout << networkData.boundaryPores[i] << ' ';
    std::cout << std::endl;

    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);


    int pore_it = 0;
    int bound_it = 0;

    for (int i = 0; i < networkData.connNumber.size(); i++) {
        for (int j = 0; j < (networkData.connNumber[i] + 1); j++) {
            if (i != networkData.boundaryPores[bound_it]) {
                triplets.emplace_back(i, networkData.porPerRow[pore_it], 7);
                pore_it++;

            } else {
                triplets.emplace_back(i, i, 1);
                pore_it += (networkData.connNumber[i] + 1);
                bound_it++;
                break;
            }
        }
    }
    matrix.setFromTriplets(triplets.begin(), triplets.end());

    for (int i = 0; i < dim; i++) {
        freeVector[i] = 0;
        guessVector[i] = 0;
        variable[i] = 0;
    }
    std::cout << std::endl;
    std::cout << matrix << std::endl;
}


