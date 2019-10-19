#include <EquationPNM.h>
#include <math.h>

EquationPNM::EquationPNM(const std::vector<double> &_propsVector,
                         const std::vector<int> &_throat_list,
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
        networkData(_throat_list, _throat_radius, _throat_length,
                    _conn_ind_in, _conn_ind_out, _pore_coord_x, _pore_coord_y,
                    _pore_coord_z, _pore_radius,
                    _pore_list, _pore_conns, _conn_number, _pore_per_row),
        dim(networkData.poreN),
        iCurr(0),
        iPrev(1),
        pIn(300000),
        pOut(200000),
        matrix(dim, dim),
        freeVector(dim),
        guessVector(dim),
        variable(dim),
        throatConns(networkData.throatN,
                    std::make_pair(-1, -1)),
        porConns(networkData.poreN),
        centralCoeff(dim, 0) {

    press.emplace_back(std::vector<double>(dim, 0));
    press.emplace_back(std::vector<double>(dim, 0));

    cfdProcedure(pIn, pOut);


    for (int i = 0; i < dim; i++)
        std::cout << press[iCurr][i] << std::endl;
}


void EquationPNM::calcThroatConns() {

    for (int i = 0; i < networkData.throatN; i++) {
        throatConns[i].first = networkData.connIndIn[i];
        throatConns[i].second = networkData.connIndOut[i];
    }
}

void EquationPNM::calcPorConns() {

    // Ugly construction, has to be rewritten later, but works well for now

    int pore_iterator = 0;
    for (int i = 0; i < networkData.poreN; i++)
        for (int j = 0; j < networkData.connNumber[i]; j++) {
            porConns[i].emplace_back(networkData.poreConns[pore_iterator]);
            pore_iterator++;
        }

    for (int i = 0; i < porConns.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++)
            for (int k = 0; k < throatConns.size(); k++) {
                if ((i == throatConns[k].first and
                     porConns[i][j] == throatConns[k].second) or
                    (porConns[i][j] == throatConns[k].first and
                     i == throatConns[k].second)) {
                    porConns[i][j] = k;
                    break;
                }
            }
}

void EquationPNM::calcMatCoeff() {

    for (int i = 0; i < networkData.throatN; i++) {

        auto tR = networkData.throatRadius[i];
        auto tL = networkData.throatLength[i];
        auto liqVisc = propsPnm.liqVisc;

        auto cond = (M_PI * tR * tR * tR * tR) / (8 * liqVisc * tL);

        connCoeff.emplace_back(cond);
    }

    for (int i = 0; i < porConns.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++)
            centralCoeff[i] += connCoeff[porConns[i][j]];
}

void EquationPNM::calculateMatrix() {

    // Matrix construction

    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);

    int pore_iterator = 0;
    int bound_it = 0;

    for (int i = 0; i < networkData.connNumber.size(); i++) {

        triplets.emplace_back(i, i, centralCoeff[i]);
        for (int j = 0; j < networkData.connNumber[i]; j++) {

            if (i != networkData.boundaryPores[bound_it]) {
                triplets.emplace_back(i,
                                      networkData.poreConns[pore_iterator],
                                      -1 *
                                      connCoeff[porConns[i][j]]);
                pore_iterator++;
            } else {
                triplets.emplace_back(i, i, 1);
                pore_iterator += networkData.connNumber[i];
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
}

void EquationPNM::calculateFreeVector(const double &pIn,
                                      const double &pOut) {

    for (int i = 0; i < networkData.boundaryPoresIn.size(); i++)
        freeVector[i] = pIn;

    for (int i = networkData.poreN - 1;
         i > networkData.poreN - 1 - networkData.boundaryPoresOut.size(); i--)
        freeVector[i] = pOut;
}

void EquationPNM::calculateGuessPress(const double &pIn,
                                      const double &pOut) {
    for (int i = 0; i < dim; i++)
        press[iCurr][i] = pIn * (dim - 1 - i) / (dim - 1) +
                          pOut * i / (dim - 1);
}

void EquationPNM::calculateGuessVector() {
    for (int i = 0; i < dim; i++)
        guessVector[i] = press[iPrev][i];
}

void EquationPNM::calculatePress() {

//    BiCGSTAB biCGSTAB;
//    biCGSTAB.compute(matrix);
//    variable = biCGSTAB.solve(freeVector, guessVector);

    SparseLU sparseLU;
    sparseLU.compute(matrix);
    variable = sparseLU.solve(freeVector);

    for (int i = 0; i < dim; i++)
        press[iCurr][i] = variable[i];
}

void EquationPNM::cfdProcedure(const double &pIn,
                               const double &pOut) {

    calcThroatConns();
    calcPorConns();
    calcMatCoeff();
    networkData.findBoundaryPores(networkData.poreCoordX);

    calculateMatrix();
    calculateFreeVector(pIn, pOut);
    calculateGuessPress(pIn, pOut);
    calculateGuessVector();

    calculatePress();
}

double EquationPNM::calculatePressRelDiff() {
    double relDiff = 0;
    for (int i = 0; i < dim; i++)
        relDiff += fabs(press[iCurr][i] - press[iPrev][i]) / press[iCurr][i] /
                   dim;
    return relDiff;
}