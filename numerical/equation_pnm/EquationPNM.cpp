#include <EquationPNM.h>
#include <math.h>

EquationPNM::EquationPNM(const std::vector<double> &propsVector,
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
                         const std::vector<int> &pore_per_row) :

        propsPNM(propsVector),
        networkData(throat_list, throat_radius, throat_length,
                    conn_ind_in, conn_ind_out, pore_coord_x, pore_coord_y,
                    pore_coord_z, pore_radius,
                    pore_list, pore_conns, conn_number, pore_per_row),
        dim(networkData.poreN),
        pIn(propsPNM.pressIn),
        pOut(propsPNM.pressOut),
        matrix(dim, dim),
        freeVector(dim),
        guessVector(dim),
        variable(dim),
        throatConns(networkData.throatN,
                    std::make_pair(-1, -1)),
        porConns(networkData.poreN),
        porConnsIsOut(networkData.poreN),
        centralCoeff(dim, 0),
        pressure(dim, 0),
        thrFlowRate(networkData.throatN, 0),
        porFlowRate(dim, 0) {

    cfdProcedure(pIn, pOut);
    calcThrFlowRate();
    getPorConnsIsOut();
    calcPorFlowRate();


    // for (int i = 0; i < dim; i++)
    //     std::cout << pressure[i] << std::endl;
    //
    // for (int i = 0; i < dim; i++)
    //     std::cout << i << ' ' << porFlowRate[i] << std::endl;
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
        auto liqVisc = propsPNM.liqVisc;

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

    auto min = std::min_element(std::begin(networkData.poreCoordX),
                                std::end(networkData.poreCoordX));
    auto max = std::max_element(std::begin(networkData.poreCoordX),
                                std::end(networkData.poreCoordX));

    for (int i = 0; i < dim; i++)
        pressure[i] = propsPNM.pressOut +
                      ((*max - networkData.poreCoordX[i]) / (*max - *min)) *
                      (pIn - pOut);
}

void EquationPNM::calculateGuessVector() {
    for (int i = 0; i < dim; i++)
        guessVector[i] = pressure[i];
}

void EquationPNM::calculatePress() {

    BiCGSTAB biCGSTAB;
    biCGSTAB.compute(matrix);
    biCGSTAB.setTolerance(propsPNM.itAccuracy);
    variable = biCGSTAB.solveWithGuess(freeVector, guessVector);

    // SparseLU sparseLU;
    // sparseLU.compute(matrix);
    // variable = sparseLU.solve(freeVector);

//    LeastSqCG leastSqCG;
//    leastSqCG.compute(matrix);
//    leastSqCG.setTolerance(propsPNM.itAccuracy);
//    variable = leastSqCG.solveWithGuess(freeVector, guessVector);

    for (int i = 0; i < dim; i++)
        pressure[i] = variable[i];
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

void EquationPNM::calcThrFlowRate() {

    for (int i = 0; i < networkData.throatN; i++)
        thrFlowRate[i] = connCoeff[i] * (pressure[throatConns[i].second] -
                                         pressure[throatConns[i].first]);

    // std::cout << "Debit" << std::endl;
    //
    // for (int i = 0; i < networkData.throatN; i++)
    //     std::cout << pressure[throatConns[i].second] << ' '
    //               << pressure[throatConns[i].first] << ' '
    //               << (pressure[throatConns[i].second] -
    //                   pressure[throatConns[i].first]) << ' ' << thrFlowRate[i]
    //               << std::endl;
    //
    // std::cout << std::endl;
}

void EquationPNM::getPorConnsIsOut() {

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++) {
            porConnsIsOut[i].emplace_back(0);
        }

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++) {

            if (throatConns[porConns[i][j]].first == i)
                porConnsIsOut[i][j] = true;
            else
                porConnsIsOut[i][j] = false;
        }
}

void EquationPNM::calcPorFlowRate() {

    for (int i = 0; i < porFlowRate.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++) {
            if (porConnsIsOut[i][j])
                porFlowRate[i] += thrFlowRate[porConns[i][j]];
            else
                porFlowRate[i] -= thrFlowRate[porConns[i][j]];
        }
}