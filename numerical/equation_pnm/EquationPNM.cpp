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
                         const std::vector<int> &pore_per_row,
                         const std::vector<bool> &pore_left_x,
                         const std::vector<bool> &pore_right_x,
                         const std::vector<double> &hydraulic_cond) :

        propsPNM(propsVector),
        networkData(throat_list, throat_radius, throat_length,
                    conn_ind_in, conn_ind_out, pore_coord_x, pore_coord_y,
                    pore_coord_z, pore_radius, pore_list, pore_conns,
                    conn_number, pore_per_row, pore_left_x, pore_right_x,
                    hydraulic_cond),
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
        porConnsIsOutByPressure(networkData.poreN),
        centralCoeff(dim, 0),
        connCoeff(networkData.throatN, 0),
        pressure(dim, 0),
        thrFlowRate(networkData.throatN, 0),
        porFlowRate(dim, 0) {


//    setInitialCond();
//
//
//    cfdProcedure(1, networkData.boundaryPores, pIn, pOut);
//
//    std::cout << "completed" << std::endl;

//    calcTotFlow(networkData.boundaryPores);

//    for (int i = 0; i < networkData.poreN; i++)
//        std::cout << pressure[i] << std::endl;
//
//    std::cout << std::endl;


//    auto min = std::min_element(std::begin(networkData.poreCoordX),
//                                std::end(networkData.poreCoordX));
//
//    auto max = std::max_element(std::begin(networkData.poreCoordX),
//                                std::end(networkData.poreCoordX));
//
//    for (int i = 0; i < connCoeff.size(); i++)
//        std::cout << connCoeff[i] << std::endl;
//
//    std::cout << std::endl;

//    auto length = (*max - *min);
//    auto area = length * length;

//    double length = 1.E-6 * 1000;
//    double area = 1.E-6 * 1000 * 1.E-6;
//
//    std::cout << "Q= " << totFlowRate << std::endl;
//    std::cout << "length " << length << std::endl;
//    std::cout << "area " << area << std::endl;
//    std::cout << "dP " << propsPNM.pressIn - propsPNM.pressOut << std::endl;
//
//    auto perm = totFlowRate * propsPNM.liqVisc * length / area /
//                (propsPNM.pressIn - propsPNM.pressOut);
//
//    std::cout << "perm =" << perm << std::endl;
    //
//    std::cout << std::endl;
//
//    std::cout << totFlowRate << std::endl;
//
//    std::cout << matrix << std::endl;
//
//    std::cout << std::endl;
//
//    cfdProcedure(0, networkData.boundaryPoresOut, pIn, pOut);
//
//    for (int i = 0; i < networkData.poreN; i++)
//        std::cout << pressure[i] << std::endl;

//    std::cout << matrix << std::endl;


//    for (int i = 0; i < inletFlow.size(); i++)
//        std::cout << inletFlow[i] << std::endl;

//    std::cout << matrix << std::endl;
//

//
//    std::cout << "porConns" << std::endl;
//    for (int i = 0; i < porConns.size(); i++) {
//        for (int j = 0; j < porConns[i].size(); j++) {
//            std::cout << porConns[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    std::cout << "throatConns" << std::endl;
//    for (int i = 0; i < networkData.throatN; i++) {
//        std::cout << throatConns[i].first << ' ' << throatConns[i].second
//                  << std::endl;
//    }

//    std::cout << "porConnsIsOut" << std::endl;
//    for (int i = 0; i < porConnsIsOutByPressure.size(); i++) {
//        for (int j = 0; j < networkData.connNumber[i]; j++) {
//            std::cout << porConnsIsOut[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << "porConnsIsOutByPressure" << std::endl;
//    for (int i = 0; i < porConnsIsOutByPressure.size(); i++) {
//        for (int j = 0; j < networkData.connNumber[i]; j++) {
//            std::cout << porConnsIsOutByPressure[i][j] << ' ';
//        }
//        std::cout << std::endl;
//    }

//    std::cout << "throatFlow" << std::endl;
//    for (int i = 0; i < networkData.throatN; i++) {
//        std::cout << i << ' ' << thrFlowRate[i] << std::endl;
//    }
//    std::cout << std::endl;


//    std::cout << "FreeVector" << std::endl;
//    std::cout << freeVector << std::endl;
//    std::cout << std::endl;
//
//    std::cout << "matrix" << std::endl;
//    std::cout << matrix << std::endl;
//    std::cout << std::endl;
//
//    std::cout << std::endl;
//

//
//
//    for (int i = 0; i < networkData.throatN; i++)
//        std::cout << thrFlowRate[i] << std::endl;

//    std::cout << std::endl;
//
//    for (int i = 0; i < dim; i++)
//        std::cout << i << ' ' << porFlowRate[i] << std::endl;
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

    for (int i = 0; i < connCoeff.size(); i++)
        connCoeff[i] = 0;
    for (int i = 0; i < centralCoeff.size(); i++)
        centralCoeff[i] = 0;

    for (int i = 0; i < networkData.throatN; i++) {

        auto tR = networkData.throatRadius[i];
        auto tL = networkData.throatLength[i];
        auto liqVisc = propsPNM.liqVisc;

        auto rI = networkData.poreRadius[throatConns[i].first];
        auto rJ = networkData.poreRadius[throatConns[i].second];

        // auto cond = (M_PI * tR * tR * tR * tR) / (8 * liqVisc * tL);
        // 3D
        // auto gIJ = (M_PI * tR * tR * tR * tR) / (8 * liqVisc * tL);
        // auto gI = (M_PI * rI * rI * rI) / (8 * liqVisc);
        // auto gJ = (M_PI * rJ * rJ * rJ) / (8 * liqVisc);
        double dz = 1.E-6;
        // 2D
        // auto gIJ = (dz * 2 * tR * tR * tR) / (3 * liqVisc * tL);
        // auto gI = (dz * 2 * rI * rI) / (3 * liqVisc);
        // auto gJ = (dz * 2 * rJ * rJ) / (3 * liqVisc);

        // 2D reversed y and z (Parallel Plates)
         auto gIJ = (tR * 2 * dz * dz * dz) / (3. * liqVisc * tL);
         auto gI = (1. / 2. * rI * 2 * dz * dz * dz) / (3. * liqVisc * rI);
         auto gJ = (1. / 2. * rJ * 2 * dz * dz * dz) / (3. * liqVisc * rJ);

        // 2D reversed y and z (Parallel Plates)
        // auto gIJ = (tR * dz * dz * dz) / (6. * liqVisc * tL);
        // auto gI = (rI * dz * dz * dz) / (12. * liqVisc * rI);
        // auto gJ = (rJ * dz * dz * dz) / (12. * liqVisc * rJ);

        // 2D reversed y and z (Ellipse Channel)
        // auto gIJ = (tR * 2 * dz * dz * dz) / (3. * liqVisc * tL);
        // auto gI = (1. / 2. * rI * 2 * dz * dz * dz) / (3. * liqVisc * rI);
        // auto gJ = (1. / 2. * rJ * 2 * dz * dz * dz) / (3. * liqVisc * rJ);


        connCoeff[i] = 1 / (1 / gI + 1 / gIJ + 1 / gJ);
        // connCoeff[i] = networkData.hydraulicCond[i];
    }

    for (int i = 0; i < porConns.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++)
            centralCoeff[i] += connCoeff[porConns[i][j]];
}

void EquationPNM::calculateMatrix(const int &boundCond,
                                  const std::vector<double> &connCoeff,
                                  const std::vector<double> &centralCoeff,
                                  const std::vector<int> &boundPores,
                                  std::vector<std::vector<int>> &inOutCoeff,
                                  const std::vector<double> &diffCoeff) {

    // Matrix construction

    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);

    int pore_iterator = 0;
    int bound_it = 0;

    for (int i = 0; i < networkData.connNumber.size(); i++) {

        triplets.emplace_back(i, i, centralCoeff[i]);
        for (int j = 0; j < networkData.connNumber[i]; j++) {

            if ((i != boundPores[bound_it] and boundCond == 0) or
                (i != boundPores[bound_it] and boundCond == 1)) {

                triplets.emplace_back(i, networkData.poreConns[pore_iterator],
                                      -1 * connCoeff[porConns[i][j]]
                                      - 1 * inOutCoeff[i][j] *
                                        diffCoeff[porConns[i][j]]);
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

void EquationPNM::calculateFreeVector(const int &boundCond,
                                      const double &pIn,
                                      const double &pOut) {

    if (boundCond == 1) {
        for (int i = 0; i < networkData.poreN; i++)
            if (networkData.poreLeftX[i])
                freeVector[i] = pIn;
//        for (int i = 0; i < networkData.boundaryPoresIn.size(); i++)
//            freeVector[i] = pIn;

    } else {
        for (int i = 0; i < networkData.poreN; i++)
            if (networkData.poreLeftX[i])
                freeVector[i] = inletFlow[i];

//        for (int i = 0; i < networkData.boundaryPoresIn.size(); i++)
//            freeVector[i] = inletFlow[i];
    }
    for (int i = 0; i < networkData.poreN; i++)
        if (networkData.poreRightX[i])
            freeVector[i] = pOut;


//    int boundPoreSize = networkData.poreN - 1 -
//                        networkData.boundaryPoresOut.size();
//
//    for (int i = networkData.poreN - 1; i > boundPoreSize; i--)
//        freeVector[i] = pOut;
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

void EquationPNM::calculatePress(const int &solverMethod) {

    if (solverMethod == 0) {

        BiCGSTAB biCGSTAB;
        biCGSTAB.compute(matrix);
        biCGSTAB.setTolerance(propsPNM.itAccuracy);
        variable = biCGSTAB.solveWithGuess(freeVector, guessVector);
    }

    if (solverMethod == 1) {
        SparseLU sparseLU;
        sparseLU.compute(matrix);
        variable = sparseLU.solve(freeVector);
    }

//    LeastSqCG leastSqCG;
//    leastSqCG.compute(matrix);
//    leastSqCG.setTolerance(propsPNM.itAccuracy);
//    variable = leastSqCG.solveWithGuess(freeVector, guessVector);

    for (int i = 0; i < dim; i++)
        pressure[i] = variable[i];
}

void EquationPNM::setInitialCond() {

    calcThroatConns();
    calcPorConns();
    networkData.findBoundaryPores(networkData.poreCoordX);

    getPorConnsIsOut();

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < networkData.connNumber[i]; j++) {
            porConnsIsOutByPressure[i].emplace_back(0);
        }
}

void EquationPNM::cfdProcedure(const int &boundCond,
                               const std::vector<int> &boundPores,
                               const double &pIn,
                               const double &pOut) {


    calcMatCoeff();

//    getPorConnsIsOutByPressure();

    std::vector<double> diffCoeff(networkData.throatN, 0);

    std::vector<int> boundPoresInput;

    for (int i = 0; i < networkData.poreN; i++)
        if (networkData.poreLeftX[i] or networkData.poreRightX[i])
            boundPoresInput.emplace_back(i);

    calculateMatrix(boundCond,
                    connCoeff,
                    centralCoeff,
                    boundPoresInput,
                    porConnsIsOutByPressure,
                    diffCoeff);

    calculateFreeVector(boundCond, pIn, pOut);

    calculateGuessPress(pIn, pOut);
    calculateGuessVector();

    calculatePress(1);

    getPorConnsIsOutByPressure();

    calcThrFlowRate();
    calcPorFlowRate();


    calcInletFlow(networkData.poreRightX.size());

    std::vector<int> boundPoresIn;

    for (int i = 0; i < networkData.poreN; i++)
        if (networkData.poreLeftX[i])
            boundPoresIn.emplace_back(i);

    calcTotFlow(boundPoresIn);

//    calculateFreeVector(pIn, pOut);
}

void EquationPNM::calcThrFlowRate() {

    for (int i = 0; i < networkData.throatN; i++)
        thrFlowRate[i] = connCoeff[i] * abs((pressure[throatConns[i].first] -
                                             pressure[throatConns[i].second]));
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

void EquationPNM::getPorConnsIsOutByPressure() {

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++) {

            if (porConnsIsOut[i][j]) {
                if (pressure[i] >= pressure[throatConns[porConns[i][j]].second])
                    porConnsIsOutByPressure[i][j] = 0;
                else
                    porConnsIsOutByPressure[i][j] = 1;

            } else if (!porConnsIsOut[i][j]) {
                if (pressure[throatConns[porConns[i][j]].first] > pressure[i])
                    porConnsIsOutByPressure[i][j] = 1;
                else
                    porConnsIsOutByPressure[i][j] = 0;
            }
        }
}

void EquationPNM::calcPorFlowRate() {

    for (int i = 0; i < porFlowRate.size(); i++)
        porFlowRate[i] = 0;

    for (int i = 0; i < porFlowRate.size(); i++)
        for (int j = 0; j < porConns[i].size(); j++) {
            if (porConnsIsOutByPressure[i][j] == 0)
                porFlowRate[i] += thrFlowRate[porConns[i][j]];
            else
                porFlowRate[i] -= thrFlowRate[porConns[i][j]];
        }
}

void EquationPNM::calcInletFlow(const int &boundPorSize) {

    for (int i = 0; i < boundPorSize; i++)
        inletFlow.emplace_back(porFlowRate[i]);
}

void EquationPNM::calcTotFlow(const std::vector<int> &boundPores) {

    totFlowRate = 0;

    for (int i = 0; i < boundPores.size(); i++)
        for (int j = 0; j < networkData.poreList.size(); j++)
            if (boundPores[i] == networkData.poreList[j]) {
                totFlowRate += porFlowRate[j];
            }
}