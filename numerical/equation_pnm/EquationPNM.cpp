#include <EquationPNM.h>

#include <vector>
#include <string>
#include <cmath>

#include <Eigen/Sparse>

#include <PropsPnm.h>
#include <NetworkData.h>


typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Matrix::InnerIterator MatrixIterator;
typedef Eigen::VectorXd Vector;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> BiCGSTAB;
typedef Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>
        LeastSqCG;
typedef Eigen::SparseLU<Eigen::SparseMatrix<double>> SparseLU;

EquationPNM::EquationPNM(PropsPnm &propsPnm,
                         NetworkData &networkData,
                         const std::string &solverMethod) :

        propsPnm(propsPnm),
        networkData(networkData),
        solverMethod(solverMethod),

        dim(networkData.poreN),
        // TODO: think about the expression below
        pIn(propsPnm.pressIn),
        pOut(propsPnm.pressOut),
        matrix(dim, dim),
        freeVector(dim),
        guessVector(dim),
        variable(dim),
        porConnsIsOut(networkData.poreN),
        gammaPnm(networkData.poreN),
        centralCoeff(dim, 0),
        connCoeff(networkData.fracturesN, 0),
        pressure(dim, 0),
        thrFlowRate(networkData.fracturesN, 0),
        porFlowRate(dim, 0) {}

EquationPNM::EquationPNM(const std::vector<double> &propsVector,
                         const std::vector<int> &fracturesList,
                         const std::vector<double> &fracturesHeights,
                         const std::vector<double> &fracturesLengths,
                         const std::vector<double> &fracturesWidths,
                         const std::vector<double> &fracsConnIndIn,
                         const std::vector<double> &fracConnIndOut,
                         const std::vector<double> &poresCoordsX,
                         const std::vector<double> &poresCoordsY,
                         const std::vector<double> &poreCoordZ,
                         const std::vector<double> &poresRadii,
                         const std::vector<int> &poresList,
                         const std::vector<bool> &poresInlet,
                         const std::vector<bool> &poresOutlet,
                         const std::vector<double> &hydraulicCond,
                         const std::string &solverMethod) :

        EquationPNM(*(new PropsPnm(propsVector)),
                    *(new NetworkData(fracturesList, fracturesHeights, fracturesLengths,
                                      fracturesWidths, fracsConnIndIn, fracConnIndOut,
                                      poresCoordsX, poresCoordsY, poreCoordZ,
                                      poresRadii, poresList,
                                      poresInlet,
                                      poresOutlet,
                                      hydraulicCond)),
                    solverMethod) {}


void EquationPNM::calcMatCoeff() {

    for (int i = 0; i < connCoeff.size(); i++)
        connCoeff[i] = 0;
    for (int i = 0; i < centralCoeff.size(); i++)
        centralCoeff[i] = 0;

    for (int i = 0; i < networkData.fracturesN; i++) {

        // TODO: include various options for calculating resistances
        // auto tR = networkData.fracturesHeights[i];
        // auto tH = tR;
        // auto tW = networkData.fracturesWidths[i];
        // auto tL = networkData.fracturesLengths[i];
        // 
        // auto liqVisc = propsPnm.liqVisc;
        // 
        // auto rI = networkData.poresRadii[throatConns[i].first];
        // auto rJ = networkData.poresRadii[throatConns[i].second];

        // Fracture conductance (equal to Yu)
        // connCoeff[i] = tH * tH * tW / 12 / liqVisc / tL;

        // auto cond = (M_PI * tR * tR * tR * tR) / (8 * liqVisc * tL);
        // 3D
        // auto gIJ = (M_PI * tR * tR * tR * tR) / (8 * liqVisc * tL);
        // auto gI = (M_PI * rI * rI * rI) / (8 * liqVisc);
        // auto gJ = (M_PI * rJ * rJ * rJ) / (8 * liqVisc);
        // double dz = 1.E-6;
        // 2D
        // auto gIJ = (dz * 2 * tR * tR * tR) / (3 * liqVisc * tL);
        // auto gI = (dz * 2 * rI * rI) / (3 * liqVisc);
        // auto gJ = (dz * 2 * rJ * rJ) / (3 * liqVisc);

        // 2D reversed y and z (Parallel Plates)
        // auto gIJ = (tR * 2 * dz * dz * dz) / (3. * liqVisc * tL);
        // auto gI = (1. / 2. * rI * 2 * dz * dz * dz) / (3. * liqVisc * rI);
        // auto gJ = (1. / 2. * rJ * 2 * dz * dz * dz) / (3. * liqVisc * rJ);

        // 2D reversed y and z (Parallel Plates)
        // auto gIJ = (tR * dz * dz * dz) / (6. * liqVisc * tL);
        // auto gI = (rI * dz * dz * dz) / (12. * liqVisc * rI);
        // auto gJ = (rJ * dz * dz * dz) / (12. * liqVisc * rJ);

        // 2D reversed y and z (Ellipse Channel)
        // auto gIJ = (tR * 2 * dz * dz * dz) / (3. * liqVisc * tL);
        // auto gI = (1. / 2. * rI * 2 * dz * dz * dz) / (3. * liqVisc * rI);
        // auto gJ = (1. / 2. * rJ * 2 * dz * dz * dz) / (3. * liqVisc * rJ);


        // connCoeff[i] = 1 / (1 / gI + 1 / gIJ + 1 / gJ);

        // Take autmatic conductance from the input file
        connCoeff[i] = networkData.hydraulicCond[i];
    }

    for (int i = 0; i < networkData.por2thrConns.size(); i++)
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++)
            centralCoeff[i] += connCoeff[networkData.por2thrConns[i][j]];
}

void EquationPNM::calculateMatrix(const std::vector<double> &connCoeff,
                                  const std::vector<double> &centralCoeff,
                                  const std::vector<bool> &boundPores,
                                  std::vector<std::vector<int>> &inOutCoeff,
                                  const std::vector<double> &diffCoeff) {

    // Matrix construction
    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);

    for (int i = 0; i < networkData.throatConns.size(); i++) {

        if (!boundPores[networkData.throatConns[i].first]) {
            triplets.emplace_back(networkData.throatConns[i].first,
                                  networkData.throatConns[i].second,
                                  -1. * connCoeff[i]);
        }

        if (!boundPores[networkData.throatConns[i].second]) {
            triplets.emplace_back(networkData.throatConns[i].second,
                                  networkData.throatConns[i].first,
                                  -1. * connCoeff[i]);
        }
    }

    for (int i = 0; i < networkData.poreN; i++) {
        if (!boundPores[i])
            triplets.emplace_back(i, i, centralCoeff[i]);
        else
            triplets.emplace_back(i, i, 1);
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());

    for (int i = 0; i < dim; i++) {
        freeVector[i] = 0;
        guessVector[i] = 0;
        variable[i] = 0;
    }
}

void EquationPNM::calculateFreeVector(const std::string &boundCond,
                                      const double &pIn,
                                      const double &pOut) {

    if (boundCond == "dirichlet") {
        for (int i = 0; i < networkData.poreN; i++)
            if (networkData.poreInlet[i])
                freeVector[i] = pIn;

    } else if (boundCond == "mixed") {
        for (int i = 0; i < networkData.poreN; i++)
            if (networkData.poreInlet[i])
                freeVector[i] = porFlowRate[i];
    }

    for (int i = 0; i < networkData.poreN; i++)
        if (networkData.poreOutlet[i])
            freeVector[i] = pOut;
}

void EquationPNM::calculateGuessPress(const double &pIn,
                                      const double &pOut) {

    auto min = std::min_element(std::begin(networkData.poresCoordsX),
                                std::end(networkData.poresCoordsX));
    auto max = std::max_element(std::begin(networkData.poresCoordsX),
                                std::end(networkData.poresCoordsX));

    for (int i = 0; i < dim; i++)
        pressure[i] = propsPnm.pressOut +
                ((*max - networkData.poresCoordsX[i]) / (*max - *min)) *
                (pIn - pOut);
}

void EquationPNM::calculateGuessVector() {
    for (int i = 0; i < dim; i++)
        guessVector[i] = pressure[i];
}

void EquationPNM::calculatePress(const std::string &solverMethod) {

    if (solverMethod == "biCGSTAB") {

        BiCGSTAB biCGSTAB;
        biCGSTAB.compute(matrix);
        biCGSTAB.setTolerance(propsPnm.itAccuracy);
        variable = biCGSTAB.solveWithGuess(freeVector, guessVector);

    } else if (solverMethod == "sparseLU") {
        SparseLU sparseLU;
        sparseLU.compute(matrix);
        variable = sparseLU.solve(freeVector);

    } else if (solverMethod == "leastSqCG") {
        LeastSqCG leastSqCG;
        leastSqCG.compute(matrix);
        leastSqCG.setTolerance(propsPnm.itAccuracy);
        variable = leastSqCG.solveWithGuess(freeVector, guessVector);
    }

    for (int i = 0; i < dim; i++)
        pressure[i] = variable[i];
}

void EquationPNM::setInitialCondPurePnm() {

    networkData.calcBoundPoresSizes();
    networkData.calcThroatConns();
    networkData.calcPor2ThrConns();
    // networkData.findBoundaryPores(networkData.poresCoordsX);

    getGammaByLabel();

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++) {
            gammaPnm[i].emplace_back(0);
        }
}

void EquationPNM::cfdProcedure(const std::string &boundCond,
                               const std::vector<bool> &boundPores,
                               const double &pIn,
                               const double &pOut) {

    calcMatCoeff();

    // getGammaByPressure();

    std::vector<double> diffCoeff(networkData.fracturesN, 0);

    calculateMatrix(connCoeff,
                    centralCoeff,
                    boundPores,
                    gammaPnm,
                    diffCoeff);

    calculateFreeVector(boundCond, pIn, pOut);

    calculateGuessPress(pIn, pOut);
    calculateGuessVector();

    calculatePress(solverMethod);

    getGammaByPressure();

    calcThrFlowRate();
    calcPorFlowRate();

    calcInletFlow(networkData.poreInlet);
    calcTotFlow(networkData.poreInlet);

//    calculateFreeVector(pIn, pOut);
}

void EquationPNM::cfdProcPurePnmDirichlet() {

    setInitialCondPurePnm();
    std::vector<bool> allBoundaryPores(networkData.poreInlet.size(), 0);

    for (int i = 0; i < allBoundaryPores.size(); i++)
        if (networkData.poreInlet[i] or networkData.poreOutlet[i])
            allBoundaryPores[i] = true;

    cfdProcedure("dirichlet", allBoundaryPores, pIn, pOut);
}

void EquationPNM::calcThrFlowRate() {

    for (int i = 0; i < networkData.fracturesN; i++)
        thrFlowRate[i] =
                connCoeff[i] * abs((pressure[networkData.throatConns[i].first] -
                                    pressure[networkData.throatConns[i].second]));
}

void EquationPNM::getGammaByLabel() {

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++) {
            porConnsIsOut[i].emplace_back(0);
        }

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++) {

            porConnsIsOut[i][j] =
                    networkData.throatConns[networkData.por2thrConns[i][j]].first ==
                    i;
        }
}

void EquationPNM::getGammaByPressure() {

    for (int i = 0; i < porConnsIsOut.size(); i++)
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++) {

            if (porConnsIsOut[i][j]) {
                if (pressure[i] >=
                    pressure[networkData.throatConns[networkData.por2thrConns[i][j]].second])
                    gammaPnm[i][j] = 0;
                else
                    gammaPnm[i][j] = 1;

            } else if (!porConnsIsOut[i][j]) {
                if (pressure[networkData.throatConns[networkData.por2thrConns[i][j]].first] >
                    pressure[i])
                    gammaPnm[i][j] = 1;
                else
                    gammaPnm[i][j] = 0;
            }
        }
}

void EquationPNM::calcPorFlowRate() {

    for (int i = 0; i < porFlowRate.size(); i++)
        porFlowRate[i] = 0;

    for (int i = 0; i < porFlowRate.size(); i++)
        for (int j = 0; j < networkData.por2thrConns[i].size(); j++) {
            if (gammaPnm[i][j] == 0)
                porFlowRate[i] += thrFlowRate[networkData.por2thrConns[i][j]];
            else
                porFlowRate[i] -= thrFlowRate[networkData.por2thrConns[i][j]];
        }
}

void EquationPNM::calcInletFlow(const std::vector<bool> &boundPorIn) {

    inletFlow.clear();
    for (int i = 0; i < boundPorIn.size(); i++)
        if (boundPorIn[i])
            inletFlow.emplace_back(porFlowRate[i]);
}

void EquationPNM::calcTotFlow(const std::vector<bool> &boundPores) {

    totFlowRate = 0;

    for (int i = 0; i < boundPores.size(); i++)
        if (boundPores[i])
            totFlowRate += porFlowRate[i];
}

const std::vector<double> EquationPNM::getPressure() const {
    return pressure;
}

double EquationPNM::getTotFlowRate() const {
    return totFlowRate;
}

//void EquationPNM::calculateMatrix(const std::vector<double> &connCoeff,
//                                  const std::vector<double> &centralCoeff,
//                                  const std::vector<bool> &boundPores,
//                                  std::vector<std::vector<int>> &inOutCoeff,
//                                  const std::vector<double> &diffCoeff) {
//
//    // Matrix construction
//    std::vector<Triplet> triplets;
//    triplets.reserve(3 * dim - 4);
//
//    int pore_iterator = 0;
//
//    for (int i = 0; i < networkData.connNumber.size(); i++) {
//
//        triplets.emplace_back(i, i, centralCoeff[i]);
//
//        for (int j = 0; j < networkData.connNumber[i]; j++) {
//
//            if (!boundPores[i]) {
//
//                triplets.emplace_back(i, networkData.pore2poreConns[pore_iterator],
//                                      -1. *
//                                      connCoeff[networkData.por2thrConns[i][j]]
//                                      + 0. * inOutCoeff[i][j] *
//                                        diffCoeff[networkData.por2thrConns[i][j]]);
//                pore_iterator++;
//
//            } else {
//                triplets.emplace_back(i, i, -1. * centralCoeff[i] + 1);
//                pore_iterator += networkData.connNumber[i];
//                break;
//            }
//        }
//    }
//
//    matrix.setFromTriplets(triplets.begin(), triplets.end());
//    std::cout << matrix << std::endl;
//
//    for (int i = 0; i < dim; i++) {
//        freeVector[i] = 0;
//        guessVector[i] = 0;
//        variable[i] = 0;
//    }
//}