#include <EquationDiffusion.h>

EquationDiffusion::EquationDiffusion(const std::vector<double> &propsVector) :
        propsDiffusion(propsVector),
        localDiffusion(propsVector),
        convectiveDiffusion(propsVector),
        dim(propsDiffusion.gridBlockN),
        conc(std::vector<std::vector<double>>()),
        time(propsDiffusion.time),
        iCurr(0), iPrev(1),
        matrix(dim, dim),
        freeVector(dim),
        guessVector(dim),
        variable(dim) {

    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);
    triplets.emplace_back(0, 0);
    for (int i = 1; i < dim - 1; i++) {
        triplets.emplace_back(i, i - 1);
        triplets.emplace_back(i, i);
        triplets.emplace_back(i, i + 1);
    }
    triplets.emplace_back(dim - 1, dim - 2);
    triplets.emplace_back(dim - 1, dim - 1);
    matrix.setFromTriplets(triplets.begin(), triplets.end());

    for (int i = 0; i < dim; i++) {
        freeVector[i] = 0;
        guessVector[i] = 0;
        variable[i] = 0;
    }

}

void EquationDiffusion::calculateMatrix() {

    MatrixIterator(matrix, 0).valueRef() = localDiffusion.alpha[0];

    for (int i = 1; i < dim - 1; ++i) {
        MatrixIterator it(matrix, i);
        double &betaLeft = convectiveDiffusion.beta[LocalDiffusion::left(i)];
        double &betaRight = convectiveDiffusion.beta[LocalDiffusion::right(i)];
        double &alpha = localDiffusion.alpha[i];
        it.valueRef() = -1 * betaLeft;
        ++it;
        it.valueRef() = alpha + betaLeft + betaRight;
        ++it;
        it.valueRef() = -1 * betaRight;
    }

    MatrixIterator it(matrix, dim - 1);
    double &betaLeft = convectiveDiffusion.beta[LocalDiffusion::left(dim - 1)];
    double &alpha = localDiffusion.alpha[dim - 1];
    it.valueRef() = -1 * betaLeft;
    ++it;
    it.valueRef() = alpha + betaLeft;
}

void EquationDiffusion::calcConcIni(const double &concIni) {

    conc.emplace_back(std::vector<double>(dim, concIni));
    conc.emplace_back(std::vector<double>(dim, concIni));
}


void EquationDiffusion::forceDirichletBound(const double &conc_in) {

    MatrixIterator it(matrix, dim - 1);
    it.valueRef() = 0;
    ++it;
    it.valueRef() = localDiffusion.alpha[dim - 1];

    freeVector[dim - 1] = localDiffusion.alpha[dim - 1] * 5 * conc_in;

}

void EquationDiffusion::calculateFreeVector(const double &conc_in) {
    freeVector[0] = localDiffusion.alpha[0] * conc_in;
    for (int i = 1; i < dim; i++)
        freeVector[i] = localDiffusion.alpha[i] * conc[iPrev][i];
}

void EquationDiffusion::calculateGuessVector() {
    for (int i = 0; i < dim; i++)
        guessVector[i] = conc[iPrev][i];
}

void EquationDiffusion::calculateConc() {

    BiCGSTAB biCGSTAB;

    biCGSTAB.compute(matrix);

    variable = biCGSTAB.solveWithGuess(freeVector, guessVector);

    for (int i = 0; i < dim; i++)
        conc[iCurr][i] = variable[i];

}

void EquationDiffusion::calcFlowRate() {
    flowRate = -convectiveDiffusion.beta[1] * (conc[iCurr][0] - conc[iCurr][1]);
}

void EquationDiffusion::calcTimeVector() {

    auto time = propsDiffusion.time;
    auto configTimeStep = propsDiffusion.timeStep;
    double division = time / configTimeStep;
    double fullStepsN;
    auto lastStep = std::modf(division, &fullStepsN);
    std::cout << lastStep << std::endl;
    auto timeSteps = std::vector<double>(fullStepsN, configTimeStep);
    if (lastStep > 0)
        timeSteps.push_back(lastStep * configTimeStep);
    timeStepsVec = timeSteps;
    for (auto &&timeStep : timeSteps)
        std::cout << timeStep << std::endl;
}

void EquationDiffusion::cfdProcedureOneStep(const std::string &boundCond,
                                            const double &concThrWall,
                                            const double &radius,
                                            const double &effRadius,
                                            const double &thrLength,
                                            const std::vector<double> &volumes,
                                            const std::vector<double> &surfaces,
                                            const double &dt) {

    std::swap(iCurr, iPrev);
    localDiffusion.calculateAlpha(propsDiffusion.timeStep,
                                  volumes);
    convectiveDiffusion.calculateBeta(radius,
                                      effRadius,
                                      thrLength,
                                      propsDiffusion.diffusivity,
                                      propsDiffusion.gridBlockN,
                                      surfaces);

    calculateGuessVector();
    calculateMatrix();
    calculateFreeVector(concThrWall);
    if (boundCond == "Dirichlet")
        forceDirichletBound(concThrWall);
    calculateConc();
    calcFlowRate();
}

void EquationDiffusion::cfdProcedure(const std::string &boundCond,
                                     const double &concThrWall,
                                     const double &radius,
                                     const double &effRadius,
                                     const double &thrLength,
                                     const std::vector<double> &volumes,
                                     const std::vector<double> &surfaces) {

    calcTimeVector();
    calcConcIni(propsDiffusion.concIni);

    for (int i = 0; i <= timeStepsVec.size(); i++) {

        cfdProcedureOneStep(boundCond, concThrWall, radius, effRadius,
                            thrLength, volumes, surfaces, timeStepsVec[i]);
    }
}

const std::vector<double> EquationDiffusion::getConc() const {
    return conc[iCurr];
}

const double EquationDiffusion::getFlowRate() const {
    return flowRate;
}


