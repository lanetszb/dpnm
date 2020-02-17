#include <EquationDiffusion.h>

EquationDiffusion::EquationDiffusion(const std::vector<double> &propsVector,
                                     const std::vector<double> &langmuirCoeff) :
        props(propsVector, langmuirCoeff),
        local(propsVector, langmuirCoeff),
        convective(propsVector, langmuirCoeff),
        dim(props.gridBlockN),
        conc(std::vector<std::vector<double>>()),
        time(props.time),
        iCurr(0),
        iPrev(1),
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

    MatrixIterator(matrix, 0).valueRef() = local.alpha[0];

    for (int i = 1; i < dim - 1; ++i) {
        MatrixIterator it(matrix, i);
        double &betaLeft = convective.beta[LocalDiffusion::left(i)];
        double &betaRight = convective.beta[LocalDiffusion::right(i)];
        double &alpha = local.alpha[i];
        it.valueRef() = -1 * betaLeft;
        ++it;
        it.valueRef() = alpha + betaLeft + betaRight;
        ++it;
        it.valueRef() = -1 * betaRight;
    }

    MatrixIterator it(matrix, dim - 1);
    double &betaLeft = convective.beta[LocalDiffusion::left(dim - 1)];
    double &alpha = local.alpha[dim - 1];
    it.valueRef() = -1 * betaLeft;
    ++it;
    it.valueRef() = alpha + betaLeft;
}

void EquationDiffusion::calcConcIni(const double &concIni) {

    conc.emplace_back(std::vector<double>(dim, concIni));
    conc.emplace_back(std::vector<double>(dim, concIni));
}

void EquationDiffusion::calculateFreeVector(const double &conc_in) {
    freeVector[0] = local.alpha[0] * conc_in;
    for (int i = 1; i < dim; i++)
        freeVector[i] = local.alpha[i] * conc[iPrev][i];
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

// Double check how the flowrate is obtained

void EquationDiffusion::calcFlowRate() {
    flowRate = -1 * (conc[iCurr][0] - conc[iCurr][1]) * convective.beta[1];
}

void EquationDiffusion::cfdProcedureOneStep(const double &concThrWall,
                                            const double &radius,
                                            const double &effRadius,
                                            const double &thrLength) {

    std::swap(iCurr, iPrev);

    local.calcVolCylindr(radius,
                         effRadius,
                         props.gridBlockN,
                         thrLength);

    local.calculateAlpha(props.timeStep,
                         local.volCylindr);

    convective.calcOmegaCylindr(thrLength);

    convective.calculateBeta(radius,
                             effRadius,
                             thrLength,
                             props.diffusivity,
                             props.gridBlockN,
                             convective.omegaCylindric);

    calculateGuessVector();
    calculateMatrix();
    calculateFreeVector(concThrWall);
    calculateConc();
    calcFlowRate();
}

void EquationDiffusion::cfdProcedure(const double &concThrWall,
                                     const double &radius,
                                     const double &effRadius,
                                     const double &thrLength) {

    calcConcIni(props.concIni);

    for (double t = props.timeStep; t <= props.time; t += props.timeStep) {
        cfdProcedureOneStep(concThrWall, radius, effRadius, thrLength);
    }
}

const std::vector<double> EquationDiffusion::getConc() const {
    return conc[iCurr];
}

const double EquationDiffusion::getFlowRate() const {
    return flowRate;
}


