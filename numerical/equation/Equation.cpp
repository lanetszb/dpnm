#include <Equation.h>

Equation::Equation(const std::vector<double> &propsVector) :
        props(propsVector),
        local(propsVector),
        convective(propsVector),
        dim(props.gridBlockN),
        conc(std::vector<std::vector<double>>()),
        time(props.time),
        iCurr(0),
        iPrev(1),
        matrix(dim, dim),
        freeVector(dim),
        guessVector(dim),
        variable(dim) {

    conc.emplace_back(std::vector<double>(dim, props.concIni));
    conc.emplace_back(std::vector<double>(dim, props.concIni));


    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);
    triplets.emplace_back(0, 0, 1);
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

void Equation::calculateMatrix() {

    MatrixIterator(matrix, 0).valueRef() = local.alpha[0];

    for (int i = 1; i < dim - 1; ++i) {
        MatrixIterator it(matrix, i);
        double &betaLeft = convective.beta[Local::left(i)];
        double &betaRight = convective.beta[Local::right(i)];
        double &alpha = local.alpha[i];
        it.valueRef() = -betaLeft;
        ++it;
        it.valueRef() = alpha + betaLeft + betaRight;
        ++it;
        it.valueRef() = -betaRight;
    }

    MatrixIterator it(matrix, dim - 1);
    double &betaLeft = convective.beta[Local::left(dim - 1)];
    double &alpha = local.alpha[dim - 1];
    it.valueRef() = -betaLeft;
    ++it;
    it.valueRef() = alpha + betaLeft;
}

void Equation::calculateFreeVector(const double &conc_in) {
    freeVector[0] = local.alpha[0] * conc_in;
    for (int i = 1; i < dim; i++)
        freeVector[i] = local.alpha[i] * conc[iPrev][i];
}

void Equation::calculateGuessVector() {
    for (int i = 0; i < dim; i++)
        guessVector[i] = conc[iPrev][i];
}

void Equation::calculateConc() {

    BiCGSTAB biCGSTAB;

    biCGSTAB.compute(matrix);

    variable = biCGSTAB.solveWithGuess(freeVector, guessVector);

    for (int i = 0; i < dim; i++)
        conc[iCurr][i] = variable[i];

}

void Equation::cfdProcedure(const double &concIn) {

    for (double t = props.timeStep; t <= props.time; t += props.timeStep) {

        std::swap(iCurr, iPrev);
        local.calculateAlpha(props.timeStep,
                             props.radius,
                             props.effRadius);
        convective.calculateBeta(props.radius,
                                 props.effRadius,
                                 props.length,
                                 props.diffusivity,
                                 props.gridBlockN);
        calculateGuessVector();
        calculateMatrix();
        calculateFreeVector(concIn);
        calculateConc();
    }

}


const std::vector<double> Equation::getConc() const {
    return conc[iCurr];
}

