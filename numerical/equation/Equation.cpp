#include <Equation.h>

Equation::Equation(const std::vector<double> &propsVector) :
        props(propsVector),
        local(propsVector),
        convective(propsVector),
        dim(props.gridBlockN),
        time(props.time),
        iCurr(0),
        iPrev(1),
        matrix(dim, dim),
        freeVector(dim),
        guessVector(dim),
        variable(dim) {


    conc.emplace_back(std::vector<double>(dim, 0));
    conc.emplace_back(std::vector<double>(dim, 0));

    std::vector<Triplet> triplets;
    triplets.reserve(3 * dim - 4);
    triplets.emplace_back(0, 0, 1);
    for (int i = 1; i < dim - 1; i++) {
        triplets.emplace_back(i, i - 1);
        triplets.emplace_back(i, i);
        triplets.emplace_back(i, i + 1);
    }
    triplets.emplace_back(dim - 1, dim - 1, 1);
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
    MatrixIterator(matrix, dim - 1).valueRef() = local.alpha[dim - 1];
}

void Equation::calculateFreeVector(const double &conc_in,
                                   const double &conc_out) {
    freeVector[0] = local.alpha[0] * conc_in;
    for (int i = 1; i < dim - 1; i++)
        freeVector[i] = local.alpha[i] * conc[iPrev][i];
    freeVector[dim - 1] = local.alpha[dim - 1] * conc_out;
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

void Equation::processConc() {

    for (int i = 1; i < 2; i++)
        for (int j = 0; j < dim; j++)

            conc_vec.emplace_back(conc[i][j]);
}

void Equation::cfdProcedure(const double &concIn,
                            const double &concOut) {

    for (double t = props.timeStep; t <= props.time; t += props.timeStep) {

        std::swap(iCurr, iPrev);
        local.calculateAlpha(props.timeStep,
                             props.radius,
                             props.effRadius,
                             props.gridBlockN);
        convective.calculateBeta(props.radius,
                                 props.effRadius,
                                 props.length,
                                 props.diffusivity,
                                 props.gridBlockN);
        calculateGuessVector();
        calculateMatrix();
        calculateFreeVector(concIn, concOut);
        calculateConc();
    }

    processConc();
}


const std::vector<double> Equation::getConc() const {

    return conc_vec;
}

