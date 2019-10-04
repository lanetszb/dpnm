#ifndef PNFLOW_EQUATION_H
#define PNFLOW_EQUATION_H

#include <vector>
#include <Eigen/Sparse>

#include <Props.h>
#include <Local.h>
#include <Convective.h>


typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Matrix::InnerIterator MatrixIterator;
typedef Eigen::VectorXd Vector;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> BiCGSTAB;


class Equation {

public:

    explicit Equation(const std::vector<double> &propsVector);

    virtual ~Equation() = default;

    void calculateMatrix();

    void calculateFreeVector(const double &conc_in);

    void calculateGuessVector();

    void calculateConc();

    void cfdProcedure(const double &concIn);

    int &dim;

    std::vector<std::vector<double>> conc;

    const std::vector<double> getConc() const;

    Props props;

    Local local;

    Convective convective;


    double &time;

    int iCurr;
    int iPrev;

    Matrix matrix;

    Vector freeVector;

    Vector guessVector;

    Vector variable;


};


#endif
