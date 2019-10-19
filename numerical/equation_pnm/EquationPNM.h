#ifndef PNFLOW_EQUATIONPNM_H
#define PNFLOW_EQUATIONPNM_H

#include <vector>
#include <Eigen/Sparse>

#include <PropsPNM.h>
#include <NetworkData.h>

typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Matrix::InnerIterator MatrixIterator;
typedef Eigen::VectorXd Vector;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> BiCGSTAB;
typedef Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>
        LeastSqCG;
typedef Eigen::SparseLU<Eigen::SparseMatrix<double>> SparseLU;

class EquationPNM {

public:

    explicit EquationPNM(
            const std::vector<double> &_propsVector,
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
            const std::vector<int> &_pore_per_row);

    virtual ~EquationPNM() = default;

    void calculateMatrix();

    void calcThroatConns();

    void calcPorConns();

    void calcMatCoeff();

    void calculateFreeVector(const double &pIn,
                             const double &pOut);

    void calculateGuessPress(const double &pIn,
                             const double &pOut);

    void calculateGuessVector();

    void calculatePress();

    void cfdProcedure(const double &pIn,
                      const double &pOut);


    std::vector<std::vector<double>> press;

    PropsPNM propsPnm;
    NetworkData networkData;

    int &dim;
    double pIn;
    double pOut;
    int iCurr;
    int iPrev;

    Matrix matrix;

    Vector freeVector;
    Vector guessVector;
    Vector variable;

    std::vector<std::pair<int, int>> throatConns;
    std::vector<std::vector<int>> porConns;
    std::vector<double> connCoeff;
    std::vector<double> centralCoeff;
    std::vector<double> pressure;

};


#endif
