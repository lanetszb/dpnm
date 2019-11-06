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
            const std::vector<double> &propsVector,
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
            const std::vector<int> &pore_per_row);

    virtual ~EquationPNM() = default;

    void calculateMatrix(const int &boundCond,
                         const std::vector<double> &connCoeff,
                         const std::vector<double> &centralCoeff,
                         const std::vector<int> &boundPores,
                         std::vector<std::vector<int>> &inOutCoeff,
                         const std::vector<double> &diffCoeff);

    void calcThroatConns();

    void calcPorConns();

    void calcMatCoeff();

    void calculateFreeVector(const int &boundCond,
                             const double &pIn,
                             const double &pOut);

    void calculateGuessPress(const double &pIn,
                             const double &pOut);

    void calculateGuessVector();

    void calculatePress(const int &solverMethod);

    void cfdProcedure(const int &boundCond,
                      const std::vector<int> &boundPores,
                      const double &pIn,
                      const double &pOut);

    void calcThrFlowRate();

    void calcPorFlowRate();

    void getPorConnsIsOut();

    void getPorConnsIsOutByPressure();

    void setInitialCond();

    void calcInletFlow(const int &boundPorSize);

    void calcTotFlow(const std::vector<int> &boundPores);


    PropsPNM propsPNM;
    NetworkData networkData;

    int &dim;
    double pIn;
    double pOut;

    Matrix matrix;

    Vector freeVector;
    Vector guessVector;
    Vector variable;

    std::vector<std::pair<int, int>> throatConns;
    std::vector<std::vector<int>> porConns;

    std::vector<std::vector<bool>> porConnsIsOut;
    std::vector<std::vector<int>> porConnsIsOutByPressure;

    std::vector<double> connCoeff;
    std::vector<double> centralCoeff;
    std::vector<double> pressure;

    std::vector<double> thrFlowRate;
    std::vector<double> porFlowRate;

    double totFlowRate;


    // manual construction, has to be automatised later

    std::vector<double> inletFlow;

};


#endif
