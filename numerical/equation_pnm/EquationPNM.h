#ifndef EQUATIONPNM_H
#define EQUATIONPNM_H

#include <vector>
#include <string>
#include <variant>
#include <map>

#include <Eigen/Sparse>
#include <NetworkData.h>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Eigen::VectorXd Vector;

class EquationPNM {

public:

    explicit EquationPNM(
            const std::map<std::string, std::variant<int, double>> &paramsPnm,
            NetworkData &networkData,
            const std::string &solverMethod);

    explicit EquationPNM(
            const std::map<std::string, std::variant<int, double>> &paramsPnm,
            const std::map<std::string, std::variant<std::vector<bool>,
                    std::vector<int>, std::vector<double>>> &paramsNetwork,
            const std::string &solverMethod);

    virtual ~EquationPNM() = default;

    void calculateMatrix(const std::vector<double> &connCoeff,
                         const std::vector<double> &centralCoeff,
                         const std::vector<bool> &boundPores,
                         std::vector<std::vector<int>> &inOutCoeff,
                         const std::vector<double> &diffCoeff);

    void calcMatCoeff();

    void calculateFreeVector(const std::string &boundCond,
                             const double &pIn,
                             const double &pOut);

    void calculateGuessPress(const double &pIn,
                             const double &pOut);

    void calculateGuessVector();

    void calculatePress(const std::string &solverMethod);

    void cfdProcedure(const std::string &boundCond,
                      const std::vector<bool> &boundPores,
                      const double &pIn,
                      const double &pOut);

    void cfdProcPurePnmDirichlet();

    void calcThrFlowRate();

    void calcPorFlowRate();

    void getGammaByLabel();

    void getGammaByPressure();

    void setInitialCondPurePnm();

    void calcInletFlow(const std::vector<bool> &boundPorIn);

    void calcTotFlow(const std::vector<bool> &boundPores);

    std::map<std::string, std::variant<int, double>> _paramsPnm;
    NetworkData &networkData;

    int dim;

    Matrix matrix;

    Vector freeVector;
    Vector guessVector;
    Vector variable;

    std::vector<std::vector<bool>> porConnsIsOut;
    std::vector<std::vector<int>> gammaPnm;

    std::vector<double> connCoeff;
    std::vector<double> centralCoeff;
    std::vector<double> pressure;

    std::vector<double> thrFlowRate;
    std::vector<double> porFlowRate;

    double totFlowRate;

    std::string solverMethod;

    std::vector<double> inletFlow;

};


#endif
