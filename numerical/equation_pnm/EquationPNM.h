#ifndef EQUATIONPNM_H
#define EQUATIONPNM_H

#include <vector>
#include <string>

#include <Eigen/Sparse>

#include <PropsPnm.h>
#include <NetworkData.h>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Eigen::VectorXd Vector;

class EquationPNM {

public:

    explicit EquationPNM(PropsPnm &propsPnm,
                         NetworkData &networkData,
                         const std::string &solverMethod);

    explicit EquationPNM(const std::vector<double> &propsVector,
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

    PropsPnm &propsPnm;
    NetworkData &networkData;

    int &dim;
    double &pIn;
    double &pOut;

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


    // manual construction, has to be automatised later

    std::vector<double> inletFlow;

    const std::vector<double> getPressure() const;

    double getTotFlowRate() const;

};


#endif
