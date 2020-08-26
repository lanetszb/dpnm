#ifndef PNFLOW_EQUATIONDIFFUSION_H
#define PNFLOW_EQUATIONDIFFUSION_H

#include <vector>
#include <Eigen/Sparse>

#include <PropsDiffusion.h>
#include <LocalDiffusion.h>
#include <ConvectiveDiffusion.h>


typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix;
typedef Matrix::InnerIterator MatrixIterator;
typedef Eigen::VectorXd Vector;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> BiCGSTAB;


class EquationDiffusion {

public:

    explicit EquationDiffusion(PropsDiffusion &propsDiffusion);

    explicit EquationDiffusion(
            const std::map<std::string, std::variant<int, double>> &params);

    virtual ~EquationDiffusion() = default;

    void calculateMatrix();

    void calculateFreeVector(const double &conc_in);

    void calculateGuessVector();

    void calculateConc();

    void calcConcIni(const double &concIni);

    void calcTimeVector();

    void forceDirichletBound(const double &concIni);

    void cfdProcedureOneStep(const std::string &boundCond,
                             const double &concThrWall,
                             const double &radius,
                             const double &effRadius,
                             const double &thrLength,
                             const std::vector<double> &volumes,
                             const std::vector<double> &surfaces,
                             const double &dt);

    void cfdProcedure(const std::string &boundCond,
                      const std::vector<double> &volumes,
                      const std::vector<double> &surfaces);


    void calcFlowRate();

    void cfdCartesian(const std::string &boundCond);

    int &dim;

    std::vector<std::vector<double>> conc;
    std::vector<double> timeStepsVec;

    PropsDiffusion &propsDiffusion;
    LocalDiffusion localDiffusion;
    ConvectiveDiffusion convectiveDiffusion;

    double &time;
    int iCurr;
    int iPrev;

    double flowRate;

    Matrix matrix;

    Vector freeVector;
    Vector guessVector;
    Vector variable;

};


#endif
