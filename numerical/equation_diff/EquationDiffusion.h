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

    explicit EquationDiffusion(const std::vector<double> &propsVector);

    virtual ~EquationDiffusion() = default;

    void calculateMatrix();

    void calculateFreeVector(const double &conc_in);

    void calculateGuessVector();

    void calculateConc();

    void calcConcIni(const double &concIni);

    void forceDirichletBound(const double &concIni);

    // ToDo DRY is not followed: radius, effRadius, thrLength, volumes, surfaces

    void cfdProcedureOneStep(const double &concThrWall,
                             const double &radius,
                             const double &effRadius,
                             const double &thrLength,
                             const std::vector<double> &volumes,
                             const std::vector<double> &surfaces);

    void cfdProcedure(const bool &boundCond,
                      const double &concThrWall,
                      const double &radius,
                      const double &effRadius,
                      const double &thrLength,
                      const std::vector<double> &volumes,
                      const std::vector<double> &surfaces);

    void cfdProcDirichlet(const double &concThrWall,
                          const double &radius,
                          const double &effRadius,
                          const double &thrLength,
                          const std::vector<double> &volumes,
                          const std::vector<double> &surfaces);

    void calcFlowRate();

    int &dim;

    std::vector<std::vector<double>> conc;

    const std::vector<double> getConc() const;

    const double getFlowRate() const;

    PropsDiffusion propsDiffusion;
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
