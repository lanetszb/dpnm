#include <CouplingIniConds.h>

CouplingIniConds::CouplingIniConds(const std::vector<double> &propsPNM,
                                   const std::vector<double> &propsDiffusion,
                                   const std::vector<int> &throatList,
                                   const std::vector<double> &throatHeight,
                                   const std::vector<double> &throatLength,
                                   const std::vector<double> &throatWidth,
                                   const std::vector<double> &connIndIn,
                                   const std::vector<double> &connIndOut,
                                   const std::vector<double> &poreCoordX,
                                   const std::vector<double> &poreCoordY,
                                   const std::vector<double> &poreCoordZ,
                                   const std::vector<double> &poreRadius,
                                   const std::vector<int> &poreList,
                                   const std::vector<int> &poreConns,
                                   const std::vector<int> &connNumber,
                                   const std::vector<int> &porePerRow,
                                   const std::vector<bool> &poreLeftX,
                                   const std::vector<bool> &poreRightX,
                                   const std::vector<double> &hydraulicCond,
                                   const std::vector<double> &langmuirCoeff,
                                   const double &matrixVolume) :

        diffusionPartMath(propsPNM, propsDiffusion, throatList, throatHeight,
                          throatLength, throatWidth, connIndIn, connIndOut,
                          poreCoordX, poreCoordY, poreCoordZ, poreRadius,
                          poreList, poreConns, connNumber, porePerRow,
                          poreLeftX, poreRightX, hydraulicCond, langmuirCoeff,
                          matrixVolume),

        equationPNM(propsPNM, throatList, throatHeight, throatLength,
                    throatWidth, connIndIn, connIndOut, poreCoordX, poreCoordY,
                    poreCoordZ, poreRadius, poreList, poreConns, connNumber,
                    porePerRow, poreLeftX, poreRightX, hydraulicCond),

        equationDiffusion(propsDiffusion) {}

std::vector<std::vector<int>> CouplingIniConds::getGamma() {

    // calculate PN no diffusion with Direchlet-Newman for finding Gamma
    std::vector<bool> boundPoresInputForGamma;

    for (int i = 0; i < equationPNM.networkData.poreN; i++) {
        boundPoresInputForGamma.emplace_back(
                !equationPNM.networkData.poreRightX[i]);
    }

    std::vector<bool> poreLeftXSaved;
    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        poreLeftXSaved.push_back(equationPNM.networkData.poreLeftX[i]);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.networkData.poreLeftX[i] = boundPoresInputForGamma[i];

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.porFlowRate[i] = 1.e-12;

    equationPNM.cfdProcedure("mixed",
                             equationPNM.networkData.poreRightX,
                             equationPNM.pIn,
                             equationPNM.pOut);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        equationPNM.networkData.poreLeftX[i] = poreLeftXSaved[i];

    std::vector<std::vector<int>> poreConnsIsOutByPressureSaved(
            equationPNM.networkData.poreN);

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        for (int j = 0; j < equationPNM.gammaPnm[i].size(); j++)
            poreConnsIsOutByPressureSaved[i].push_back(
                    equationPNM.gammaPnm[i][j]);

    return poreConnsIsOutByPressureSaved;
}

void CouplingIniConds::getInletFlow() {

    // calculate PN no diffusion with Direchlet
    std::vector<bool> boundPoresInput;

    for (int i = 0; i < equationPNM.networkData.poreN; i++)
        boundPoresInput.emplace_back(equationPNM.networkData.poreLeftX[i] or
                                     equationPNM.networkData.poreRightX[i]);

    equationPNM.cfdProcedure("dirichlet",
                             boundPoresInput,
                             equationPNM.pIn,
                             equationPNM.pOut);

    // calculate PN no diffusion with mixed Dirichlet-Newman

    equationPNM.cfdProcedure("mixed",
                             equationPNM.networkData.poreRightX,
                             equationPNM.pIn,
                             equationPNM.pOut);

}

void CouplingIniConds::setInitialCondCoupledMod() {

    diffusionPartMath.densityConst = diffusionPartMath.calcDensConst();

    // equationPNM.calcPorConns();
    equationPNM.calcThroatConns();

    diffusionPartMath.calcRockVolume();
    diffusionPartMath.calcMatrixWidth();
    diffusionPartMath.calcMatricesOmega();
    diffusionPartMath.calcMatricesVolume();

    // equationPNM.calculateGuessPress(equationPNM.propsPNM.pressIn,
    //                                equationPNM.propsPNM.pressOut);

    equationDiffusion.calcConcIni(diffusionPartMath.conc_ini);

    for (int i = 0; i < equationPNM.networkData.throatN; i++)
        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++)
            matrixConc[i][j] = diffusionPartMath.conc_ini;
}