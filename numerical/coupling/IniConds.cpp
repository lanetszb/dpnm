#include <IniConds.h>

IniConds::IniConds(EquationPNM &equationPNM,
                   EquationDiffusion &equationDiffusion,
                   DiffusionMath &diffusionMath,
                   const std::vector<double> &langmuirCoeff,
                   const double &matrixVolume) :

        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),
        matrixConc(equationPNM.networkData.throatN,
                   std::vector<double>(
                           equationDiffusion.propsDiffusion.gridBlockN, 0)) {}

std::vector<std::vector<int>> IniConds::getGamma() {

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

void IniConds::getInletFlow() {

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

void IniConds::setInitialCondCoupledMod() {

    auto gammaByPressureSaved = getGamma();

    getInletFlow();

    equationPNM.gammaPnm = gammaByPressureSaved;

    diffusionMath.densityConst = diffusionMath.calcDensConst();

    // equationPNM.calcPorConns();
    equationPNM.calcThroatConns();

    diffusionMath.calcRockVolume();
    diffusionMath.calcMatrixWidth();
    diffusionMath.calcMatricesOmega();
    diffusionMath.calcMatricesVolume();

    // equationPNM.calculateGuessPress(equationPNM.propsPnm.pressIn,
    //                                equationPNM.propsPnm.pressOut);

    equationDiffusion.calcConcIni(diffusionMath.conc_ini);

    for (int i = 0; i < equationPNM.networkData.throatN; i++) {
        for (int j = 0; j < equationDiffusion.propsDiffusion.gridBlockN; j++) {
            matrixConc[i][j] = diffusionMath.conc_ini;
        }
    }
}