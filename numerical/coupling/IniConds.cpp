#include <IniConds.h>
#include "NetworkData.h"

IniConds::IniConds(NetworkData &networkData, EquationPNM &equationPNM,
                   EquationDiffusion &equationDiffusion,
                   DiffusionMath &diffusionMath,
                   const std::vector<double> &langmuirCoeff,
                   const double &matrixVolume) :

        networkData(networkData),
        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),

        gridBlockN(equationDiffusion.propsDiffusion.gridBlockN),
        matrixConc(networkData.fracturesN,
                   std::vector<double>(gridBlockN, 0)) {}

std::vector<std::vector<int>> IniConds::getGamma() {

    // calculate PN no diffusion with Direchlet-Newman for finding Gamma
    std::vector<bool> boundPoresInputForGamma;

    for (int i = 0; i < networkData.poreN; i++) {
        boundPoresInputForGamma.emplace_back(
                !equationPNM.networkData.poreOutlet[i]);
    }

    std::vector<bool> poreLeftXSaved;
    for (int i = 0; i < networkData.poreN; i++)
        poreLeftXSaved.push_back(networkData.poreInlet[i]);

    for (int i = 0; i < networkData.poreN; i++)
        networkData.poreInlet[i] = boundPoresInputForGamma[i];

    for (int i = 0; i < networkData.poreN; i++)
        equationPNM.porFlowRate[i] = 1.e-12;

    equationPNM.cfdProcedure("mixed",
                             networkData.poreOutlet,
                             equationPNM.pIn,
                             equationPNM.pOut);

    for (int i = 0; i < networkData.poreN; i++)
        networkData.poreInlet[i] = poreLeftXSaved[i];

    std::vector<std::vector<int>> poreConnsIsOutByPressureSaved(
            networkData.poreN);

    for (int i = 0; i < networkData.poreN; i++)
        for (int j = 0; j < equationPNM.gammaPnm[i].size(); j++)
            poreConnsIsOutByPressureSaved[i].push_back(
                    equationPNM.gammaPnm[i][j]);

    return poreConnsIsOutByPressureSaved;
}

void IniConds::getInletFlow() {

    // calculate PN no diffusion with Direchlet
    std::vector<bool> boundPoresInput;

    for (int i = 0; i < networkData.poreN; i++)
        boundPoresInput.emplace_back(networkData.poreInlet[i] or
                                     networkData.poreOutlet[i]);

    equationPNM.cfdProcedure("dirichlet",
                             boundPoresInput,
                             equationPNM.pIn,
                             equationPNM.pOut);

    // calculate PN no diffusion with mixed Dirichlet-Newman

    equationPNM.cfdProcedure("mixed",
                             networkData.poreOutlet,
                             equationPNM.pIn,
                             equationPNM.pOut);

}

void IniConds::setInitialCondCoupledMod() {

    auto gammaByPressureSaved = getGamma();

    getInletFlow();

    equationPNM.gammaPnm = gammaByPressureSaved;

    diffusionMath.densityConst = diffusionMath.calcDensConst();

    // equationPNM.calcPor2ThrConns();
    networkData.calcThroatConns();

    diffusionMath.calcRockVolume();
    diffusionMath.calcMatrixWidth();
    diffusionMath.calcMatricesOmega();
    diffusionMath.calcMatricesVolume();

    // equationPNM.calculateGuessPress(equationPNM.propsPnm.pressIn,
    //                                equationPNM.propsPnm.pressOut);

    equationDiffusion.calcConcIni(diffusionMath.conc_ini);

    for (int i = 0; i < networkData.fracturesN; i++) {
        for (int j = 0; j < gridBlockN; j++) {
            matrixConc[i][j] = diffusionMath.conc_ini;
        }
    }
}