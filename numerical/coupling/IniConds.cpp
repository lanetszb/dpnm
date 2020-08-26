#include <IniConds.h>
#include "NetworkData.h"

IniConds::IniConds(
        const std::map<std::string, std::variant<int, double>> &paramsPnm,
        NetworkData &networkData, EquationPNM &equationPNM,
        EquationDiffusion &equationDiffusion,
        DiffusionMath &diffusionMath,
        const std::vector<double> &langmuirCoeff,
        const double &matrixVolume) :

        _paramsPnm(paramsPnm),
        networkData(networkData),
        equationPNM(equationPNM),
        equationDiffusion(equationDiffusion),
        diffusionMath(diffusionMath),

        gridBlockN(equationDiffusion.propsDiffusion.gridBlockN),
        matrixConc(networkData.fracturesN,
                   std::vector<double>(gridBlockN, diffusionMath.conc_ini)) {}

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

    networkData.poreInlet = boundPoresInputForGamma;

    for (int i = 0; i < networkData.poreN; i++)
        equationPNM.porFlowRate[i] = 1.e-12;

    auto &pressIn = std::get<double>(_paramsPnm["pressIn"]);
    auto &pressOut = std::get<double>(_paramsPnm["pressOut"]);

    equationPNM.cfdProcedure("mixed",
                             networkData.poreOutlet,
                             pressIn,
                             pressOut);

    networkData.poreInlet = poreLeftXSaved;
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

    auto &pressIn = std::get<double>(_paramsPnm["pressIn"]);
    auto &pressOut = std::get<double>(_paramsPnm["pressOut"]);

    equationPNM.cfdProcedure("dirichlet",
                             boundPoresInput,
                             pressIn,
                             pressOut);

    // calculate PN no diffusion with mixed Dirichlet-Newman

    equationPNM.cfdProcedure("mixed",
                             networkData.poreOutlet,
                             pressIn,
                             pressOut);

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
}