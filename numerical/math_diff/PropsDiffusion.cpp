#include <PropsDiffusion.h>
#include <vector>


PropsDiffusion::PropsDiffusion(
        const std::map<std::string, std::variant<int, double>> &params) :

        _params(params),
        time(std::get<double>(_params["time"])),
        timeStep(std::get<double>(_params["time_step"])),
        gridBlockN(std::get<int>(_params["grid_block_n"])),
        length(std::get<double>(_params["length"])),
        radius(std::get<double>(_params["radius"])),
        effRadius(std::get<double>(_params["eff_radius"])),
        concIni(std::get<double>(_params["conc_ini"])),
        diffusivity(std::get<double>(_params["diffusivity"])),
        iterativeAccuracy(std::get<double>(_params["it_accuracy"])) {}


void PropsDiffusion::printParams() {
    for (auto &ent : _params) {
        std::cout << ent.first << ": ";
        if (std::get_if<int>(&ent.second))
            std::cout << std::get<int>(ent.second) << std::endl;
        else if (std::get_if<double>(&ent.second))
            std::cout << std::get<double>(ent.second) << std::endl;
    }
}


