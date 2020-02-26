#include <boost/python.hpp>

#include <string>

#include <PropsDiffusion.h>
#include <LocalDiffusion.h>
#include <ConvectiveDiffusion.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>
#include <PropsPNM.h>
#include <NetworkData.h>
#include <DiffusionPNM.h>


#include <PythonConversions.h>

namespace p = boost::python;

BOOST_PYTHON_MODULE (cfd) {

    providePythonThings();

    p::class_<PropsDiffusion>("PropsDiffusion",
                              p::init<std::vector<double>>(
                                      p::args("props_array")))
            .def("__str__", __str__<PropsDiffusion>)

            .add_property("time",
                          &PropsDiffusion::getTime,
                          &PropsDiffusion::setTime)
            .add_property("time_step",
                          &PropsDiffusion::getTimeStep,
                          &PropsDiffusion::setTimeStep)
            .add_property("length",
                          &PropsDiffusion::getLength,
                          &PropsDiffusion::setLength)
            .add_property("radius",
                          &PropsDiffusion::getRadius,
                          &PropsDiffusion::setRadius)
            .add_property("eff_radius",
                          &PropsDiffusion::getEffRadius,
                          &PropsDiffusion::setEffRadius)
            .add_property("grid_block_n",
                          &PropsDiffusion::getGridBlockN,
                          &PropsDiffusion::setGridBlockN)
            .add_property("concIni",
                          &PropsDiffusion::getConcentration,
                          &PropsDiffusion::setConcentration)
            .add_property("diffusivity",
                          &PropsDiffusion::getDiffusivity,
                          &PropsDiffusion::setDiffusivity)
            .add_property("iterative_accuracy",
                          &PropsDiffusion::getIterativeAccuracy,
                          &PropsDiffusion::setIterativeAccuracy)
            .add_property("props_vector",
                          &PropsDiffusion::getPropsVector)

            .def("print_props_vector",
                 &PropsDiffusion::printPropsVector);


    p::class_<LocalDiffusion>("LocalDiffusion",
                              p::init<std::vector<double>>(
                                      p::args("props_array")))

            .def("calc_vol_cylinder",
                 &LocalDiffusion::calcVolCylindr,
                 p::args("radius",
                         "eff_radius",
                         "grid_block_n",
                         "throat_length"))

            .def("calc_vol_cartesian",
                 &LocalDiffusion::calcVolCartesian,
                 p::args("frac_height",
                         "matrix_width",
                         "frac_length",
                         "frac_width"))

            .def("calculate_alpha",
                 &LocalDiffusion::calculateAlpha,
                 p::args("dt",
                         "volume_list"))

            .add_property("radius_curr",
                          &LocalDiffusion::getRadCurr)

            .add_property("vol_cylindr",
                          &LocalDiffusion::getVolCylindr)

            .add_property("vol_cartes",
                          &LocalDiffusion::getVolCartes)

            .add_property("alpha",
                          &LocalDiffusion::getAlpha);



//
    p::class_<ConvectiveDiffusion>("ConvectiveDiffusion",
                                   p::init<std::vector<double>>(
                                           p::args("props_array")))

            .def("calc_diffusivityList",
                 &ConvectiveDiffusion::calc_diffusivityList,
                 p::args("grid_block_n",
                         "concIni"))

            .def("calc_omega_cylindr",
                 &ConvectiveDiffusion::calcOmegaCylindr,
                 p::args("length"))

            .def("calc_omega_cartes",
                 &ConvectiveDiffusion::calcOmegaCartes,
                 p::args("frac_height",
                         "frac_length"))

            .def("calculate_beta",
                 &ConvectiveDiffusion::calculateBeta,
                 p::args("radius",
                         "effRadius",
                         "length",
                         "diffusivity",
                         "grid_block_n",
                         "omega"))

            .def("get_beta",
                 &ConvectiveDiffusion::getBeta)

            .add_property("omega_cylindr",
                          &ConvectiveDiffusion::getOmegaCylindr)

            .add_property("omega_cartes",
                          &ConvectiveDiffusion::getOmegaCartes);
//
    p::class_<EquationDiffusion>("EquationDiffusion",
                                 p::init<std::vector<double>>(
                                         p::args("props_array")))

            .def("cfd_procedure",
                 &EquationDiffusion::cfdProcedure,
                 p::args("boundCond",
                         "conc_in",
                         "radius",
                         "effective_radius",
                         "thr_length",
                         "volumes",
                         "surfaces"))

            .def("getConc",
                 &EquationDiffusion::getConc)

            .def("getFlowRate",
                 &EquationDiffusion::getFlowRate);
//
//    // Wrapper for PNM
//
    p::class_<PropsPNM>("PropsPNMCpp",
                        p::init<std::vector<double>>(
                                p::args("props_array")))
            .def("__str__", __str__<PropsPNM>)

            .add_property("a_gas_dens",
                          &PropsPNM::getAGasDens,
                          &PropsPNM::setAGasDens)
            .add_property("b_gas_dens",
                          &PropsPNM::getBGasDens,
                          &PropsPNM::setBGasDens)
            .add_property("gas_visc",
                          &PropsPNM::getGasVisc,
                          &PropsPNM::setGasVisc)
            .add_property("liq_dens",
                          &PropsPNM::getLiqDens,
                          &PropsPNM::setLiqDens)
            .add_property("liq_visc",
                          &PropsPNM::getLiqVisc,
                          &PropsPNM::setLiqVisc)
            .add_property("props_vector",
                          &PropsPNM::getPropsVector)

            .def("print_props_vector",
                 &PropsPNM::printPropsVector);
//
//
    p::class_<NetworkData>("NetworkDataCpp",
                           p::init<std::vector<int>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<bool>,
                                   std::vector<bool>,
                                   std::vector<double>>(
                                   p::args("throat_list",
                                           "throat_radius",
                                           "throat_length",
                                           "conn_ind_in",
                                           "conn_ind_out",
                                           "pore_coord_x",
                                           "pore_coord_y",
                                           "pore_coord_z",
                                           "pore_radius",
                                           "pore_list",
                                           "pore_conns",
                                           "conn_number",
                                           "pore_per_row",
                                           "pore_left_x",
                                           "pore_right_x",
                                           "hydr_cond")))
            .def("__str__", __str__<NetworkData>)

            .add_property("throat_radius",
                          &NetworkData::getThroatRadius)
            .add_property("throat_length",
                          &NetworkData::getThroatLength)
            .add_property("pore_list",
                          &NetworkData::getPoreList)

            .def("print_throat_radius",
                 &NetworkData::printThroatRadius)
            .def("print_throat_radius",
                 &NetworkData::printThroatLength)
            .def("print_throat_radius",
                 &NetworkData::printPoreList);
//
    p::class_<EquationPNM>("EquationPNM",
                           p::init<std::vector<double>,
                                   std::vector<int>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<double>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<bool>,
                                   std::vector<bool>,
                                   std::vector<double>>(
                                   p::args(
                                           "props_array",
                                           "throat_list",
                                           "throat_radius",
                                           "throat_length",
                                           "conn_ind_in",
                                           "conn_ind_out",
                                           "pore_coord_x",
                                           "pore_coord_y",
                                           "pore_coord_z",
                                           "pore_radius",
                                           "pore_list",
                                           "pore_conns",
                                           "conn_number",
                                           "pore_per_row",
                                           "pore_left_x",
                                           "pore_right_x",
                                           "hydr_cond")));
//            .def("__str__", __str__<EquationPNM>);
}

//    p::class_<DiffusionPNM>("DiffusionPNM",
//                            p::init<std::vector<double>,
//                                    std::vector<int>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<int>,
//                                    std::vector<int>,
//                                    std::vector<int>,
//                                    std::vector<int>,
//                                    std::vector<double>,
//                                    std::vector<double>,
//                                    std::vector<bool>,
//                                    std::vector<bool>,
//                                    std::vector<double>>(
//                                    p::args("props_vector",
//                                            "throat_list",
//                                            "throat_radius",
//                                            "throat_length",
//                                            "conn_ind_in",
//                                            "conn_ind_out",
//                                            "pore_coord_x",
//                                            "pore_coord_y",
//                                            "pore_coord_z",
//                                            "pore_radius",
//                                            "pore_list",
//                                            "pore_conns",
//                                            "conn_number",
//                                            "pore_per_row",
//                                            "props_diffusion",
//                                            "langmuir_coeff",
//                                            "pore_left_x",
//                                            "pore_right_x",
//                                            "hydr_cond")))
//
//            .def("get_pressure_av",
//                 &DiffusionPNM::getPressureAverage)
//
//            .def("get_conc_av",
//                 &DiffusionPNM::getConcAverage)
//
//            .def("get_flow_pores_out",
//                 &DiffusionPNM::getTotalFlowPoresOut)
//
//            .def("get_flow_pores_in",
//                 &DiffusionPNM::getTotalFlowPoresIn)
//
//            .def("get_flow_diff",
//                 &DiffusionPNM::getTotalFlowDiff)
//
//            .def("get_pressure_by_pore",
//                 &DiffusionPNM::getPorePressure);
//
//}




