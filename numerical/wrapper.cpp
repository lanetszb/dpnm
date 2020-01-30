#include <boost/python.hpp>

#include <string>

#include <PropsDiffusion.h>
#include <Local.h>
#include <Convective.h>
#include <Equation.h>
#include <EquationPNM.h>
#include <PropsPNM.h>
#include <NetworkData.h>
#include <DiffusionPNM.h>


#include <PythonConversions.h>

namespace p = boost::python;

BOOST_PYTHON_MODULE (cfd) {

    providePythonThings();

    p::class_<PropsDiffusion>("PropsDiffusion",
                              p::init<std::vector<double>,
                             std::vector<double>>(
                             p::args("props_array",
                                     "langm_coeff")))
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
            .add_property("langmuir_coeff",
                          &PropsDiffusion::getLangmuirCoeff)

            .def("print_props_vector",
                 &PropsDiffusion::printPropsVector)
            .def("print_langmuir_coeffs",
                 &PropsDiffusion::printLangmuirCoeff);


    p::class_<Local>("LocalCpp",
                     p::init<std::vector<double>,
                             std::vector<double>>(
                             p::args("props_array",
                                     "langmuir_coeff")))

            .def("calculate_alpha",
                 &Local::calculateAlpha,
                 p::args("dt",
                         "radius",
                         "effRadius",
                         "throat_length"))

            .add_property("alpha",
                          &Local::getAlpha)

            .add_property("radius_curr",
                          &Local::getRadCurr);
//
//
    p::class_<Convective>("ConvectiveCpp",
                          p::init<std::vector<double>,
                                  std::vector<double>>(
                                  p::args("props_array",
                                          "langmuir_coeff")))


            .def("calc_diffusivityList",
                 &Convective::calc_diffusivityList,
                 p::args("grid_block_n",
                         "concIni"))

            .def("calculateBeta",
                 &Convective::calculateBeta,
                 p::args("radius",
                         "effRadius",
                         "length",
                         "diffusivity"))

            .def("get_beta",
                 &Convective::getBeta);
//
    p::class_<Equation>("EquationCpp",
                        p::init<std::vector<double>,
                                std::vector<double>>(
                                p::args("props_array",
                                        "langmuir_coeff")))

            .def("cfdProcedure",
                 &Equation::cfdProcedure,
                 p::args("conc_in",
                         "radius",
                         "effective_radius",
                         "thr_length"))

            .def("getConc",
                 &Equation::getConc)

            .def("getFlowRate",
                 &Equation::getFlowRate);
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

    p::class_<DiffusionPNM>("DiffusionPNM",
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
                                    std::vector<double>,
                                    std::vector<double>,
                                    std::vector<bool>,
                                    std::vector<bool>,
                                    std::vector<double>>(
                                    p::args("props_vector",
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
                                            "props_diffusion",
                                            "langmuir_coeff",
                                            "pore_left_x",
                                            "pore_right_x",
                                            "hydr_cond")))

            .def("get_pressure_av",
                 &DiffusionPNM::getPressureAverage)

            .def("get_conc_av",
                 &DiffusionPNM::getConcAverage)

            .def("get_flow_pores_out",
                 &DiffusionPNM::getTotalFlowPoresOut)

            .def("get_flow_pores_in",
                 &DiffusionPNM::getTotalFlowPoresIn)

            .def("get_flow_diff",
                 &DiffusionPNM::getTotalFlowDiff)

            .def("get_pressure_by_pore",
                 &DiffusionPNM::getPorePressure);

}




