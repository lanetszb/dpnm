#include <boost/python.hpp>

#include <string>

#include <PropsDiffusion.h>
#include <LocalDiffusion.h>
#include <ConvectiveDiffusion.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>
#include <PropsPnm.h>
#include <Aggregator.h>


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
            .add_property("conc_ini",
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
                         "throat_length"))

            .def("calc_vol_cartesian",
                 &LocalDiffusion::calcVolCartesian,
                 p::args("frac_height",
                         "matrix_width",
                         "frac_length",
                         "frac_width"))

            .def("calc_matr_coord_curr",
                 &LocalDiffusion::calcMatrCoordCurr,
                 p::args("radius",
                         "eff_radius"))

            .add_property("radius_curr",
                          &LocalDiffusion::getRadCurr)

            .add_property("vol_cylindr",
                          &LocalDiffusion::getVolCylindr)

            .add_property("vol_cartes",
                          &LocalDiffusion::getVolCartes);


    p::class_<ConvectiveDiffusion>("ConvectiveDiffusion",
                                   p::init<std::vector<double>>(
                                           p::args("props_array")))

            .def("calc_omega_cylindr",
                 &ConvectiveDiffusion::calcOmegaCylindr,
                 p::args("length"))

            .def("calc_omega_cartes",
                 &ConvectiveDiffusion::calcOmegaCartes,
                 p::args("frac_height",
                         "frac_length"))

            .add_property("omega_cylindr",
                          &ConvectiveDiffusion::getOmegaCylindr)

            .add_property("omega_cartes",
                          &ConvectiveDiffusion::getOmegaCartes);
//

    p::class_<EquationDiffusion>("EquationDiffusion",
                                 p::init<std::vector<double>>(
                                         p::args("props_array")))

            .def("cfd_procedure_one_step",
                 &EquationDiffusion::cfdProcedureOneStep,
                 p::args("boundCond",
                         "conc_in",
                         "radius",
                         "effective_radius",
                         "thr_length",
                         "volumes",
                         "surfaces"))

            .def("cfd_procedure",
                 &EquationDiffusion::cfdProcedure,
                 p::args("boundCond",
                         "volumes",
                         "surfaces"))

            .def("getConc",
                 &EquationDiffusion::getConc)

            .def("getFlowRate",
                 &EquationDiffusion::getFlowRate);
//
//    // Wrapper for PNM
//
    p::class_<PropsPnm>("PropsPNMCpp",
                        p::init<std::vector<double>>(
                                p::args("props_array")))
            .def("__str__", __str__<PropsPnm>)

            .add_property("a_gas_dens",
                          &PropsPnm::getAGasDens,
                          &PropsPnm::setAGasDens)
            .add_property("b_gas_dens",
                          &PropsPnm::getBGasDens,
                          &PropsPnm::setBGasDens)
            .add_property("gas_visc",
                          &PropsPnm::getGasVisc,
                          &PropsPnm::setGasVisc)
            .add_property("liq_dens",
                          &PropsPnm::getLiqDens,
                          &PropsPnm::setLiqDens)
            .add_property("liq_visc",
                          &PropsPnm::getLiqVisc,
                          &PropsPnm::setLiqVisc)
            .add_property("props_vector",
                          &PropsPnm::getPropsVector)

            .def("print_props_vector",
                 &PropsPnm::printPropsVector);

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
                                   std::vector<double>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<int>,
                                   std::vector<bool>,
                                   std::vector<bool>,
                                   std::vector<double>,
                                   std::string>(
                                   p::args("props_array",
                                           "throat_list",
                                           "throat_height",
                                           "throat_length",
                                           "throat_width",
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
                                           "hydr_cond",
                                           "solver_method")))

            .def("cfd_proc_pure_pnm_dirichlet",
                 &EquationPNM::cfdProcPurePnmDirichlet)

            .add_property("tot_flow_rate",
                          &EquationPNM::getTotFlowRate)

            .add_property("pressure",
                          &EquationPNM::getPressure);
//            .def("__str__", __str__<EquationPNM>);

    p::class_<Aggregator>("Aggregator",
                          p::init<std::vector<double>,
                                  std::vector<double>,
                                  std::vector<int>,
                                  std::vector<double>,
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
                                  std::vector<double>,
                                  std::vector<double>,
                                  double,
                                  std::string>(
                                  p::args("props_PNM",
                                          "props_diffusion",
                                          "throat_list",
                                          "throat_height",
                                          "throat_length",
                                          "throat_width",
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
                                          "hydr_cond",
                                          "langmuir_coeff",
                                          "matrix_volume",
                                          "solver_method")))

            .def("cfd_procedure_pnm_diff",
                 &Aggregator::cfdProcedurePnmDiff)

            .def("get_pressure_av",
                 &Aggregator::getPressureAverage)

            .def("get_matrix_mass_total",
                 &Aggregator::getMatrixMassTotal)

            .def("get_flow_pores_out",
                 &Aggregator::getTotalFlowPoresOut)

            .def("get_flow_pores_in",
                 &Aggregator::getTotalFlowPoresIn)

            .def("get_flow_diff",
                 &Aggregator::getTotalFlowDiff)

            .def("get_pressure_by_pore",
                 &Aggregator::getPorePressure)

            .def("get_inlet_pressure",
                 &Aggregator::getInletPressure)

            .def("get_time_steps_vec",
                 &Aggregator::getTimeStepsVec);
}




