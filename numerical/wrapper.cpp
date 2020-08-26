#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include <PropsDiffusion.h>
#include <EquationDiffusion.h>
#include <EquationPNM.h>
#include <Aggregator.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(dpnm, m) {

    py::class_<PropsDiffusion>(m, "PropsDiffusion")
            .def(py::init<const std::map<std::string, std::variant<int, double>> &>(),
                 "params"_a)

            .def_readwrite("params", &PropsDiffusion::_params)
            .def("print_params", &PropsDiffusion::printParams);


    py::class_<EquationDiffusion>(m, "EquationDiffusion")
            .def(py::init<const std::map<std::string, std::variant<int, double>> &>(),
                 "params"_a)

            .def("cfd_cartesian", &EquationDiffusion::cfdCartesian,
                 "bound_cond"_a)

            .def_readwrite("flow_rate", &EquationDiffusion::flowRate)
            .def_readwrite("conc", &EquationDiffusion::conc);

// Wrapper for PNM

    py::class_<EquationPNM>(m, "EquationPNM")
            .def(py::init<const std::map<std::string, std::variant<int, double>> &,
                         const std::map<std::string, std::variant<std::vector<bool>,
                                 std::vector<int>, std::vector<double>>> &, const std::string &>(),
                 "params_pnm"_a, "params_network"_a, "solver_method"_a)

            .def("cfd_pnm_dirichlet", &EquationPNM::cfdProcPurePnmDirichlet)

            .def_readwrite("flow_rate", &EquationPNM::totFlowRate)
            .def_readwrite("pressure", &EquationPNM::pressure);


    py::class_<Aggregator>(m, "Aggregator")
            .def(py::init<const double &, const std::string &, const std::vector<double> &,
                         const std::map<std::string, std::variant<int, double>> &,
                         const std::map<std::string, std::variant<int, double>> &,
                         const std::map<std::string, std::variant<std::vector<bool>,
                                 std::vector<int>, std::vector<double>>> &>(),
                 "matrixVolume"_a, "solver_method"_a, "langmuirCoeffs"_a,
                 "params"_a, "params_pnm"_a, "params_network"_a)

            .def("cfd_procedure_pnm_diff", &Aggregator::cfdProcedurePnmDiff)

            .def_property("get_flow_pores_in",
                          &Aggregator::getTotalFlowPoresIn,
                          &Aggregator::getTotalFlowPoresIn)

            .def_property("get_pressure_av",
                          &Aggregator::getPressureAverage,
                          &Aggregator::getPressureAverage)

            .def_property("get_matrix_mass_total",
                          &Aggregator::getMatrixMassTotal,
                          &Aggregator::getMatrixMassTotal)

            .def_property("get_inlet_pressure",
                          &Aggregator::getInletPressure,
                          &Aggregator::getInletPressure)

            .def_property("get_time_steps_vec",
                          &Aggregator::getTimeStepsVec,
                          &Aggregator::getTimeStepsVec)

            .def_property("get_flow_diff",
                          &Aggregator::getTotalFlowDiff,
                          &Aggregator::getTotalFlowDiff);


}
//            .def("__str__", __str__<EquationPNM>);

//     p::class_<Aggregator>("Aggregator",
//                           p::init<std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<int>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   std::vector<int>,
//                                   std::vector<bool>,
//                                   std::vector<bool>,
//                                   std::vector<double>,
//                                   std::vector<double>,
//                                   double,
//                                   std::string>(
//                                   p::args("props_PNM",
//                                           "props_diffusion",
//                                           "fractures_list",
//                                           "fractures_heights",
//                                           "fractures_lengths",
//                                           "fractures_widths",
//                                           "fracs_conn_ind_in",
//                                           "fracs_conn_ind_out",
//                                           "pores_coords_x",
//                                           "pores_coords_y",
//                                           "pores_coords_z",
//                                           "pores_radii",
//                                           "pores_list",
//                                           "pores_inlet",
//                                           "pores_outlet",
//                                           "hydr_cond",
//                                           "langm_coeffs",
//                                           "matrix_volume",
//                                           "solver_method")))
//
//             .def("cfd_procedure_pnm_diff",
//                  &Aggregator::cfdProcedurePnmDiff)
//
//             .def("get_pressure_av",
//                  &Aggregator::getPressureAverage)
//
//             .def("get_matrix_mass_total",
//                  &Aggregator::getMatrixMassTotal)
//
//             .def("get_flow_pores_out",
//                  &Aggregator::getTotalFlowPoresOut)
//
//             .def("get_flow_pores_in",
//                  &Aggregator::getTotalFlowPoresIn)
//
//             .def("get_flow_diff",
//                  &Aggregator::getTotalFlowDiff)
//
//             .def("get_pressure_by_pore",
//                  &Aggregator::getPorePressure)
//
//             .def("get_inlet_pressure",
//                  &Aggregator::getInletPressure)
//
//             .def("get_time_steps_vec",
//                  &Aggregator::getTimeStepsVec);
// }




