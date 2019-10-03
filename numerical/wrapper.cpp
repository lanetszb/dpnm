#include <boost/python.hpp>

#include <string>

#include <Props.h>
#include <Local.h>
#include <Convective.h>

#include <PythonConversions.h>

namespace p = boost::python;

BOOST_PYTHON_MODULE (cfd) {

    providePythonThings();

    p::class_<Props>("PropsCpp",
                     p::init<std::vector<double>>(
                             p::args("props_array")))
            .def("__str__", __str__<Props>)

            .add_property("length",
                          &Props::getLength,
                          &Props::setLength)
            .add_property("radius",
                          &Props::getRadius,
                          &Props::setRadius)
            .add_property("eff_radius",
                          &Props::getEffRadius,
                          &Props::setEffRadius)
            .add_property("grid_block_n",
                          &Props::getGridBlockN,
                          &Props::setGridBlockN)
            .add_property("concentration",
                          &Props::getConcentration,
                          &Props::setConcentration)
            .add_property("diffusivity",
                          &Props::getDiffusivity,
                          &Props::setDiffusivity)
            .add_property("iterative_accuracy",
                          &Props::getIterativeAccuracy,
                          &Props::setIterativeAccuracy)
            .add_property("props_vector",
                          &Props::getPropsVector)

            .def("print_props_vector",
                 &Props::printPropsVector);


    p::class_<Local>("LocalCpp",
                     p::init<std::vector<double>>(
                             p::args("props_array")))


            .def("calc_concListIni",
                 &Local::calc_concListIni,
                 p::args("grid_block_n",
                         "concentration"))


            .def("calculateAlpha",
                 &Local::calculateAlpha,
                 p::args("dt",
                         "radius",
                         "effRadius",
                         "gridBlockN"))

            .def("get_alpha",
                 &Local::getAlpha);


    p::class_<Convective>("ConvectiveCpp",
                     p::init<std::vector<double>>(
                             p::args("props_array")))


            .def("calc_diffusivityList",
                 &Convective::calc_diffusivityList,
                 p::args("grid_block_n",
                         "concentration"))

            .def("calculateBeta",
                 &Convective::calculateBeta,
                 p::args("radius",
                         "effRadius",
                         "length",
                         "diffusivity",
                         "gridBlockN"))

            .def("get_beta",
                 &Convective::getBeta);
}



