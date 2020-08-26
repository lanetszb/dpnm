cmake_minimum_required(VERSION 3.15)
project(dpnm)

set(CMAKE_CXX_STANDARD 17)

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
        numerical/coupling
        numerical/equation_pnm
        numerical/math_pnm
        numerical/equation_diff
        numerical/math_diff)

pybind11_add_module(${PROJECT_NAME}
        numerical/wrapper.cpp
        numerical/coupling/Aggregator.cpp
        numerical/coupling/IniConds.cpp
        numerical/coupling/ParamsOut.cpp
        numerical/coupling/MatrixSolver.cpp
        numerical/coupling/DiffusionMath.cpp
        numerical/coupling/DiffusionFlow.cpp
        numerical/equation_pnm/EquationPNM.cpp
        numerical/math_pnm/NetworkData.cpp
        numerical/equation_diff/EquationDiffusion.cpp
        numerical/math_diff/PropsDiffusion.cpp
        numerical/math_diff/LocalDiffusion.cpp
        numerical/math_diff/ConvectiveDiffusion.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} ${EIGEN3_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../numerical)
