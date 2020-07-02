cmake_minimum_required(VERSION 3.15)

project(PNflow)

set(CMAKE_CXX_STANDARD 17)

set(PYTHON_LIBRARY "/usr/local/Cellar/python@3.8/3.8.3/Frameworks/Python.framework/Versions/3.8/lib/libpython3.8.dylib")
set(PYTHON_INCLUDE_DIR "/usr/local/Cellar/python@3.8/3.8.3/Frameworks/Python.framework/Versions/3.8/include/python3.8")
find_package(PythonLibs)

set(Boost_NO_BOOST_CMAKE TRUE)

find_package(Boost REQUIRED COMPONENTS python38 numpy38)

find_package(Eigen3 REQUIRED)


include_directories(
        numerical/coupling
        numerical/equation_pnm
        numerical/math_pnm
        numerical/equation_diff
        numerical/math_diff
        numerical/utilities)


python_add_module(cfd
        numerical/wrapper.cpp
        numerical/coupling/CouplingIniConds.cpp
        numerical/coupling/CouplingParamsOut.cpp
        numerical/coupling/CouplingMatrix.cpp
        numerical/coupling/DiffusionPNM.cpp
        numerical/coupling/DiffusionPartMath.cpp
        numerical/coupling/DiffusionPartFlow.cpp
        numerical/equation_pnm/EquationPNM.cpp
        numerical/math_pnm/PropsPNM.cpp
        numerical/math_pnm/NetworkData.cpp
        numerical/equation_diff/EquationDiffusion.cpp
        numerical/math_diff/PropsDiffusion.cpp
        numerical/math_diff/LocalDiffusion.cpp
        numerical/math_diff/ConvectiveDiffusion.cpp)

target_compile_definitions(cfd PRIVATE BOOST_PYTHON_MAX_ARITY=30)

target_include_directories(cfd PUBLIC
        ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})

target_link_libraries(cfd
        ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${EIGEN3_LIBRARIES})

set_target_properties(cfd PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../numerical)
