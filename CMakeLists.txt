cmake_minimum_required(VERSION 3.15)

project(PNflow)

set(CMAKE_CXX_STANDARD 17)

find_package(PythonLibs 3.7 REQUIRED)

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS python37 numpy37)

find_package(Eigen3 REQUIRED)


include_directories(
        numerical/equation
        numerical/math
        numerical/utilities)


python_add_module(cfd
        numerical/wrapper.cpp
        numerical/equation/Equation.cpp
        numerical/math/Props.cpp
        numerical/math/Local.cpp
        numerical/math/Convective.cpp)


target_include_directories(cfd PUBLIC
        ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})

target_link_libraries(cfd
        ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${EIGEN3_LIBRARIES})

set_target_properties(cfd PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../numerical)
