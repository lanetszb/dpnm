#ifndef PYTHONCONVERSIONS_H
#define PYTHONCONVERSIONS_H


#include <vector>
#include <cstdlib>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>


namespace p = boost::python;
namespace np = boost::python::numpy;


template<typename T>
np::ndarray convertVectorCppToNpArray(const std::vector<T> &vectorCpp) {	
	p::list listP;
	for (auto element : vectorCpp)
		listP.append(p::object(T(element)));	
    return np::array(listP);
}

template<typename T>
struct VectorCppToNpArray {
    static PyObject *convert(const std::vector<T> &vectorCpp) {
        return p::incref(convertVectorCppToNpArray<T>(vectorCpp).ptr());
    }
};

template<typename T>
struct custom_vector_from_seq {
    custom_vector_from_seq() {
        p::converter::registry::push_back(&convertible, &construct,
                                          p::type_id<std::vector<T>>());
    }

    static void *convertible(PyObject *obj) {
        return PySequence_Check(obj) ? obj : nullptr;
    }

    static void construct(PyObject *obj,
                          p::converter::rvalue_from_python_stage1_data *data) {
        void *storage =
                ((p::converter::rvalue_from_python_storage<std::vector<T>> *)
                        (data))->storage.bytes;
        new(storage) std::vector<T>();
        std::vector<T> *v = (std::vector<T> *) (storage);
        int size = PySequence_Size(obj);
        if (size < 0)
            abort();
        v->reserve(size);
        for (int i = 0; i < size; i++) {
            v->push_back(p::extract<T>(PySequence_GetItem(obj, i)));
        }
        data->convertible = storage;
    }
};


void providePythonThings() {
	
	setenv("PYTHONPATH", ".", 1);
	
    Py_Initialize();
    np::initialize();
	
    p::to_python_converter<std::vector<bool>, VectorCppToNpArray<bool>>();
    custom_vector_from_seq<bool>();
	 
    p::to_python_converter<std::vector<int>, VectorCppToNpArray<int>>();
    custom_vector_from_seq<int>();
	
    p::to_python_converter<std::vector<double>, VectorCppToNpArray<double>>();
    custom_vector_from_seq<double>();
	
    p::to_python_converter<std::vector<std::string>, VectorCppToNpArray<std::string>>();
	custom_vector_from_seq<std::string>();
}


template<class T>
std::string __str__(T const &t) {
    std::stringstream stream;
    stream << t;
    return stream.str();
}


#endif // PYTHONCONVERSIONS_H
