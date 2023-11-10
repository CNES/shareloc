#include <pybind11/pybind11.h>
#include "hello_world.cpp"

namespace py = pybind11;


PYBIND11_MODULE(pbhelloworld, m) {
    py::class_<HW>(m, "HW")
        .def(py::init<>())
        .def("hellow_world", &HW::hellow_world);
    /*m.doc() = "Pybind hello world"; // optional module docstring
    m.def("main", &main, "Print hello world str");*/
}

//c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) bind_helloworld.cpp -o pbhelloworld$(python3-config --extension-suffix)