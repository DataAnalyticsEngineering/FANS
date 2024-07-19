pybind11_add_module(PyFANS micro.hpp micro.cpp)
target_link_libraries(PyFANS PRIVATE FANS::FANS)