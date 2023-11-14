cmake_minimum_required(VERSION 3.22)
project(pbrpc)

add_library(GeoModelTemplate SHARED GeoModelTemplate.cpp)
target_link_libraries(pbrpc GeoModelTemplate)
#add_executable(main main.cpp Student.cpp)

find_package(pybind11 REQUIRED)
pybind11_add_module(pbrpc bind_helloworld.cpp)
