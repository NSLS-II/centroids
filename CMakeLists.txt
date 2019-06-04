cmake_minimum_required(VERSION 3.14)
project(centroids)

set(CMAKE_CXX_STANDARD 14)

find_package(pybind11 REQUIRED)
pybind11_add_module(pybind11 src/binding.cpp)

add_library(centroids src/photons.cpp src/photons.h)
