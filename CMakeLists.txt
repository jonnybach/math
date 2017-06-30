cmake_minimum_required(VERSION 2.8.12)
project(sie_math)

#Add versioning
set (PRODUCT_VERSION_MAJOR 0)
set (PRODUCT_VERSION_MINOR 1)
set (PRODUCT_VERSION_BUILD 0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

add_subdirectory(./src)