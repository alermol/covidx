cmake_minimum_required(VERSION 3.20)
project(covidx_consumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(covidx CONFIG REQUIRED)

add_executable(covidx_consumer main.cpp)
target_link_libraries(covidx_consumer PRIVATE covidx::covidx)
