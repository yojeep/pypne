cmake_minimum_required(VERSION 4.2)
# set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
# set(CMAKE_CXX_COMPILER "g++") # only support g++ or MinGW-w64

project(pypne_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # 必须使用 C++20，否则报错

set(pybind11_DIR ${CMAKE_SOURCE_DIR}/thirdparty/system/pybind11/share/cmake/pybind11)
find_package(pybind11 REQUIRED)
pybind11_add_module(pypne_cpp
${CMAKE_SOURCE_DIR}/src/pnextract/medialSurf.cpp
${CMAKE_SOURCE_DIR}/src/pnextract/blockNet.cpp
${CMAKE_SOURCE_DIR}/src/pnextract/pypne.cpp
${CMAKE_SOURCE_DIR}/src/pnextract/main.cpp
)

target_include_directories(pypne_cpp
PRIVATE
${CMAKE_SOURCE_DIR}/src/include
)

target_include_directories(pypne_cpp
SYSTEM PRIVATE
${CMAKE_SOURCE_DIR}/thirdparty/system
${CMAKE_SOURCE_DIR}/thirdparty/system/pybind11/include
${CMAKE_SOURCE_DIR}/thirdparty/system/eigen-5.0.1
${CMAKE_SOURCE_DIR}/thirdparty/system/eigen-5.0.1/unsupported
)


target_compile_definitions(pypne_cpp PRIVATE BOOST_DISABLE_ASSERTS)
if(MSVC)
    target_compile_options(pypne_cpp PRIVATE /O2 /W1 /GL /LTCG /MT)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pypne_cpp PRIVATE -O3 -Wall)
    target_link_options(pypne_cpp PRIVATE)
endif()

install(TARGETS pypne_cpp LIBRARY DESTINATION pypne/libcpp)