cmake_minimum_required(VERSION 4.2)
# set(CMAKE_C_COMPILER clang)
# set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_COMPILER "g++") # only support g++ or MinGW-w64

project(pypne_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # 必须使用 C++20，否则报错

# add_subdirectory(src/pnm/pnextract)
set(pybind11_DIR ${CMAKE_SOURCE_DIR}/thirdparty/pybind11/share/cmake/pybind11)
find_package(pybind11 REQUIRED)
pybind11_add_module(pypne_cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/blockNet.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/pypne.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/medialSurf.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/writers_vtk.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/writers_vxl.cpp
${CMAKE_SOURCE_DIR}/src/libvoxel/voxelImage.cpp)

target_include_directories(pypne_cpp
PRIVATE
${CMAKE_SOURCE_DIR}/src/include
${CMAKE_SOURCE_DIR}/src/libvoxel
)

target_include_directories(pypne_cpp
SYSTEM PRIVATE
${CMAKE_SOURCE_DIR}/thirdparty/boost
${CMAKE_SOURCE_DIR}/thirdparty/thread_pool
${CMAKE_SOURCE_DIR}/thirdparty/pybind11/include

)

target_compile_definitions(pypne_cpp PRIVATE BOOST_DISABLE_ASSERTS)
if(MSVC)
    target_compile_options(pypne_cpp PRIVATE /O2 /W1 /GL /LTCG /MT)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pypne_cpp PRIVATE -O3 -Wall)
endif()

install(TARGETS pypne_cpp LIBRARY DESTINATION pypne/libcpp)



