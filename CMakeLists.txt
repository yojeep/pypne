cmake_minimum_required(VERSION 3.18)
# set(CMAKE_C_COMPILER clang)
# set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_COMPILER "g++") # only support g++ or MinGW-w64
project(pypne_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # 必须使用 C++17，否则报错
# set(CMAKE_BUILD_TYPE Release)

# set(CMAKE_BUILD_RPATH $ORIGIN) # 编译时从可执行文件所在目录查找
# set(CMAKE_INSTALL_RPATH $ORIGIN) # 运行时从可执行文件所在目录查找
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)


# add_subdirectory(src/pnm/pnextract)
find_package(pybind11 REQUIRED)
pybind11_add_module(pypne_cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/blockNet.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/pypne.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/medialSurf.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/writers_vtk.cpp
${CMAKE_SOURCE_DIR}/src/pnm/pnextract/writers_vxl.cpp
${CMAKE_SOURCE_DIR}/src/libvoxel/voxelImage.cpp)
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

target_include_directories(pypne_cpp
PRIVATE
${CMAKE_SOURCE_DIR}/src/include
${CMAKE_SOURCE_DIR}/src/libvoxel
${CMAKE_SOURCE_DIR}/thirdparty/pybind11-3.0.1/include
${CMAKE_SOURCE_DIR}/thirdparty/boost)

if(MSVC)
    target_compile_options(pypne_cpp PRIVATE /O2 /W1 /GL /LTGC /MT)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pypne_cpp PRIVATE -Ofast -Wall)
endif()


###### OpenMP Parallelization ######
option(USE_OMP "Enable OpenMP parallelization" ON)
if(USE_OMP)
    find_package(OpenMP REQUIRED)
    if(MSVC)
        # MSVC 编译器不支持 OpenMP4.0
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        # GCC/Clang 使用 -fopenmp
        target_compile_options(pypne_cpp PRIVATE ${OpenMP_CXX_FLAGS} -DOpenMP)
        target_link_libraries(pypne_cpp PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()
install(TARGETS pypne_cpp LIBRARY DESTINATION pypne/libcpp)
