cmake_minimum_required(VERSION 3.18)
project(pedHap)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE DEBUG)

set(CMAKE_EXPORT_COMPILECOMMANDS ON)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb -fkeep-inline-functions")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")
aux_source_directory(./src PED_HAP_SRC)
#aux_source_directory(./vcflib VCF_LIB)
add_executable(pedHap main.cpp ${PED_HAP_SRC} src/trioPhaser.cpp src/trioPhaser.h src/phaseInfo.cpp src/phaseInfo.h)

target_link_libraries(pedHap -lbz2 -lz -lhts)
install(TARGETS pedHap DESTINATION bin)
