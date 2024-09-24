project(NPP_Roe CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_EXE_LINKER_FLAGS "-static")

cmake_minimum_required(VERSION 3.0)

include_directories(./include)

aux_source_directory(./src DIRSRCS)

add_executable(NPP_Roe main.cpp ${DIRSRCS})

