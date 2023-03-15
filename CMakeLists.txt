cmake_minimum_required(VERSION 3.0)

project(RabbitTrim)

set(CMAKE_INSTALL_PREFIX ..)
set(EXECUTABLE_OUTPUT_PATH .)

set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(COMMON_FLAGS "-O3 -march=native -fopt-info-vec -g -msse2 -mavx -mavx2")

option(USE_OMP "Specify whether to use openMP" OFF)
option(TRIM_USE_VEC "Specify whether to use vectorization" OFF)

# cmake build type
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "release")
    # set(CMAKE_BUILD_TYPE "debug")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "")

if("${CMAKE_BUILD_TYPE}" STREQUAL "release")
    set(COMMON_FLAGS "${COMMON_FLAGS} -w")
    message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}, compile ${CMAKE_PROJECT_NAME} with ${COMMON_FLAGS} flag.")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "debug")
    set(COMMON_FLAGS "${COMMON_FLAGS} -g -Wall -fsanitize=address")
    message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}, compile ${CMAKE_PROJECT_NAME} with ${COMMON_FLAGS} flag.")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "release")

if(USE_OMP)
  set(COMMON_FLAGS "${COMMON_FLAGS} -fopenmp")
endif()
if(TRIM_USE_VEC)
  set(COMMON_FLAGS "${COMMON_FLAGS} -DTRIM_USE_VEC")
endif()



if(DEFINED IGZIP_PATH)
  message(STATUS "Using igzip : ${IGZIP_PATH}")
  include_directories(${IGZIP_PATH}/include/isa-l)
  link_directories(${IGZIP_PATH}/lib64)
  link_directories(${IGZIP_PATH}/lib32)
  link_directories(${IGZIP_PATH}/lib)
  set(COMMON_FLAGS "${COMMON_FLAGS} -DUSE_IGZIP")
endif()

# compile flags
set(CMAKE_C_FLAGS "${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS "${COMMON_FLAGS}")
message(STATUS "CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")

# include files
set(INCLUlDE_PATH ${CMAKE_SOURCE_DIR}/include)
set(IO_INCLUlDE_PATH ${CMAKE_SOURCE_DIR}/src/io)
set(PIGZ_INCLUlDE_PATH ${CMAKE_SOURCE_DIR}/src/pigz)
include_directories(${INCLUlDE_PATH} ${IO_INCLUlDE_PATH} ${PIGZ_INCLUlDE_PATH}) 

# source files
set(SOURCE_PATH ${CMAKE_SOURCE_DIR}/src)
set(IO_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/io)
set(PIGZ_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/pigz)
aux_source_directory(${SOURCE_PATH}/handler DIR_SRC_HANDLER)
aux_source_directory(${SOURCE_PATH}/trimmer DIR_SRC_TRIM)
aux_source_directory(${SOURCE_PATH}/trimlog DIR_SRC_LOG)
aux_source_directory(${SOURCE_PATH}/utils DIR_SRC_UTIL)
add_library(rtrim STATIC ${DIR_SRC_HANDLER} ${DIR_SRC_TRIM} ${DIR_SRC_LOG} ${DIR_SRC_UTIL})
aux_source_directory(${IO_SOURCE_PATH} DIR_SRC_IO)
add_library(io STATIC ${DIR_SRC_IO})
aux_source_directory(${PIGZ_SOURCE_PATH} DIR_SRC_PIGZ)
add_library(pigz STATIC ${DIR_SRC_PIGZ})

# executable
aux_source_directory(${SOURCE_PATH} DIR_SRC)
add_executable(${PROJECT_NAME} ${DIR_SRC})

# libraries
if(DEFINED IGZIP_PATH) 
  target_link_libraries(${PROJECT_NAME} PUBLIC rtrim io isal pigz z pthread)
else()
  target_link_libraries(${PROJECT_NAME} PUBLIC rtrim io pigz z pthread)
endif()
  

