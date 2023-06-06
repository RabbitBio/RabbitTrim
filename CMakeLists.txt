cmake_minimum_required(VERSION 3.20)

project(RabbitTrim)

set(CMAKE_INSTALL_PREFIX ..)
set(EXECUTABLE_OUTPUT_PATH .)

set(CMAKE_C_COMPILER "/home/user_home/public/application/gcc-10.1.0/bin/gcc")
set(CMAKE_CXX_COMPILER "/home/user_home/public/application/gcc-10.1.0/bin/g++")
message(${CMAKE_CXX_COMPILER})

# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(COMMON_FLAGS "-O3 -march=native -fopt-info-vec -g -msse2 -mavx -mavx2")
# set(COMMON_FLAGS "-g -DTRIM_DEBUG -O3 -march=native -fconstexpr-ops-limit=99000100 ")
set(COMMON_FLAGS "-DTRIM_DEBUG -O3 -march=native -fconstexpr-ops-limit=99000100 ")

option(USE_OMP "Specify whether to use openMP" OFF)
option(TRIM_USE_VEC "Specify whether to use vectorization" ON)

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
set(PRAGZIP_INCLUlDE_PATH ${CMAKE_SOURCE_DIR}/src/pragzip)
include_directories(${INCLUlDE_PATH} ${IO_INCLUlDE_PATH} ${PIGZ_INCLUlDE_PATH} ${PRAGZIP_INCLUlDE_PATH})

# source files
set(SOURCE_PATH ${CMAKE_SOURCE_DIR}/src)
set(IO_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/io)
set(PIGZ_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/pigz)


aux_source_directory(${IO_SOURCE_PATH} DIR_SRC_IO)
add_library(io STATIC ${DIR_SRC_IO})
aux_source_directory(${PIGZ_SOURCE_PATH} DIR_SRC_PIGZ)
add_library(pigz STATIC ${DIR_SRC_PIGZ})

add_subdirectory(${CMAKE_SOURCE_DIR}/src/pragzip)

aux_source_directory(${SOURCE_PATH} DIR_SRC)
foreach (tmp IN ITEMS handler trimmer trimlog utils)
  aux_source_directory(${SOURCE_PATH}/${tmp} DIR_TMP)
  list(APPEND DIR_SRC ${DIR_TMP})
endforeach()

# executable
add_executable(${PROJECT_NAME} ${DIR_SRC})

# libraries
if(DEFINED IGZIP_PATH) 
  target_link_libraries(${PROJECT_NAME} PUBLIC io isal pigz libpragzip z pthread stdc++fs)
else()
  target_link_libraries(${PROJECT_NAME} PUBLIC io pigz libpragzip z pthread)
endif()


