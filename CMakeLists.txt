cmake_minimum_required(VERSION 3.0)

project(RabbitTrim)

set(CMAKE_INSTALL_PREFIX ..)
set(EXECUTABLE_OUTPUT_PATH .)

set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(COMMON_FLAGS "-O3")

# cmake build type
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "debug")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "")

if("${CMAKE_BUILD_TYPE}" STREQUAL "release")
    set(COMMON_FLAGS "${COMMON_FLAGS} ")
    message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}, compile ${CMAKE_PROJECT_NAME} with ${COMMON_FLAGS} flag.")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "debug")
    set(COMMON_FLAGS "${COMMON_FLAGS} -g -Wall")
    message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}, compile ${CMAKE_PROJECT_NAME} with ${COMMON_FLAGS} flag.")
endif("${CMAKE_BUILD_TYPE}" STREQUAL "release")

# compile flags
set(CMAKE_C_FLAGS "${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS "${COMMON_FLAGS}")
message(STATUS "CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")

aux_source_directory(${CMAKE_SOURCE_DIR}/src DIR_SRC)

# include files
set(INCLUlDE_PATH ${CMAKE_SOURCE_DIR}/include)
set(IO_INCLUlDE_PATH ${CMAKE_SOURCE_DIR}/src/io)
# include_directories(${INCLUlDE_PATH} ${INCLUlDE_PATH}/core ${INCLUlDE_PATH}/trimmer ${INCLUlDE_PATH}/trimlog)
include_directories(${INCLUlDE_PATH})
include_directories(${IO_INCLUlDE_PATH}) 

# source files
set(SOURCE_PATH ${CMAKE_SOURCE_DIR}/src)
set(IO_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/io)
aux_source_directory(${SOURCE_PATH}/core DIR_SRC_CORE)
# add_library(core STATIC ${DIR_SRC_CORE})
aux_source_directory(${SOURCE_PATH}/trimmer DIR_SRC_TRIM)
# add_library(trim STATIC ${DIR_SRC_TRIM})
aux_source_directory(${SOURCE_PATH}/trimlog DIR_SRC_LOG)
# add_library(log STATIC ${DIR_SRC_LOG})
add_library(rtrim STATIC ${DIR_SRC_CORE} ${DIR_SRC_TRIM} ${DIR_SRC_LOG})
aux_source_directory(${IO_SOURCE_PATH} DIR_SRC_IO)
add_library(io STATIC ${DIR_SRC_IO})

# executable
aux_source_directory(${SOURCE_PATH} DIR_SRC)
add_executable(${PROJECT_NAME} ${DIR_SRC})

# libraries
# target_link_libraries(${PROJECT_NAME} PUBLIC core trim log io)
target_link_libraries(${PROJECT_NAME} PUBLIC rtrim io)
