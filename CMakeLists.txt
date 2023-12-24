cmake_minimum_required(VERSION 2.8)
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(mvrp)
SET(PROJECT_VERSION 1.2.3.4)


set(
  CURVE
  "BN128"
  CACHE
  STRING
  "Default curve: one of ALT_BN128, BN128, EDWARDS, MNT4, MNT6"
)

add_definitions(
  -DCURVE_BN128
)
if(${CURVE} STREQUAL "BN128")
  add_definitions(
    -DBN_SUPPORT_SNARK=1
  )
endif()

# set(CMAKE_BUILD_TYPE "Debug" )
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -Wendif-labels -g2 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -Wendif-labels")

# ADD_DEFINITIONS(-DPROFILE_OP_COUNTS)

enable_testing()

add_subdirectory(depends)
add_subdirectory(mvrp)
# add_subdirectory(compiler)