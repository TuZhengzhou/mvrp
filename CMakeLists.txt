cmake_minimum_required(VERSION 2.8)
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(mvrp)
SET(PROJECT_VERSION 1.2.3.4)

add_definitions(
  -DCURVE_BN128
)
if(${CURVE} STREQUAL "BN128")
  add_definitions(
    -DBN_SUPPORT_SNARK=1
  )
endif()

set(CMAKE_BUILD_TYPE "Debug" )
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -Wendif-labels -g2 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -Wendif-labels")

add_subdirectory(depends)
add_subdirectory(mvrp)