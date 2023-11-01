cmake_minimum_required(VERSION 2.8)
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(mvrp)
SET(PROJECT_VERSION 1.2.3.4)

set(CMAKE_BUILD_TYPE "Debug" )
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

add_subdirectory(our_cred)
# add_subdirectory(bulletproof)
# add_subdirectory(libvectcom)
add_subdirectory(libsnark)
# add_subdirectory(depends)