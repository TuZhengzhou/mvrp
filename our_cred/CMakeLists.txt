set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DCRED_DEBUG)
endif()

find_package(OpenSSL REQUIRED)

include_directories(/usr/include)
include_directories(/usr/include/x86_64-linux-gnu)
include_directories(../depends/libff)
include_directories(includes)
include_directories(src)

link_directories(/usr/local/lib)
link_directories(/usr/lib)

add_executable(
  cred 

  ./src/main.cpp
  ./src/structs.cpp
  ./src/prover.cpp
  ./src/verifier.cpp
  ./src/ipa.cpp
  ./src/kzg.cpp
  ./src/circuit_generator.cpp
)

target_link_libraries(
  cred 

  ff
  gmpxx
  OpenSSL::Crypto
)

add_executable(
  test_naive 

  ./src/test_naive.cpp
  ./src/structs.cpp
  ./src/prover.cpp
  ./src/verifier.cpp
  ./src/ipa.cpp
  ./src/kzg.cpp
  ./src/circuit_generator.cpp
)

target_link_libraries(
  test_naive 

  ff
  gmpxx
  OpenSSL::Crypto
)

add_subdirectory(src)
