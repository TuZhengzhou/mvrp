set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DCRED_DEBUG)
endif()

find_package(OpenSSL REQUIRED)

file(GLOB_RECURSE SRC_FILES
  "./src/*.cpp"
)

set(LINK_LIBRARIES

  mvrp
  ff
  gmpxx
  OpenSSL::Crypto
)

include_directories(/usr/include)
include_directories(/usr/include/x86_64-linux-gnu)
include_directories(../depends/libff)
include_directories(../compiler/src)
include_directories(includes)
include_directories(src)
include_directories(tests)

link_directories(/usr/local/lib)
link_directories(/usr/lib)

add_library(
  mvrp
  STATIC

  ${SRC_FILES}
)

target_link_libraries(
  mvrp

  ff
  gmpxx
  OpenSSL::Crypto
)

install(
  DIRECTORY "" DESTINATION "include/mvrp"
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.tcc"
  PATTERN "tests" EXCLUDE
  PATTERN "examples" EXCLUDE
)

install(
  TARGETS mvrp DESTINATION lib  # /usr/local/lib
)

add_executable(
  test_range 

  ./tests/test_range.cpp
)

target_link_libraries(
  test_range 

  # mvrp
  # ff
  # gmpxx
  # OpenSSL::Crypto
  ${LINK_LIBRARIES}
)

add_executable(
  test_circuit

  ./tests/test_circuit.cpp
)
target_link_libraries(
  test_circuit

  mvrp
  ff
  gmpxx
  OpenSSL::Crypto
)

add_executable(
  test_circuit_gen 

  ./tests/test_circuit_gen.cpp
)

target_link_libraries(
  test_circuit_gen 

  mvrp
  ff
  gmpxx
  OpenSSL::Crypto
)

add_executable(
  microbenchmark 

  ./tests/microbenchmark.cpp
)

target_link_libraries(
  microbenchmark 

  mvrp
  ff
  gmpxx
  OpenSSL::Crypto
)

add_executable(
  test_poly_merge

  ./tests/test_poly_merge.cpp
)

target_link_libraries(
  test_poly_merge 

  mvrp
  ff
  gmpxx
  OpenSSL::Crypto
)

add_executable(
  test_matrix_mult

  ./tests/test_matrix_mult.cpp
)

target_link_libraries(
  test_matrix_mult 

  mvrp
  ff
  gmpxx
  OpenSSL::Crypto
)

add_executable(
  test_linear_combination 
    
  ./tests/test_linear_combination.cpp 
)

target_link_libraries(
    test_linear_combination

    mvrp
    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(
  test_shift_and_reverse

  ./tests/test_shift_and_reverse.cpp
)

target_link_libraries(
  test_shift_and_reverse

  mvrp
  ff
  gmpxx
  OpenSSL::Crypto
)

add_test(
  NAME test_circuit
  COMMAND test_circuit
)

add_test(
  NAME test_circuit_gen
  COMMAND test_circuit_gen
)

add_test(
  NAME test_range
  COMMAND test_range
)

add_test(
  NAME test_poly_merge
  COMMAND test_poly_merge
)

# add_test(
#   NAME test_matrix_mult
#   COMMAND test_matrix_mult
# )

add_test(
  NAME test_linear_combination
  COMMAND test_linear_combination
)

add_test(
  NAME test_shift_and_reverse
  COMMAND test_shift_and_reverse
)
