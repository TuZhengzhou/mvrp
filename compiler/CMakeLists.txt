include_directories(./src)
include_directories(./test)
include_directories(../mvrp/src)
include_directories(../mvrp/includes)

find_package(OpenSSL REQUIRED)

add_executable(
    test_linear_combination 
    
    ./test/test_linear_combination.cpp 
    ./src/linear_combination.cpp
)

target_link_libraries(
    test_linear_combination

    ff
    gmpxx
    OpenSSL::Crypto
)

add_executable(
    test_shift_and_reverse

    ./test/test_shift_and_reverse.cpp
    ../mvrp/src/basic_types.cpp
    ../mvrp/src/math_operations.cpp
)

target_link_libraries(
    test_shift_and_reverse

    ff
    gmpxx
    OpenSSL::Crypto
)