cmake_minimum_required(VERSION 3.12)
project(CMakeTemplate VERSION 1.0.0 LANGUAGES C CXX DESCRIPTION "monero-bulletproof project")
set (CMAKE_CXX_COMPILER "/usr/bin/g++")
# Specified the language standard
set(CMAKE_CXX_STANDARD 14)

##--------------------- Compile Options ------------------------------------##
# Configure compile options
add_compile_options(-Wextra)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -D_FORCE_INLINES -std=c++14")

# Set build type
# set(CMAKE_BUILD_TYPE Debug)  # Use `cmake -DCMAKE_BUILD_TYPE=Debug` more better
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Compile options for Debug variant
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Compile options for Release variant
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

message(STATUS "Compile options for c: ${CMAKE_C_FLAGS}")
message(STATUS "Compile options for c++: ${CMAKE_CXX_FLAGS}")

##--------------------- Global Macros --------------------------------------##
add_definitions(-DDEBUG -DREAL_COOL_ENGINEER)


# link_directories(/usr/local/lib/)

add_subdirectory(./easylogging++)
add_subdirectory(./epee)
add_subdirectory(./common)
add_subdirectory(./crypto)
add_subdirectory(./ringct)
# add_subdirectory(./serialization)

##--------------------- Include directories --------------------------------##
include_directories(.)
include_directories(./epee)
# set(LIB_DIR "./lib/ /usr/local/lib/")
# link_libraries(${LIB_DIR})
link_directories(./lib/)
link_directories(/usr/local/lib/)

# Build demo executable
add_executable(demo main.cpp)
target_link_libraries(
    demo 
    
    b_ringct
    b_crypto
    b_common 
    b_epee
    b_easylogging
     
    # b_serialization
    boost_thread
    boost_chrono

    libpthread.so # 如果不增加libpthread.so也是不能通过编译，是因为linux系统有时候默认找不到pthread
    )

# https://blog.csdn.net/hunterhe/article/details/120591117
# cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
# project(lidar_demo)
# find_package(OpenCV REQUIRED)
# set(lib_DIR "～/tools/Neuvition_tools/neusdk-x86_64-2.5.5/sdk/lib")
# set(include_DIR "～/tools/Neuvition_tools/neusdk-x86_64-2.5.5/sdk/include")
# include_directories(${OpenCV_INCLUDE_DIRS} "～/tools/Neuvition_tools/neusdk-x86_64-2.5.5/sdk/include")
# link_directories(${OpenCV_LIBRARY_DIRS} "～/tools/Neuvition_tools/neusdk-x86_64-2.5.5/sdk/lib")
# add_executable (lidar_demo main.cpp)
# target_link_libraries (lidar_demo pthread curl uuid libpthread.so pthread libneusdk.a ${OpenCV_LIBS})
