cmake_minimum_required(VERSION 3.12)

include_directories(./)
include_directories(../epee/)
include_directories(../easylogging++/)

aux_source_directory(. COMMON_SRC)
message("COMMON_SRC: ${COMMON_SRC}")

add_library(b_common STATIC ${COMMON_SRC})
install(TARGETS b_common
    LIBRARY DESTINATION ./lib/)
# cd common
# gcc -c -fpic *.c *.cpp -I ./ -I ../epee/ -I ../easylogging++    ./common 目录下
# gcc -shared *.o -o ../lib/libcommon.so

#  common 
# epee  easylogging++