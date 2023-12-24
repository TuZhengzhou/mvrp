cmake_minimum_required(VERSION 3.12)

include_directories(./)
include_directories(../)
include_directories(../epee/)
include_directories(../easylogging++/)

aux_source_directory(. RINGCT)
message("RINGCT: ${RINGCT}")

add_library(b_ringct STATIC ${RINGCT})
install(TARGETS b_ringct
    LIBRARY DESTINATION ./lib/)
# cd ringct
# gcc -c gcc -c -fpic ./src/c/math/add.c ./src/c/math/minus.c -I ./src/c/ 
# gcc -shared *.o -o ./lib/libmath.so