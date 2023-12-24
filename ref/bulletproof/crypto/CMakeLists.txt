cmake_minimum_required(VERSION 3.12)

include_directories(./)
include_directories(../epee/)
include_directories(../)

aux_source_directory(. CRYPTO_SRC)
aux_source_directory(./crypto_ops_builder CRYPTO_OPS_BUILDER_SRC)
message("CRYPTO_SRC: ${CRYPTO_SRC}")
message("CRYPTO_OPS_BUILDER_SRC: ${CRYPTO_OPS_BUILDER_SRC}")

add_library(b_crypto STATIC ${CRYPTO_SRC} ${CRYPTO_OPS_BUILDER_SRC})
install(TARGETS b_crypto
    LIBRARY DESTINATION ./lib/)
# cd crypto
# gcc -c -fpic *.c *.cpp ./crypto_ops_builder/*.c -I ./crypto_ops_builder/ -I ./ -I ../epee/ -I ../  ./crypto 目录下
# gcc -shared *.o -o ../lib/libcrypto.so