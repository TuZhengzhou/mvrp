cmake_minimum_required(VERSION 3.12)

aux_source_directory(. EASYLOGGING_SRC)

message("EASYLOGGING_SRC: ${EASYLOGGING_SRC}")

add_library(b_easylogging STATIC ${EASYLOGGING_SRC})
install(TARGETS b_easylogging
    LIBRARY DESTINATION ./lib/)