cmake_minimum_required(VERSION 3.12)

include_directories(./)
include_directories(../easylogging++/)
include_directories(./storages/)

aux_source_directory(. EPEE_SRC)
message("EPEE_SRC: ${EPEE_SRC}")

add_library(b_epee STATIC ${EPEE_SRC})
install(TARGETS b_epee
    LIBRARY DESTINATION ./lib/)