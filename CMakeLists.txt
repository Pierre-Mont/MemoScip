cmake_minimum_required(VERSION 3.14)
project(MemoSCIP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Modes (exclusifs) :
#   défaut (les deux OFF)     → performant (-O3, NDEBUG), pas orienté perf
#   -DMEMOSCIP_PROFILE=ON     → -O2 -g -fno-omit-frame-pointer, NDEBUG (perf, etc.)
#   -DMEMOSCIP_DEBUG=ON       → -O0 -g, sanitizers, pas de NDEBUG


#cmake ..

#cmake .. -DMEMOSCIP_PROFILE=ON

#cmake .. -DMEMOSCIP_DEBUG=ON

#cmake .. -DSCIP_DIR=/path/to/scip-master

option(MEMOSCIP_PROFILE "Profiling-friendly build (-O2 -g, frame pointers)" OFF)
option(MEMOSCIP_DEBUG "Debug build with AddressSanitizer / UBSan" OFF)

set(SCIP_DIR "$ENV{HOME}/local/scip-master" CACHE PATH
    "Préfixe d'installation de SCIP pour find_package (ex. -DSCIP_DIR=/opt/scip)")
if(NOT SCIP_DIR)
    message(FATAL_ERROR
        "SCIP_DIR est vide. Indiquez l'emplacement d'installation de SCIP, par exemple : "
        "cmake -B <build> -DSCIP_DIR=/chemin/vers/scip")
endif()
list(INSERT CMAKE_PREFIX_PATH 0 "${SCIP_DIR}")
message(STATUS "SCIP_DIR : ${SCIP_DIR}")
message(STATUS "CMAKE_PREFIX_PATH (SCIP en tête) : ${CMAKE_PREFIX_PATH}")

find_package(SCIP REQUIRED CONFIG)
find_package(LibLZMA REQUIRED)

add_executable(memoscipl
    src/main.cpp
    src/memoscip_option.cpp
    src/memoscip_types.cpp
    src/memoscip_dual.cpp
    src/memoscip_events.cpp
    src/memoscip_reader.cpp
    src/sepa_mycuts.cpp
)

target_link_libraries(memoscipl PRIVATE libscip LibLZMA::LibLZMA)

target_include_directories(memoscipl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(MEMOSCIP_DEBUG)
    target_compile_options(memoscipl PRIVATE
        -g
        -O0
        -Wall
        -fsanitize=address
        -fsanitize=undefined
        -fno-omit-frame-pointer
    )
    target_link_options(memoscipl PRIVATE
        -fsanitize=address
        -fsanitize=undefined
    )
elseif(MEMOSCIP_PROFILE)
    target_compile_options(memoscipl PRIVATE
        -O2
        -g
        -fno-omit-frame-pointer
    )
    target_compile_definitions(memoscipl PRIVATE NDEBUG)
else()
    target_compile_options(memoscipl PRIVATE -O3)
    target_compile_definitions(memoscipl PRIVATE NDEBUG)
endif()
