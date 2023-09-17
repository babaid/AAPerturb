cmake_minimum_required(VERSION 3.22)
project(AA_Sidechain_Perturbation_CPP)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -fsanitize=address -O2 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
set(CMAKE_CXX_STANDARD_REQUIRED  ON)
set(CMAKE_CXX_STANDARD 23)

include_directories(include)

include(FetchContent)
FetchContent_Declare(
        argparse
        GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

#FetchContent_Declare(
#        googletest
#        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
#)
# For Windows: Prevent overriding the parent project's compiler/linker settings
#set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#FetchContent_MakeAvailable(googletest)


add_executable(aaperturb
        include/geometry.h
        include/pdbparser.h
        include/constants.h
        src/geometry.cpp
        src/pdbparser.cpp
        src/main.cpp
        include/io.h
        src/io.cpp
        include/montecarlo.h
        src/montecarlo.cpp
        include/molecules.h
        src/molecules.cpp
        include/fancy.h
        src/fancy.cpp
        include/threadpool.h
)

add_executable(iffind
        include/geometry.h
        include/pdbparser.h
        include/constants.h
        src/geometry.cpp
        src/pdbparser.cpp
        include/io.h
        src/io.cpp
        include/molecules.h
        src/molecules.cpp
        src/iffind.cpp
)
add_executable(pdbcleaner
        include/pdbparser.h
        src/pdbparser.cpp
        include/io.h
        src/io.cpp
        include/molecules.h
        src/molecules.cpp
        include/fancy.h
        src/fancy.cpp
        src/pdbcleaner.cpp
)
target_link_libraries(iffind argparse)
target_link_libraries(aaperturb argparse)
target_link_libraries(pdbcleaner argparse)
