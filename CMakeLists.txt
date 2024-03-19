cmake_minimum_required(VERSION 3.22)


project(AAPerturb VERSION 1.2.0)

option(BUILD_TESTS "Enable building the tests" OFF)
set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -fsanitize=address -O2 -Wsign-conversion -Wconversion -pedantic")
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

add_subdirectory(Submodules/spdlog)


if(NOT TARGET spdlog)
    # Stand-alone build
    find_package(spdlog REQUIRED)
endif()

file(GLOB_RECURSE SOURCES "source/*.cpp" "include/*.h")

list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/main.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/find_interfaces.cpp")

add_executable(aaperturb source/main.cpp ${SOURCES})
add_executable(find_interfaces source/find_interfaces.cpp ${SOURCES})



target_link_libraries(aaperturb argparse spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)
target_link_libraries(find_interfaces argparse spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)


install(TARGETS aaperturb DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS find_interfaces DESTINATION ${CMAKE_INSTALL_PREFIX})

if (BUILD_TESTS)
    add_subdirectory(unittest)
endif ()