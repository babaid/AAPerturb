cmake_minimum_required(VERSION 3.26)

project(unittests)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

include(GoogleTest)

set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest")

add_executable(
        unittest_Geometry
        source/unittests_Geometry.cpp
        ${SOURCES}
)
target_link_libraries(
        unittest_Geometry
        GTest::gtest_main
        spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>
)


gtest_discover_tests(unittest_Geometry)