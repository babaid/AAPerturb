enable_testing()

add_executable(
        run_test
        run_test.cc
)
target_link_libraries(
        run_test
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(run_test)