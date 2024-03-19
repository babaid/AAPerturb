#include "geometry.h"
#include <gtest/gtest.h>

#include<cmath>

TEST(unittest_Geometry, BasicAssertions) {
// Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
// Expect equality.
    EXPECT_EQ(7 * 6, 42);
}

TEST(unittest_Geometry, VectorRotation2D)
{
    Vector3 axis{0., 0., 1.};
    Vector3 vec{1., 0., 0.};
    Vector3 pivot{0., 0., 0.};
    rotateCoordinatesAroundAxis(vec, pivot, axis, 90.0);
    Vector3 expected{0., 1., 0.};
    auto stddev = sum(vec - expected);
    EXPECT_NEAR(stddev, 0.0, 0.001);

    vec = {1., 0., 2.};
    expected = {0., 1., 2.};
    pivot = {0., 0., 1.};
    rotateCoordinatesAroundAxis(vec, pivot, axis, 90.);
    stddev = sum(vec - expected);
    EXPECT_NEAR(stddev, 0.0, 0.001);
}


TEST(unittest_Geometry, VectorRotation3D)
{
    Vector3 axis{0., 0., 1.};
    Vector3 vec{1., 0., 2.};
    Vector3 expected{0., 1., 2.};
    Vector3 pivot{0., 0., 1.};
    rotateCoordinatesAroundAxis(vec, pivot, axis, 90.);
    auto stddev = sum(vec - expected);
    EXPECT_NEAR(stddev, 0.0, 0.001);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

