//
// Created by babaid on 11.09.23.
//
#include <gtest/gtest.h>
#include<map>
#include<vector>
#include "molecules.h"
#include "../include/pdbparser.h"
// Demonstrate some basic assertions.
TEST(ParsePDBTEST, BasicAssertions) {
    // Expect two strings not to be equal.
    std::map<char, std::vector<Residue*>> pdb = parsePDB()
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}



