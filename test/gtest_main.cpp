// gtest_main.cpp
#include <stdio.h>
#include "gtest/gtest.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from gtest_main.cc\n");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
