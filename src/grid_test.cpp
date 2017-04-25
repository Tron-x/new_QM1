/**
 * \file grid_test.cpp
 * \brief Test grid
 */

#include "grid.hpp"

using namespace qm;
int main() {
  Grid grid;
  printf("setting up\n");
  grid.setup();

  printf("counting\n");
  printf("Total number of points:%d\n", grid.get_n());
}
