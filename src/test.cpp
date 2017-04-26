#include <iostream>
#include <getopt.h>
#include <dirent.h>
#include <cstring>
#include <time.h>
#include <cmath>

#include "eft_calculator.hpp"
#include "test.h"

using namespace std;


int main(int argc, char* argv[]) {
  int order;
  string database_file = "grid.dat";
  string input_dir = "random";
  string energy_ref_name = "energy";
  string force_ref_name = "force";
  string torque_ref_name = "torque";

  for (int order = 1; order <=3; order++) {
      qm::EFTCalculator calculator(order);
      calculator.setup(database_file.c_str());

      string name;
      name = energy_ref_name + ".o" + to_string(order) + ".txt";
      ifstream ener_fs(name.c_str());
      name = force_ref_name + ".o" + to_string(order) + ".txt";
      ifstream force_fs(name.c_str());
      name = torque_ref_name + ".o" + to_string(order) + ".txt";
      ifstream torque_fs(name.c_str());

      auto dir = opendir(input_dir.c_str());
      assert(dir != nullptr);
      auto entity = readdir(dir);
      while (entity != nullptr) {
        int len = strlen(entity->d_name);
        if (entity->d_type == DT_REG && len >= 4 &&
            strcmp(entity->d_name + len - 4, ".inp") == 0) {
          vector<double> res = test_inp(input_dir + "/" + string(entity->d_name), calculator);
          vector<double> ref_res(7);
          string line;
          getline(ener_fs, line);
          ref_res[0] = stod(common::split(line)[1]);
          getline(force_fs, line);
          auto fields = common::split(line);
          for (size_t i = 0; i < 3; i++) {
              ref_res[i+1] = stod(fields[i+1]);
          }
          getline(torque_fs, line);
          fields = common::split(line);
          for (size_t i = 0; i < 3; i++) {
              ref_res[i+4] = stod(fields[i+1]);
          }

          for (size_t i = 0; i < 7; i++) {
              assert((abs(res[i] - ref_res[i])) < 1.E-7);
          }
        }
        entity = readdir(dir);
      }
      ener_fs.close();
      force_fs.close();
      torque_fs.close();
      cout  << "Passed corretness test for order "  << order << endl;
  }
  return 0;
}
