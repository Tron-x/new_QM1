#include <iostream>
#include <getopt.h>
#include <dirent.h>
#include <cstring>
#include <time.h>

#include "eft_calculator.hpp"
#include "test.h"

using std::string;


static void usage() {
  using std::cerr;
  using std::endl;
  cerr << endl;
  cerr << "About: energy force calculation, current only support waterbox"
       << endl;
  cerr << "Usage: qm [option]" << endl;
  cerr << "Options:" << endl;
  cerr << "  -o, --order  the order of the interpolant" << endl;
  cerr << "  -d, --database  the database file(gz file)" << endl;
  cerr << "  -i, --input     input dir which contains .inp files" << endl;
  cerr << "  -e, --ouput_energy     output file" << endl;
  cerr << "  -f, --ouput_force     output file" << endl;
  cerr << "  -t, --ouput_torque     output file" << endl;

  cerr << "Example:" << endl;
  cerr << "  qm -o 2 -d grid.dat -i random -e energy.txt -f force.txt -t torque.txt"
       << endl;

  exit(1);
}

int main(int argc, char* argv[]) {
  int order;
  string database_file;
  string input_dir;
  string output_energy_file;
  string output_force_file;
  string output_torque_file;
  string model;
  long start,finish;
  double duration;

  int c;
  static struct option loptions[] = {
      {"help", no_argument, NULL, 'h'},
      {"order", required_argument, NULL, 'o'},
      {"database", required_argument, NULL, 'd'},
      {"input", required_argument, NULL, 'i'},
      {"ouput_energy", required_argument, NULL, 'e'},
      {"ouput_force", required_argument, NULL, 'f'},
      {"ouput_torque", required_argument, NULL, 't'},
      {"ouput", required_argument, NULL, 'o'},
      {NULL, 0, NULL, 0}};

  while ((c = getopt_long(argc, argv, "h?o:d:i:e:f:t:?", loptions, NULL)) != -1) {
    switch (c) {
      case 'o':
        order = std::stoi(optarg);
      case 'd':
        database_file = optarg;
        break;
      case 'i':
        input_dir = optarg;
        break;
      case 'e':
        output_energy_file = optarg;
        break;
      case 'f':
        output_force_file = optarg;
        break;
      case 't':
        output_torque_file = optarg;
        break;
      case 'm':
        model = optarg;
        break;
      default:
        usage();
    }
  }
  if (optind < 8)  // 4 args and qm
  {
    std::cerr << "Error: you input less argv, see:" << std::endl;
    usage();
  }

  qm::EFTCalculator calculator(order);
  calculator.setup(database_file.c_str());
  start = clock();

  auto dir = opendir(input_dir.c_str());
  assert(dir != nullptr);

  FILE* fp_energy_out = fopen(output_energy_file.c_str(), "w");
  FILE* fp_force_out = fopen(output_force_file.c_str(), "w");
  FILE* fp_torque_out = fopen(output_torque_file.c_str(), "w");
  auto entity = readdir(dir);
  while (entity != nullptr) {
    int len = std::strlen(entity->d_name);
    if (entity->d_type == DT_REG && len >= 4 &&
        strcmp(entity->d_name + len - 4, ".inp") == 0) {
      string filename = input_dir + "/" + std::string(entity->d_name);
      auto eft = test_inp(filename, calculator);
      fprintf(fp_energy_out, "%s %12.7f\n", filename.c_str(), eft[0]);
      fprintf(fp_force_out, "%s %12.7f %12.7f %12.7f\n", filename.c_str(), eft[1], eft[2],
              eft[3]);
      fprintf(fp_torque_out, "%s %12.7f %12.7f %12.7f\n", filename.c_str(), eft[4], eft[5],
              eft[6]);
    }

    entity = readdir(dir);
  }

  fclose(fp_energy_out);
  fclose(fp_force_out);
  fclose(fp_torque_out);
  finish =clock();
  duration = (double)(finish-start) / CLOCKS_PER_SEC;
  std::cout<<"the computation time is : "<<duration<<std::endl;
  return 0;
}
