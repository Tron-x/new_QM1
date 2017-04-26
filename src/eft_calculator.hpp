/**
 * \file grid.hpp
 * \brief grid function
 */
#pragma once

#include "common.h"
#include "grid.hpp"
#include "tools.h"

#include <vector>
using std::vector;

namespace qm {

/**
 * A class that holds information related to the atomic structure of a water
 * molecule. It also includes several methods that carries out operations
 * related to the atomic coordinates.
 */
class Water {
 public:
  explicit Water() {
    m_mass = {15.99900, 1.00800, 1.00800};
    m_ref_coor = {{-0.06556939, 0., 0.},
                  {0.52035943, -0.76114632, 0.},
                  {0.52035943, 0.76114632, 0.}};
    m_refl_axes = {1, 2};
  }

  vector<vector<double>> get_R(vector<vector<double>> coors) {
    assert(coors.size() == 3);
    auto x = coors[0];
    for (auto& coor : coors) {
      coor = coor - x;
    }

    vector<double> xvec = coors[1] + coors[2];
    vector<double> zvec = common::cross(coors[1], coors[2]);
    vector<double> yvec = common::cross(zvec, xvec);
    xvec = xvec / common::norm(xvec);
    yvec = yvec / common::norm(yvec);
    zvec = zvec / common::norm(zvec);

    //! transpose [xvec, yvec, zvec]
    vector<vector<double>> res;
    for (size_t idx = 0; idx < xvec.size(); ++idx) {
      res.push_back(vector<double>{xvec[idx], yvec[idx], zvec[idx]});
    }
    return res;
  }

  vector<double> get_com(const vector<vector<double>>& coors) {
    vector<double> res = {0, 0, 0};
    double sum_mass = sum(m_mass);
    for (size_t i = 0; i < coors.size(); ++i) {
      for (int idx = 0; idx < 3; ++idx) {
        res[idx] += coors[i][idx] * m_mass[i];
      }
    }
    for (int idx = 0; idx < 3; ++idx) {
      res[idx] /= sum_mass;
    }
    return res;
  }

  vector<int>& get_ref_axes() {
    return m_refl_axes;
  }

 private:
  vector<double> m_mass;
  vector<vector<double>> m_ref_coor;
  vector<int> m_refl_axes;
};

/**
 * A class that carries the logic of evaluating the energy, force and torque
 * of a pair of rigid molecules. The coordinates of each molecule are given
 * in the form of Xcom and q, with Xcom being the Cartesian coordinates of the
 * center of mass, q being the quaternion representation of its orientation
 * wrt a reference pose. The class evaluates the EFTs for a pair of such
 * coordinates by
 *   1. Apply translational and rotational operations to the pair to align the
 *      COM of the first molecule with the origin and its orientation to the
 *      reference pose.
 *   2. Convert the modified Xcom and q of the second molecule into spherical
 *      coordinates.
 *   3. Use the resulted six-dimensional coordinate to query a six-dimensional
 *      grid that stores precomputed EFTs.
 *   4. Unapply rotation in step 1 to obtain correctly oriented forces and
 * torques
 */
class EFTCalculator {
 public:
  explicit EFTCalculator(int order = 2) : m_order(order) {}

  void setup(const char* filename = nullptr) {
    if (filename) {
      m_grid.load(filename);
    } else {
      m_grid.setup();
    }
  }

  vector<double> eval(const vector<vector<double>>& coors) {
    vector<vector<double>> coors0(coors.begin(), coors.begin() + 3);
    vector<vector<double>> coors1(coors.begin() + 3, coors.end());
    vector<double> x_com0 = m_mol.get_com(coors0);
    vector<double> x_com1 = m_mol.get_com(coors1);

    vector<vector<double>> R0 = m_mol.get_R(coors0);
    vector<double> q0 = tools::R2q(R0);
    vector<vector<double>> R1 = m_mol.get_R(coors1);
    vector<double> q1 = tools::R2q(R1);

    //! move COM if mol0 to origin
    vector<double> X = dot(x_com1 - x_com0, R0);
    vector<double> q = tools::qdiv(q1, q0);
    //! Use mirror symmetry of mol0 to move mol1 such that its COM has positive y and z values
    vector<int> reflections;

    //! reserve q[0] value
    double tmp = q[0];
    for (int i : m_mol.get_ref_axes()) {
      if (X[i] < 0) {
        X[i] = -X[i];
        q[i+1] = -q[i+1];
        q = q * -1.0;
        reflections.push_back(i);
      }
    }
    q[0] = tmp;

    if (q[0] < 0) {
      q = q * -1.0;
    }
    if (q[1] < 0) {
      q = {-q[1], q[0], q[3], -q[2]};
    }
    auto xyz2sperical = tools::xyz2spherical(X);
    auto q2spherical = tools::q2spherical(q);

    vector<double> coor;
    coor.insert(coor.end(), xyz2sperical.begin(), xyz2sperical.end());
    coor.insert(coor.end(), q2spherical.begin(), q2spherical.end());

    // use the grid to obtain results
    // res.size() == 7
    // energy -> 0
    // force -> 1:4
    // torque -> 4:7
    data_t res = m_grid.interpolate(coor, m_order);

    double energy = res[0];
    vector<double> force(res.begin() + 1, res.begin()+4);
    vector<double> torque(res.begin()+4, res.end());
    for (auto i : reflections) {
      force[i] = -force[i];
      torque[i] = -torque[i];
      torque = torque * -1.;
    }

    // reverse the reorientation applied to align mol0 with refcoor
    vector<double> new_force;
    vector<double> new_torque;
    for (const vector<double>& item : R0) {
      new_force.push_back(common::dot(force, item));
      new_torque.push_back(common::dot(torque, item));
    }
    return {energy, new_force[0], new_force[1], new_force[2], new_torque[0], new_torque[1], new_torque[2]};
  }

 private:
  int m_order;
  Water m_mol;
  Grid m_grid;
};

class QMInterpolation {
 public:
  explicit QMInterpolation(EFTCalculator& calculator)
      : m_calculator(calculator) {}

  vector<string> process(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.good()) {
      std::cerr << "ERROR: Opening file '" << filename << "' failed.\n";
      exit(EXIT_FAILURE);
    }

    vector<string> lines;
    string line;
    while (!ifs.eof()) {
      std::getline(ifs, line);
      common::trim(line);

      if (line == "" || line == "$END" || line == "$end") {
        break;
      }
      lines.push_back(line);
    }

    vector<vector<double>> coors;
    for (int idx = lines.size() - 6; idx < lines.size(); idx++) {
      auto fields = common::tokenizer(lines[idx], ' ');
      assert(fields.size() == 5);
      coors.push_back(
          {std::stod(fields[2]), std::stod(fields[3]), std::stod(fields[4])});
    }
    ifs.close();

    auto eft = m_calculator.eval(coors);
    std::vector<std::string> result;
    char res[100];
    sprintf(res, "%s %12.7f", filename.c_str(), eft[0]);
    result.push_back(res);
    sprintf(res, "%s %12.7f %12.7f %12.7f", filename.c_str(), eft[1], eft[2],
            eft[3]);
    result.push_back(res);
    sprintf(res, "%s %12.7f %12.7f %12.7f", filename.c_str(), eft[4], eft[5],
            eft[6]);
    result.push_back(res);

    return result;
  }

 private:
  EFTCalculator m_calculator;
};

}  // qm
