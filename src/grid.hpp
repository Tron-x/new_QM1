/**
 * \file grid.hpp
 * \brief grid function
 */
#pragma once

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>

#include "common.h"

namespace qm {

using std::vector;
using std::array;

using namespace common;

typedef array<double, 7> data_t;


/**
* A structure for a grid that discretizes six dimensions one by one.
* The six dimensions are in the sequence of R, PHI, THETA for the
* translational space, followed by PHI1, PHI2, THETA for the rotational
* space
* represented by the quaternion q.
*/
struct Node {
  vector<double> xs;
  data_t y;
  vector<Node*> next;

  explicit Node(const vector<double>& _xs = {}) : xs(_xs) {}

  int count() {
    if (next.empty())
      return 1;
    else {
      int res = 0;
      for (auto node : next) {
        res += node->count();
      }
      return res;
    }
  }

};

class Grid {
 public:
  // Set the parameters that control the density of the grid points
  // For the meaning of ang_params and ori_params, check the _get_dl() method
  explicit Grid() {
    m_rs = {2.,         2.10714286, 2.21428571,  2.32142857,  2.42857143,
            2.53571429, 2.64285714, 2.75,        2.85714286,  2.96428571,
            3.07142857, 3.17857143, 3.28571429,  3.39285714,  3.5,
            3.8125,     4.125,      4.4375,      4.75,        5.0625,
            5.375,      5.6875,     6.,          6.85714286,  7.71428571,
            8.57142857, 9.42857143, 10.28571429, 11.14285714, 12.};

    m_ang_params = {0.20, 0.45, 3.8, 0.4};
    m_ori_params = {0.23, 0.45, 3.8, 0.4};
  }

  ~Grid() {
    // FIXME memory leak
    //free_node(m_head_node);
  }

  // Setup the grid structure without assigning the y values
  void setup() {
    Node* node_r = new Node(m_rs);

    for (double r : m_rs) {
      vector<double> phis = discretize_phi(r);
      Node* node_phi = new Node(phis);
      node_r->next.push_back(node_phi);
      for (double phi : phis) {
        vector<double> thetas = discretize_theta(r, phi);
        Node* node_theta = new Node(thetas);
        node_phi->next.push_back(node_theta);

        for (double theta : thetas) {
          vector<double> ophi1s = discretize_ophi1(r);
          Node* node_ophi1 = new Node(ophi1s);
          node_theta->next.push_back(node_ophi1);

          for (double ophi1 : ophi1s) {
            vector<double> ophi2s = discretize_ophi2(r, ophi1);
            Node* node_ophi2 = new Node(ophi2s);
            node_ophi1->next.push_back(node_ophi2);

            for (double ophi2 : ophi2s) {
              vector<double> othetas = discretize_otheta(r, ophi1, ophi2);
              Node* node_otheta = new Node(othetas);
              node_ophi2->next.push_back(node_otheta);

              for (double otheta : othetas) {
                node_otheta->next.push_back(new Node());
              }
            }
          }
        }
      }
    }

    m_head_node = node_r;
    m_n = node_r->count();
    std::cout << "The grid consists of "  << m_n << " points" << std::endl;
  }

  void load(const char* filename) {
    printf("Loading y values for each grid from file\n");

    std::ifstream ifs(filename);
    if (!ifs.good()) {
      std::cerr << "ERROR: Opening file '" << filename << "' failed.\n";
      exit(EXIT_FAILURE);
    }

    std::string line;
    if (!ifs.eof()) {
      //! Read first 3 lines
      //! First is ANGPRM
      //! Second is ORIPRM
      //! Third is RS
      std::getline(ifs, line);
      auto fields = common::split(line);
      assert(fields[0] == "ANGPRM" && fields.size() == 5);
      for (size_t idx = 1; idx < fields.size(); ++idx) {
        m_ang_params[idx-1] = std::stod(fields[idx]);
      }

      std::getline(ifs, line);
      fields = common::split(line);
      assert(fields[0] == "ORIPRM" && fields.size() == 5);
      for (size_t idx = 1; idx < fields.size(); ++idx) {
        m_ori_params[idx-1] = std::stod(fields[idx]);
      }

      std::getline(ifs, line);
      fields = common::split(line);
      assert(fields[0] == "RS");
      m_rs.clear();
      for (size_t idx = 1; idx < fields.size(); ++idx) {
        m_rs.push_back(std::stod(fields[idx]));
      }

      setup();

      fill_grid(m_head_node, ifs);
    }
  }

  data_t interpolate(const vector<double>& coor, int order = 2) {
    return interpolate_help(coor.begin(), m_head_node, order);
  }

  int get_n() const { return m_n; }

 private:
  void free_node(Node* head) {
    if (head == nullptr) return;
    for (Node* node : head->next) {
      free_node(node);
    }
    delete head;
  }

  void fill_grid(Node* head, std::ifstream& ifs) {
    if (head->next.empty()) {
      std::string line;
      std::getline(ifs, line);
      if (line == "") {
        return;
      }
      auto fields = common::split(line);
      assert(fields.size() == 7);
      for (size_t i = 0; i < 7; i++) {
        head->y[i] = std::stod(fields[i]);
      }
    }
    int idx = 0;
    for (Node* node : head->next) {
      fill_grid(node, ifs);
      idx += 1;
    }
  }

  data_t interpolate_help(vector<double>::const_iterator icoor, Node* node,
                                  int order = 2) {
    if (node->next.empty()) {
      return node->y;
    }

    double my_x = *icoor++;

    array<int, 4> neighbors;
    int n_nbr;
    switch (order) {
      case 1:
        n_nbr = find_neighbors2(node->xs, my_x, neighbors);
        break;
      case 2:
        n_nbr = find_neighbors3(node->xs, my_x, neighbors);
        break;
      case 3:
        n_nbr = find_neighbors4(node->xs, my_x, neighbors);
        break;
      default:
        QMException("Invalid order for interpolate . Choose from 1, 2, and 3");
        break;
    }

    array<double, 4> xs;
    array<data_t, 4> ys;
    for (size_t i = 0; i < n_nbr; i++) {
      xs[i] = node->xs[neighbors[i]];
      ys[i] = interpolate_help(icoor, node->next[neighbors[i]], order);
    }

    data_t my_y = interp_1D(xs, ys, my_x, n_nbr);
    return my_y;
  }

  // Two points for linear interpolation
  int find_neighbors2(const vector<double>& xs, double my_x, array<int, 4>& neighbors) {
    int n = xs.size();
    if (n <= 2) {
      for (int i = 0; i < n; i++) {
        neighbors[i] = i;
      }
      return n;
    }

    int idx = 1;
    while (idx < n) {
      if (xs[idx] > my_x) break;
      ++idx;
    }

    if (idx == n) {
      QMException("X value out of range");
    }
    neighbors[0] = idx - 1;
    neighbors[1] = idx;
    //neighbors = {idx - 1, idx};
    return 2;
  }

  // Three points for 2nd order interpolation
  int find_neighbors3(const vector<double>& xs, double my_x, array<int, 4>& neighbors) {
    int n = xs.size();
    if (n <= 3) {
      for (int i = 0; i < n; i++) {
        neighbors[i] = i;
      }
      return n;
    }

    int idx = 1;
    while (idx < n) {
      if (xs[idx] > my_x) break;
      ++idx;
    }

    if (idx == n) {
      QMException("X value out of range");
    }
    if (idx > n - 2) {
      neighbors[0] = idx - 2;
      neighbors[1] = idx - 1;
      neighbors[2] = idx;
    }
    else {
    neighbors[0] = idx - 1;
    neighbors[1] = idx;
    neighbors[2] = idx + 1;
    }
    return 3;
  }

  // Four points for 3rd order interpolation
  int find_neighbors4(const vector<double>& xs, double my_x, array<int, 4>& neighbors) {
    int n = xs.size();
    if (n < 4) {
      for (int i = 0; i < n; i++) {
        neighbors[i] = i;
      }
      return n;
    }

    int idx = 1;
    while (idx < n) {
      if (xs[idx] > my_x) break;
      ++idx;
    }

    if (idx == n) {
      QMException("X value out of range");
    }
    if (idx == 1) {
      neighbors[0] = 0;
      neighbors[1] = 1;
      neighbors[2] = 2;
      return 3;
    }
    if (idx == n - 1) {
      neighbors[0] = n - 3;
      neighbors[1] = n - 2;
      neighbors[2] = n - 1;
      return 3;
    }
    neighbors[0] = idx - 2;
    neighbors[1] = idx - 1;
    neighbors[2] = idx;
    neighbors[3] = idx + 1;
    return 4;
  }

  data_t interp_1D(const array<double, 4>& xs,
                           const array<data_t, 4>& ys, double my_x,
                           size_t size) {
    if (size == 1) return ys[0];
    if (size == 2) {
      double x0 = xs[0], x1 = xs[1];
      const data_t y0 = ys[0];
      const data_t y1 = ys[1];
      double x01 = x1 - x0;
      double m0 = my_x - x0;
      data_t rst;
      for (size_t i = 0; i < 7; i++) {
        rst[i] = y0[i] + m0 * (y1[i] - y0[i]) / x01;
      }
      return rst;
    }
    if (size == 3) {
      double x0 = xs[0], x1 = xs[1], x2 = xs[2];
      const data_t y0 = ys[0];
      const data_t y1 = ys[1];
      const data_t y2 = ys[2];
      data_t rst;
      double x01 = x1 - x0;
      double x02 = x2 - x0;
      double x12 = x2 - x1;
      double m0 = my_x - x0;
      double m1 = my_x - x1;
      double a1, a2;
      for (size_t i = 0; i < 7; i++) {
          a1 = (y1[i] - y0[i]) / x01;
          a2 = (y2[i] - y0[i] - a1 * x02) / x02 / x12;
          rst[i] = y0[i] + a1 * m0 + a2 * m0 * m1;
      }
      return rst;
    }

    if (size == 4) {
      double x0 = xs[0], x1 = xs[1], x2 = xs[2], x3 = xs[3];
      const data_t y0 = ys[0];
      const data_t y1 = ys[1];
      const data_t y2 = ys[2];
      const data_t y3 = ys[3];
      data_t rst;
      double x01 = x1 - x0;
      double x02 = x2 - x0;
      double x03 = x3 - x0;
      double x12 = x2 - x1;
      double x13 = x3 - x1;
      double x23 = x3 - x2;
      double m0 = my_x - x0;
      double m1 = my_x - x1;
      double m2 = my_x - x2;
      double a1, a2, a3;
      for (size_t i = 0; i < 7; i++) {
        a1 = (y1[i] - y0[i]) / x01;
        a2 = (y2[i] - y0[i] - a1 * x02) / x02 / x12;
        a3 = (y3[i] - y0[i] - a1 * x03 - a2 * x03 * x13) / x03 / x13 / x23;
        rst[i] = y0[i] + a1 * m0 + a2 * m0 * m1 + a3 * m0 * m1 * m2;
      }
      return rst;
    }
  }
  /**
   * The following methods discretizes each dimension.
   * The density depends only on r.
   * attention: the discretizing range for each dimension fits the reflection
   * planes of water's refCoor, need to consider propers ways to generalize
   */
  vector<double> discretize_phi(double r) {
    double dl = get_dl(r, m_ang_params);
    double mini = 0;
    double maxi = M_PI / 2;
    double lvalue = maxi - mini;
    int n = std::max(1, static_cast<int>(std::ceil(lvalue / dl)));

    return common::linspace(mini, maxi, n);
  }

  vector<double> discretize_theta(double r, double phi) {
    double dl = get_dl(r, m_ang_params);
    double mini = 0;
    double maxi = M_PI;
    double lvalue = (maxi - mini) * std::cos(phi);
    int n = std::max(1, static_cast<int>(std::ceil(lvalue / dl)));

    return common::linspace(mini, maxi, n);
  }

  vector<double> discretize_ophi1(double r) {
    double dl = get_dl(r, m_ori_params);
    double mini = 0;
    double maxi = M_PI / 2;
    double lvalue = maxi - mini;
    int n = std::max(1, static_cast<int>(std::ceil(lvalue / dl)));

    return common::linspace(mini, maxi, n);
  }

  vector<double> discretize_ophi2(double r, double ophi1) {
    double dl = get_dl(r, m_ori_params);
    double mini = 0;
    double maxi = M_PI / 2;
    double lvalue = (maxi - mini) * std::cos(ophi1);
    int n = std::max(1, static_cast<int>(std::ceil(lvalue / dl)));

    return common::linspace(mini, maxi, n);
  }

  vector<double> discretize_otheta(double r, double ophi1, double ophi2) {
    double dl = get_dl(r, m_ori_params);
    double mini = -M_PI;
    double maxi = M_PI;
    double lvalue = (maxi - mini) * std::cos(ophi1) * std::cos(ophi2);
    int n = std::max(1, static_cast<int>(std::ceil(lvalue / dl)));

    return common::linspace(mini, maxi, n);
  }

  /**
   * Use a 4-parameter sigmoidal function to calculate the approximate distance
   * between two points, on the 3D sphere for angular DOFs or 4D sphere for
   * rotational DOFs.The four parameters specify the lower limit, upper limit,
   * transition position and transition steepness of the sigmoidal function
   */
  double get_dl(double r, const array<double, 4>& params) {
    return (params[1] - params[0]) /
               (1 + std::exp((params[2] - r) / params[3])) +
           params[0];
  }

  Node* m_head_node = nullptr;
  int m_n = 0;
  vector<double> m_rs;
  array<double, 4> m_ang_params;
  array<double, 4> m_ori_params;
};

}  // qm
