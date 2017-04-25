/**
 * \file common.h
 * \brief contain the helper functions
 */
#pragma once

#include "common.h"

#include <vector>
#include <array>
#include <cmath>

using std::vector;
using std::array;

using namespace common;

namespace tools {

static inline int sign(double v) {
  if (v == 0) return 0;
  if (v < 0) return -1;
  return 1;
}

static inline vector<double> qmult(const vector<double>& a,
                                   const vector<double>& b) {
  assert(a.size() == b.size());
  assert(a.size() == 4);
  vector<double> q = {0, 0, 0, 0};
  q[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  q[1] = a[0] * b[1] + a[1] * b[0] - a[2] * b[3] + a[3] * b[2];
  q[2] = a[0] * b[2] + a[2] * b[0] + a[1] * b[3] - a[3] * b[1];
  q[3] = a[0] * b[3] + a[3] * b[0] - a[1] * b[2] + a[2] * b[1];
  return q;
}

static inline vector<double> qinv(const vector<double>& q) {
  vector<double> res = q;
  res[0] = -res[0];
  return res;
}


static inline vector<double> qdiv(const vector<double>& a,
                                   const vector<double>& b) {
  assert(a.size() == b.size());
  assert(a.size() == 4);
  return qmult(a, qinv(b));
}

vector<double> R2q(const vector<vector<double>>& R) {
  double q0 = (R[0][0] + R[1][1] + R[2][2] + 1) / 4;
  double q1 = (R[0][0] - R[1][1] - R[2][2] + 1) / 4;
  double q2 = (-R[0][0] + R[1][1] - R[2][2] + 1) / 4;
  double q3 = (-R[0][0] - R[1][1] + R[2][2] + 1) / 4;

  q0 = std::max(q0, 0.);
  q1 = std::max(q1, 0.);
  q2 = std::max(q2, 0.);
  q3 = std::max(q3, 0.);
  q0 = std::sqrt(q0);
  q1 = std::sqrt(q1);
  q2 = std::sqrt(q2);
  q3 = std::sqrt(q3);

  if (q0 >= q1 && q0 >= q2 && q0 >= q3) {
    q0 *= 1;
    q1 *= sign(R[2][1] - R[1][2]);
    q2 *= sign(R[0][2] - R[2][0]);
    q3 *= sign(R[1][0] - R[0][1]);
  } else if (q1 >= q0 && q1 >= q2 && q1 >= q3) {
    q0 *= sign(R[2][1] - R[1][2]);
    q1 *= 1;
    q2 *= sign(R[1][0] + R[0][1]);
    q3 *= sign(R[0][2] + R[2][0]);
  } else if (q2 >= q0 && q2 >= q1 && q2 >= q3) {
    q0 *= sign(R[0][2] - R[2][0]);
    q1 *= sign(R[1][0] + R[0][1]);
    q2 *= 1;
    q3 *= sign(R[2][1] + R[1][2]);
  } else if (q3 >= q0 && q3 >= q1 && q3 >= q2) {
    q0 *= sign(R[1][0] - R[0][1]);
    q1 *= sign(R[2][0] + R[0][2]);
    q2 *= sign(R[2][1] + R[1][2]);
    q3 *= 1;
  } else {
    QMException("coding error");
  }

  vector<double> res = {q0, q1, q2, q3};
  res = res / common::norm(res);
  return res;
}

array<double, 3> xyz2spherical(const vector<double>& X) {
  double x = X[0], y = X[1], z = X[2];
  double r = std::sqrt(x * x + y * y + z * z);
  double phi = std::atan2(z, std::sqrt(x * x + y * y));
  double theta = std::atan2(y, x);

  return {r, phi, theta};
}

array<double, 3> q2spherical(vector<double>& q) {
  if (q[0] < 0) q = q * -1.;

  double q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
  double theta = std::atan2(q2, q3);
  double phi2 = std::atan2(q1, std::sqrt(q2 * q2 + q3 * q3));
  double phi1 = std::atan2(q0, std::sqrt(q1 * q1 + q2 * q2 + q3 * q3));

  return {phi1, phi2, theta};
}

}  // tools
