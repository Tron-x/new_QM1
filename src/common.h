/**
 * \file common.h
 * \brief contain the helper functions
 */
#pragma once

#include <cstddef>
#include <cassert>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <string>
#include <sstream>
#include <iostream>

using std::vector;
using std::string;
using std::stringstream;

namespace common {

#define QMException(x)         \
  fprintf(stderr, "%s\n", #x); \
  exit(1);

vector<double> linspace(double min, double max, size_t n) {
  if (n == 1) return {min};
  double step = (max - min) / (n - 1);

  vector<double> res;
  double start = min;
  for (size_t idx = 0; idx < n; ++idx) {
    res.push_back(start);
    start += step;
  }
  return res;
}

vector<int> range(size_t n) {
  vector<int> res;
  for (size_t idx = 0; idx < n; ++idx) {
    res.push_back(idx);
  }
  return res;
}

double dot(const vector<double>& v1, const vector<double>& v2) {
  assert(v1.size() == v2.size());
  double res = 0.0;
  for (unsigned int i = 0; i < v1.size(); i++) {
    res += v1[i] * v2[i];
  }
  return res;
}

/**
 * XXX The dot function is not the combine of the dot function above
 */
vector<double> dot(const vector<double>& v1, const vector<vector<double>>& v2) {
  vector<double> res(v1.size(), 0);

  for (size_t i = 0; i < v2.size(); ++i) {
    for (size_t j = 0; j < v1.size(); ++j) {
      res[j] += v2[i][j] * v1[i];
    }
  }

  return res;
}


template <typename T>
T sum(const vector<T>& vec) {
  T res = 0;
  for (auto&& item : vec) {
    res += item;
  }
  return res;
}

/**
 * TODO: only for vector of length 3
 */
vector<double> get_normal(const vector<double>& vec1,
                          const vector<double>& vec2) {
  assert(vec1.size() == 3);
  assert(vec2.size() == 3);
  vector<double> res(3, 0.0);
  res[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  res[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  res[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
  return res;
}

double length(const vector<double>& vec)
{
  return sqrt(dot(vec,vec));
}

double norm(const vector<double>& vec) {
  return length(vec);
}

vector<double> cross(const vector<double>& vec1, const vector<double>& vec2) {
  return get_normal(vec1, vec2);
}

double angle(const vector<double>& vec1, const vector<double>& vec2) {
  double dotprod = dot(vec1, vec2);
  return acos(dotprod / ((length(vec1) * length(vec2))));
}

//////////////////////////////helper function////////////////
// trim from start (in place)
static inline void ltrim(std::string& s) {
  s.erase(s.begin(),
          std::find_if(s.begin(), s.end(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
          s.end());
}

// trim from both ends (in place)
static inline void trim(std::string& s) {
  ltrim(s);
  rtrim(s);
}

// trim from start (copying)
static inline std::string ltrimmed(std::string s) {
  ltrim(s);
  return s;
}

// trim from end (copying)
static inline std::string rtrimmed(std::string s) {
  rtrim(s);
  return s;
}

template <typename Out>
void split(string s, Out result) {
  string item;
  string::iterator pos;

  trim(s);
  while (s.size()) {
    pos = std::find_if(s.begin(), s.end(), std::ptr_fun<int, int>(std::isspace));
    item = string(s.begin(), pos);
    *(result++) = item;
    s = string(pos, s.end());
    ltrim(s);
  }
}

vector<string> split(const string& s) {
  vector<string> elems;
  split(s, back_inserter(elems));
  return elems;
}

vector<string> tokenizer(const string& s, char delim) {
  const char* str = s.c_str();
  vector<string> result;
  do {
    const char* begin = str;

    while (*str != delim && *str) str++;

    result.push_back(string(begin, str));

    // skip delim
    while (*str == delim) str++;
  } while (0 != *str);

  return result;
}


}  // common

namespace std {

////////////////////////// Vecop ///////////////////////
template <typename T>
vector<T> vec_op(const vector<T>& a, const vector<T>& b,
                 std::function<T(T, T)> op) {
  assert(a.size() == b.size());

  vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), op);
  return result;
}

template <typename T>
vector<T> vec_op(const vector<T>& a, T b, std::function<T(T, T)> op) {
  vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), std::back_inserter(result),
                 std::bind(op, std::placeholders::_1, b));
  return result;
}

#define VECOP(op, opfun)                                          \
  template <typename T>                                           \
  vector<T> operator op(const vector<T>& a, const vector<T>& b) { \
    return vec_op<T>(a, b, std::opfun<T>());                      \
  }                                                               \
  template <typename T>                                           \
  vector<T> operator op(const vector<T>& a, T b) {                \
    return vec_op<T>(a, b, std::opfun<T>());                      \
  }                                                               \
  template <typename T>                                           \
  vector<T> operator op(T b, const vector<T>& a) {                \
    return vec_op<T>(a, b, std::opfun<T>());                      \
  }

VECOP(+, plus)
VECOP(-, minus)
VECOP(*, multiplies)
VECOP(/, divides)

} // std
