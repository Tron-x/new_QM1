/**
 * \file common_test.cpp
 * \brief Contain the test functions for common.h functions
 */

#include <assert.h>
#include <cmath>
#include <vector>

#include "common.h"

using namespace std;

#define EPS 1e-6
template<typename T>
void assert_equal(T actual, T expect) {
  assert(actual == expect);
}

template<>
void assert_equal<double>(double actual, double expect) {
  assert(std::abs(actual - expect) <= EPS);
}
template<>
void assert_equal<float>(float actual, float expect) {
  assert(std::abs(actual - expect) <= EPS);
}

template<typename T>
void assert_equal(vector<T>&& actual, vector<T>&& expect) {
  assert(actual.size() == expect.size());
  for (size_t idx = 0; idx < actual.size(); ++idx) {
    assert_equal(actual[idx], expect[idx]);
  }
}

void test_linspace() {
  vector<double> expect = { 1.3,  1.9,  2.5};
  assert_equal(common::linspace(1.3, 2.5, 3), expect);
}

void test_vecop() {
  vector<double> expect = { 1.3,  1.9,  2.5};
  vector<double> a = {0.3, 0.9, 1.5};
  vector<double> b = {1.0, 1.0, 1.0};
  assert_equal(a + 1.0, expect);
  assert_equal(a + b, expect);
}

void test_tokenizer() {
  using namespace common;
  string s = "12  32  4354";
  vector<string> tokens = tokenizer(s, ' ');
  vector<string> expect = {"12", "32", "4354"};

  assert_equal(tokens, expect);
}

int main() {
  test_linspace();
  test_vecop();
  test_tokenizer();

  printf("Test passed!!!\n");
  return 0;
}
