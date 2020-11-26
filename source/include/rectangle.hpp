#pragma once
#include <iostream>
#include <climits>
#include <algorithm>
#include <utility>
#include <vector>

template<size_t N>
class rectangle {

  typedef std::pair<float, float> interval;
  
  interval bounds[N];

 public:

  typedef interval* iterator;
  
  iterator begin();
  iterator end();
  
  rectangle& operator=(rectangle<N> &rect);
  interval& operator[](size_t idx);

  float get_area();
  void reset();
  void adjust(rectangle<N> &rect);
  
};

template<size_t N>
typename rectangle<N>::iterator rectangle<N>::begin() {
  return bounds;
}

template<size_t N>
typename rectangle<N>::iterator rectangle<N>::end() {
  return bounds + N;
}

template<size_t N>
rectangle<N>& rectangle<N>::operator=(rectangle<N> &rect) {
  std::copy(rect.begin(), rect.end(), begin());
  return *this;
}

template<size_t N>
typename rectangle<N>::interval& rectangle<N>::operator[](size_t idx) {
  return bounds[idx];
}

template<size_t N>
float rectangle<N>::get_area() {
  float area = 1;
  for(size_t i = 0; i < N; ++i) {
    area *= ((*this)[i].second - (*this)[i].first);
  }
  return area;
}

template<size_t N>
void rectangle<N>::reset() {
  for(size_t i = 0; i < N; ++i) {
    (*this)[i].first = LONG_MAX;
    (*this)[i].second = LONG_MIN;
  }
}

template<size_t N>
void rectangle<N>::adjust(rectangle<N> &rect) {
  for(size_t i = 0; i < N; ++i) {
    (*this)[i].first = std::min((*this)[i].first, rect[i].first);
    (*this)[i].second = std::max((*this)[i].second, rect[i].second);
  }
}

