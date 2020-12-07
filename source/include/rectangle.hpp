#pragma once
#include <iostream>
#include <climits>
#include <algorithm>
#include <utility>

template<size_t N>
class Rectangle {

  typedef std::pair<float, float> interval;
  
  interval bounds[N];

 public:

  typedef interval* iterator;
  
  iterator begin();
  iterator end();
  
  Rectangle& operator=(Rectangle<N> &rect);
  interval& operator[](size_t idx);

  float get_area();
  void reset();
  void adjust(Rectangle<N> &rect);
  float get_overlap(Rectangle<N> &rect);
  
};

template<size_t N>
typename Rectangle<N>::iterator Rectangle<N>::begin() {
  return bounds;
}

template<size_t N>
typename Rectangle<N>::iterator Rectangle<N>::end() {
  return bounds + N;
}

template<size_t N>
Rectangle<N>& Rectangle<N>::operator=(Rectangle<N> &rect) {
  std::copy(rect.begin(), rect.end(), begin());
  return *this;
}

template<size_t N>
typename Rectangle<N>::interval& Rectangle<N>::operator[](size_t idx) {
  return bounds[idx];
}

template<size_t N>
float Rectangle<N>::get_area() {
  float area = 1;
  for(size_t i = 0; i < N; ++i) {
    area *= ((*this)[i].second - (*this)[i].first);
  }
  return area;
}

template<size_t N>
void Rectangle<N>::reset() {
  for(size_t i = 0; i < N; ++i) {
    (*this)[i].first = LONG_MAX;
    (*this)[i].second = LONG_MIN;
  }
}

template<size_t N>
void Rectangle<N>::adjust(Rectangle<N> &rect) {
  for(size_t i = 0; i < N; ++i) {
    (*this)[i].first = std::min((*this)[i].first, rect[i].first);
    (*this)[i].second = std::max((*this)[i].second, rect[i].second);
  }
}

template<size_t N>
float Rectangle<N>::get_overlap(Rectangle<N> &rect) {
  float area = 1;
  for(size_t i = 0; i < N; ++i) {
    float left = std::max((*this)[i].first, rect[i].first);
    float right = std::min((*this)[i].second, rect[i].second);
    area *= std::max(float(0), right - left);
  }
  return area;
}
