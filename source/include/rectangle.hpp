#pragma once
#include <iostream>
#include <algorithm>
#include <utility>
#include <limits>
#include <vector>

template<size_t N>
class Rectangle {

  typedef std::pair<float, float> interval;
  
  interval bounds[N];

 public:

  typedef interval* iterator;
  typedef const interval* const_iterator;
  
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  
  Rectangle& operator=(const Rectangle<N> &rect);
  interval& operator[](size_t idx);
  interval operator[](size_t idx) const;

  float get_area();
  float get_margin();
  void reset();
  void adjust(const Rectangle<N> &rect);
  float get_overlap(const Rectangle<N> &rect);
  
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
typename Rectangle<N>::const_iterator Rectangle<N>::begin() const {
  return bounds;
}

template<size_t N>
typename Rectangle<N>::const_iterator Rectangle<N>::end() const {
  return bounds + N;
}

template<size_t N>
Rectangle<N>& Rectangle<N>::operator=(const Rectangle<N> &rect) {
  std::copy(rect.begin(), rect.end(), begin());
  return *this;
}

template<size_t N>
typename Rectangle<N>::interval& Rectangle<N>::operator[](size_t idx) {
  return bounds[idx];
}

template<size_t N>
typename Rectangle<N>::interval Rectangle<N>::operator[](size_t idx) const {
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
float Rectangle<N>::get_margin() {
  float margin = 0;
  for(size_t i = 0; i < N; ++i) {
    margin += ((*this)[i].second - (*this)[i].first);
  }
  margin *= (1 << (N - 1));
  return margin;
}

template<size_t N>
void Rectangle<N>::reset() {
  for(size_t i = 0; i < N; ++i) {
    (*this)[i].first = std::numeric_limits<float>::max();
    (*this)[i].second = std::numeric_limits<float>::min();
  }
}

template<size_t N>
void Rectangle<N>::adjust(const Rectangle<N> &rect) {
  for(size_t i = 0; i < N; ++i) {
    (*this)[i].first = std::min((*this)[i].first, rect[i].first);
    (*this)[i].second = std::max((*this)[i].second, rect[i].second);
  }
}

template<size_t N>
float Rectangle<N>::get_overlap(const Rectangle<N> &rect) {
  float area = 1;
  for(size_t i = 0; i < N; ++i) {
    float left = std::max((*this)[i].first, rect[i].first);
    float right = std::min((*this)[i].second, rect[i].second);
    area *= std::max(float(0), right - left);
  }
  return area;
}
