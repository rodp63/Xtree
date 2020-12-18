#pragma once

#include <utility>
#include <limits>
#include <vector>

template<size_t N>
class Point;

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

  float MINDIST(const Point<N> &pt);
  float MINMAXDIST(const Point<N> &pt);
  
};

template<size_t N>
class Point {

  float coords[N];

 public:

  typedef float* iterator;
  typedef const float* const_iterator;
  
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

  Point& operator=(const Point<N> &rect);
  float& operator[](size_t idx);
  float operator[](size_t idx) const;

  Rectangle<N> get_rect() const;
  std::vector<float> get_vector();
};

// Rectangle Implementation

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
    if ((*this)[i].second != (*this)[i].first)
      area *= ((*this)[i].second - (*this)[i].first);
  }
  return area == 1 ? 0 : area;
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
    if (right != left)
      area *= std::max(float(0), right - left);
  }
  return area == 1 ? 0 : area;
}

template<size_t N>
float Rectangle<N>::MINDIST(const Point<N> &pt) {
  float ans = 0, dist;
  for (size_t i = 0; i < N; ++i) {
    if (pt[i] < bounds[i].first) {
      dist = pt[i] - bounds[i].first;
    }
    else if (pt[i] > bounds[i].second) {
      dist = pt[i] - bounds[i].second;
    }
    else {
      dist = 0;
    }
    ans += dist * dist;
  }
  return ans;
}

template<size_t N>
float Rectangle<N>::MINMAXDIST(const Point<N> &pt) {
  float S = 0, dist;
  for (size_t i = 0; i < N; ++i) {
    if (pt[i] >= (bounds[i].first + bounds[i].second) / 2) {
      dist = pt[i] - bounds[i].first;
    }
    else {
      dist = pt[i] - bounds[i].second;
    }
    S += dist * dist;
  }
  float ans = std::numeric_limits<float>::max(), cur;
  for (size_t i = 0; i < N; ++i) {
    if (pt[i] >= (bounds[i].first + bounds[i].second) / 2) {
      dist = pt[i] - bounds[i].first;
    }
    else {
      dist = bounds[i].second - pt[i];
    }
    cur = S - (dist * dist);
    if (pt[i] <= (bounds[i].first + bounds[i].second) / 2) {
      dist = pt[i] - bounds[i].first;
    }
    else {
      dist = bounds[i].second - pt[i];
    }
    cur += dist * dist;
    ans = std::min(ans, cur);
  }
  return ans;
}

// Point Implementation

template<size_t N>
typename Point<N>::iterator Point<N>::begin() {
  return coords;
}

template<size_t N>
typename Point<N>::iterator Point<N>::end() {
  return coords + N;
}

template<size_t N>
typename Point<N>::const_iterator Point<N>::begin() const {
  return coords;
}

template<size_t N>
typename Point<N>::const_iterator Point<N>::end() const {
  return coords + N;
}

template<size_t N>
Point<N>& Point<N>::operator=(const Point<N> &p) {
  std::copy(p.begin(), p.end(), begin());
  return *this;
}

template<size_t N>
float& Point<N>::operator[](size_t idx) {
  return coords[idx];
}

template<size_t N>
float Point<N>::operator[](size_t idx) const {
  return coords[idx];
}

template<size_t N>
Rectangle<N> Point<N>::get_rect() const {
  Rectangle<N> r;
  for (size_t i = 0; i < N; ++i) {
    r[i].first = r[i].second = coords[i];
  }
  return r;
}

template<size_t N>
std::vector<float> Point<N>::get_vector() {
  std::vector<float> vect(N);
  std::copy(begin(), end(), vect.begin());
  return vect;
}
