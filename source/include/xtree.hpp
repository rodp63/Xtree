#pragma once

#include <vector>
#include "rectangle.hpp"

using namespace std;

template<typename T, size_t N, size_t M, size_t m = M/2>
class xtree {
  
 public:

  struct node;
  
  struct cell {
    
    rectangle<N> MBR;
    std::shared_ptr<node> child;
    std::shared_ptr<T> data;
    
  };

  struct node {

    typedef typename std::vector<cell>::iterator iterator;

    iterator begin();
    iterator end();

    cell& operator[](size_t index);

    bool is_leaf();
    
    std::vector<cell> entries;
    size_t max_size;
    
  };
  

  xtree();
  ~xtree();
  size_t dimension() const;
  size_t size() const;
  bool empty() const;
  

 private:

  size_t entries;
  std::shared_ptr<node> root;
  
};

template<typename T, size_t N, size_t M, size_t m>
typename xtree<T, N, M, m>::node::iterator
xtree<T, N, M, m>::node::begin() {
  return entries.begin();
}

template<typename T, size_t N, size_t M, size_t m>
typename xtree<T, N, M, m>::node::iterator
xtree<T, N, M, m>::node::end() {
  return entries.end();
}

template<typename T, size_t N, size_t M, size_t m>
xtree<T, N, M, m>::xtree() {
  
}

template<typename T, size_t N, size_t M, size_t m>
xtree<T, N, M, m>::~xtree() {
  
}

