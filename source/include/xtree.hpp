#pragma once

#include <vector>
#include "rectangle.hpp"

template<typename T, size_t N, size_t M, size_t m = 40>
class Xtree {
  
 public:

  struct Node;
  
  struct Cell {
    Rectangle<N> MBR;
    std::shared_ptr<Node> child;
    std::shared_ptr<T> data; 
  };

  struct Node {
    typedef typename std::vector<Cell>::iterator iterator;

    iterator begin();
    iterator end();
    Cell& operator[](size_t index);
    
    Node(size_t mult);
    
    bool is_leaf();
    size_t max_size();
    size_t min_size();
    std::shared_ptr<Node> insert(const Cell &new_entry);
    
    std::vector<Cell> entry;
    size_t multiplier;
    size_t size;
  };
  
  Xtree();
  ~Xtree();
  
  size_t dimension() const;
  size_t size() const;
  bool empty() const;

  void insert(Rectangle<N> box, const T data);
  
  std::shared_ptr<Node> choose_subtree(const std::shared_ptr<Node> current_node,
                                       Rectangle<N> &box,
                                       const std::shared_ptr<T> data_ptr);
  
  std::shared_ptr<Node> choose_leaf_node(const std::shared_ptr<Node> current_node,
                                         Rectangle<N> &box,
                                         size_t &cell_pos);
  
  std::shared_ptr<Node> choose_dir_node(const std::shared_ptr<Node> current_node,
                                        Rectangle<N> &box,
                                        size_t &cell_pos);
  
  std::shared_ptr<Node> adjust_tree(const std::shared_ptr<Node> &parent,
                                    const std::shared_ptr<Node> &left,
                                    const std::shared_ptr<Node> &right,
                                    size_t cell_pos);

  //std::vector<T> KNNquery(Rectangle<N> &box, size_t k);

  //private:

  size_t tree_entries;
  std::shared_ptr<Node> root;
  
};

// Node !!

template<typename T, size_t N, size_t M, size_t m>
Xtree<T, N, M, m>::Node::Node(size_t mult) : multiplier(mult) {}


template<typename T, size_t N, size_t M, size_t m>
typename Xtree<T, N, M, m>::Node::iterator
Xtree<T, N, M, m>::Node::begin() {
  return entry.begin();
}

template<typename T, size_t N, size_t M, size_t m>
typename Xtree<T, N, M, m>::Node::iterator
Xtree<T, N, M, m>::Node::end() {
  return entry.end();
}

template<typename T, size_t N, size_t M, size_t m>
typename Xtree<T, N, M, m>::Cell&
Xtree<T, N, M, m>::Node::operator[](size_t index) {
  return entry[index];
}

template<typename T, size_t N, size_t M, size_t m>
bool Xtree<T, N, M, m>::Node::is_leaf() {
  if (size && entry[0].child) {
    return false;
  }
  return true;
}

template<typename T, size_t N, size_t M, size_t m>
size_t Xtree<T, N, M, m>::Node::max_size() {
  return M * multiplier;
}

template<typename T, size_t N, size_t M, size_t m>
size_t Xtree<T, N, M, m>::Node::min_size() {
  return max_size() * (m / 100.0);
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::Node::insert(const Cell &new_entry) {
  if (size < max_size()) {
    entry.push_back(new_entry);
    size++;
    return nullptr;
  }
  return nullptr;
}

// Xtree!!

template<typename T, size_t N, size_t M, size_t m>
Xtree<T, N, M, m>::Xtree() : root(std::make_shared<Node>(1)), tree_entries(0) {}

template<typename T, size_t N, size_t M, size_t m>
Xtree<T, N, M, m>::~Xtree() {
  root.reset();
}

template<typename T, size_t N, size_t M, size_t m>
size_t Xtree<T, N, M, m>::dimension() const{
  return M;
}

template<typename T, size_t N, size_t M, size_t m>
size_t Xtree<T, N, M, m>::size() const{
  return tree_entries;
}

template<typename T, size_t N, size_t M, size_t m>
bool Xtree<T, N, M, m>::empty() const{
  return !tree_entries;
}

template<typename T, size_t N, size_t M, size_t m>
void Xtree<T, N, M, m>::insert(Rectangle<N> box, const T data) {
  std::shared_ptr<T> data_ptr = std::make_shared<T>(data);
  std::shared_ptr<Node> splitted_node = choose_subtree(root, box, data_ptr);
  if (!splitted_node)
    return;
  // Split the root !
  std::shared_ptr<Node> new_root = std::make_shared<Node>(root->multiplier);
  (*new_root)[0].child = root;
  (new_root->size)++;
  adjust_tree(new_root, root, splitted_node, 0);
  root = new_root;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::choose_subtree(const std::shared_ptr<Node> current_node,
                                  Rectangle<N> &box,
                                  const std::shared_ptr<T> data_ptr) {
  // R* choose_subtree!!
  if (!current_node->is_leaf()) {
    size_t cell_pos;
    std::shared_ptr<Node> next_node;
    if((*current_node)[0].child->is_leaf())
      next_node = choose_leaf_node(current_node, box, cell_pos);
    else
      next_node = choose_dir_node(current_node, box, cell_pos);
    std::shared_ptr<Node> splitted_node = choose_subtree(next_node, box, data_ptr);
    return adjust_tree(current_node, next_node, splitted_node, cell_pos);
  }
  Cell new_entry;
  new_entry.MBR = box;
  new_entry.data = data_ptr;
  new_entry.child = nullptr;
  tree_entries++;
  return current_node->insert(new_entry);
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::choose_dir_node(const std::shared_ptr<Node> current_node,
                                   Rectangle<N> &box,
                                   size_t &cell_pos) {
  Rectangle<N> enlarged_box = (*current_node)[0].MBR;
  enlarged_box.adjust(box);
  float minimum_area = (*current_node)[0].MBR.get_area();
  float minimum_enlargement = enlarged_box.get_area() - minimum_area;
  std::shared_ptr<Node> node = (*current_node)[0].child;
  float enlargement, area;

  cell_pos = 0;
  for (size_t i = 1; i < current_node->size; ++i) {
    Cell &current_entry = (*current_node)[i];
    area = current_entry.MBR.get_area();
    enlarged_box = current_entry.MBR;
    enlarged_box.adjust(box);
    enlargement = enlarged_box.get_area() - area;

    if (enlargement < minimum_enlargement ||
        (enlargement == minimum_enlargement && area < minimum_area)) {
      minimum_enlargement = enlargement;
      minimum_area = area;
      node = current_entry.child;
      cell_pos = i;
    }
  }
  return node;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::choose_leaf_node(const std::shared_ptr<Node> current_node,
                                    Rectangle<N> &box,
                                    size_t &cell_pos) {
  Rectangle<N> enlarged_box = (*current_node)[0].MBR;
  enlarged_box.adjust(box);
  float minimum_enlargement = enlarged_box.get_area() - (*current_node)[0].MBR.get_area();
  float minimum_overlap = 0;
  for (size_t i = 1; i < current_node->size; ++i) {
    minimum_overlap += enlarged_box.get_overlap((*current_node)[i].MBR);
  }
  float enlargement, overlap;
  std::shared_ptr<Node> node = (*current_node)[0].child;

  cell_pos = 0;
  for (size_t i = 1; i < current_node->size; ++i) {
    Cell &current_entry = (*current_node)[i];
    enlarged_box = current_entry.MBR;
    enlarged_box.adjust(box);
    enlargement = enlarged_box.get_area() - current_entry.MBR.get_area();
    overlap = 0;
    for (size_t j = 0; j < current_node->size; ++j) {
      if (i == j) continue;
      overlap += enlarged_box.get_overlap((*current_node)[j].MBR);
    }
    if (overlap < minimum_overlap ||
        (overlap == minimum_overlap && enlargement < minimum_enlargement)) {
      minimum_overlap = overlap;
      minimum_enlargement = enlargement;
      node = current_entry.child;
      cell_pos = i;
    }
  }
  return node;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::adjust_tree(const std::shared_ptr<Node> &parent,
                               const std::shared_ptr<Node> &left,
                               const std::shared_ptr<Node> &right,
                               size_t cell_pos) {
  (*parent)[cell_pos].MBR.reset();
  for (Cell &current_entry : *left) {
    (*parent)[cell_pos].MBR.adjust(current_entry.MBR);
  }
  if (!right) {
    return nullptr;
  }
  Cell new_entry;
  new_entry.MBR.reset();
  for (Cell &current_entry : *right) {
    new_entry.MBR.adjust(current_entry.MBR);
  }
  new_entry.child = right;
  return parent->insert(new_entry);
}
