#pragma once

#include "rectangle.hpp"

#define MAX_OVERLAP 100

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
    
    size_t choose_split_axis(const Cell &new_entry,
                             std::vector<size_t> &axis_order);
    
    size_t choose_split_index(const Cell &new_entry,
                              const size_t axis,
                              std::vector<size_t> &axis_order);
    
    std::shared_ptr<Node> topological_split(const Cell &new_entry);
    std::shared_ptr<Node> overlap_min_split(const Cell &new_entry);
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

  void insert(const Rectangle<N> box, const T data);
  
  std::shared_ptr<Node> choose_subtree(const std::shared_ptr<Node> current_node,
                                       const Rectangle<N> &box,
                                       const std::shared_ptr<T> data_ptr);
  
  std::shared_ptr<Node> choose_leaf_node(const std::shared_ptr<Node> current_node,
                                         const Rectangle<N> &box,
                                         size_t &cell_pos);
  
  std::shared_ptr<Node> choose_dir_node(const std::shared_ptr<Node> current_node,
                                        const Rectangle<N> &box,
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

// Node implementation!!

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
size_t Xtree<T, N, M, m>::Node::choose_split_axis(const Cell &new_entry,
                                                  std::vector<size_t> &axis_order) {
  Rectangle<N> MBR_1, MBR_2;
  float minimum_margin = std::numeric_limits<float>::max();
  float margin;
  size_t M_val = max_size();
  size_t m_val = min_size();
  size_t chosen_axis;
  
  for (size_t axis = 0; axis < N; ++axis) {
    std::vector<std::pair<std::pair<float, float>, size_t> > order;
    for (size_t i = 0; i < size; ++i) {
      order.push_back({{entry[i].MBR[axis].first, entry[i].MBR[axis].second}, i});
    }
    order.push_back({{new_entry.MBR[axis].first, new_entry.MBR[axis].second}, size});
    sort(order.begin(), order.end());
    // Distributions
    margin = 0;
    for (size_t i = 0; i < (M_val - 2*m_val + 2); ++i) {
      MBR_1.reset();
      MBR_2.reset();
      for (size_t ff = 0; ff < m_val + i; ++ff) {
        if (order[ff].second == size)
          MBR_1.adjust(new_entry.MBR);
        else
          MBR_1.adjust(entry[order[ff].second].MBR);
      }
      for (size_t ss = m_val + i; ss < M_val + 1; ++ss) {
        if (order[ss].second == size)
          MBR_2.adjust(new_entry.MBR);
        else
          MBR_2.adjust(entry[order[ss].second].MBR);
      }
      margin += MBR_1.get_margin();
      margin += MBR_2.get_margin();
    }
    if (margin < minimum_margin) {
      minimum_margin = margin;
      chosen_axis = axis;
      for (size_t i = 0; i <= size; ++i) {
        axis_order[i] = order[i].second;
      }
    }
  }
  return chosen_axis;
}

template<typename T, size_t N, size_t M, size_t m>
size_t Xtree<T, N, M, m>::Node::choose_split_index(const Cell &new_entry,
                                                   const size_t axis,
                                                   std::vector<size_t> &axis_order) {
  Rectangle<N> MBR_1, MBR_2;
  float minimum_overlap = std::numeric_limits<float>::max();
  float minimum_area = std::numeric_limits<float>::max();
  float overlap, area;
  size_t M_val = max_size();
  size_t m_val = min_size();
  size_t index;
  
  for (size_t i = 0; i < (M_val - 2*m_val + 2); ++i) {
    MBR_1.reset();
    MBR_2.reset();
    for (size_t ff = 0; ff < m_val + i; ++ff) {
      if (axis_order[ff] == size)
        MBR_1.adjust(new_entry.MBR);
      else
        MBR_1.adjust(entry[axis_order[ff]].MBR);
    }
    for (size_t ss = m_val + i; ss < M_val + 1; ++ss) {
      if (axis_order[ss] == size)
        MBR_2.adjust(new_entry.MBR);
      else
        MBR_2.adjust(entry[axis_order[ss]].MBR);
    }
    area = MBR_1.get_area() + MBR_2.get_area();
    overlap = MBR_1.get_overlap(MBR_2);
    if (overlap < minimum_overlap ||
        (overlap == minimum_overlap && area < minimum_area)) {
      minimum_overlap = overlap;
      minimum_area = area;
      index = i;
    }
  }
  return minimum_overlap > MAX_OVERLAP ? -1 : index;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::Node::topological_split(const Cell &new_entry) {
  std::vector<size_t> axis_order(size + 1);
  size_t axis = choose_split_axis(new_entry, axis_order);
  size_t index = choose_split_index(new_entry, axis, axis_order);
  
  if (index != -1) {
    std::shared_ptr<Node> new_node = std::make_shared<Node>(multiplier);
    size_t M_val = max_size();
    size_t m_val = min_size();
    std::vector<Cell> tmp = entry;
    tmp.push_back(new_entry);
    entry.clear();
    for (size_t ff = 0; ff < m_val + index; ++ff)
      entry.push_back(tmp[axis_order[ff]]);
    for (size_t ss = m_val + index; ss < M_val + 1; ++ss)
      new_node->entry.push_back(tmp[axis_order[ss]]);
    size = entry.size();
    new_node->size = new_node->entry.size();
    return new_node;
  }
  return nullptr;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::Node::overlap_min_split(const Cell &new_entry) {
  return nullptr;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::Node::insert(const Cell &new_entry) {
  if (size < max_size()) {
    entry.push_back(new_entry);
    size++;
    return nullptr;
  }
  std::shared_ptr<Node> new_node;
  // Try R* split
  new_node = topological_split(new_entry);
  if(new_node)
    return new_node;
  // Try split-history
  new_node = overlap_min_split(new_entry);
  if(new_node)
    return new_node;
  //Create supernode
  multiplier++;
  entry.push_back(new_entry);
  size++;
  return nullptr;
}

// Xtree implementation!!

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
void Xtree<T, N, M, m>::insert(const Rectangle<N> box, const T data) {
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
                                  const Rectangle<N> &box,
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
                                   const Rectangle<N> &box,
                                   size_t &cell_pos) {
  Rectangle<N> enlarged_box;
  std::shared_ptr<Node> node = nullptr;
  float minimum_area = std::numeric_limits<float>::max();
  float minimum_enlargement = std::numeric_limits<float>::max();
  float enlargement, area;
  cell_pos = 0;
  
  for (size_t i = 0; i < current_node->size; ++i) {
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
                                    const Rectangle<N> &box,
                                    size_t &cell_pos) {
  Rectangle<N> enlarged_box;
  std::shared_ptr<Node> node = (*current_node)[0].child;
  float minimum_enlargement = std::numeric_limits<float>::max();
  float minimum_overlap = std::numeric_limits<float>::max();
  float enlargement, overlap;
  cell_pos = 0;
  
  for (size_t i = 0; i < current_node->size; ++i) {
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
