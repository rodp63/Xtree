#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

#include "rectangle.hpp"

#define MAX_OVERLAP 0.2

template<typename T, size_t N, size_t M, size_t m = 40>
struct Xtree {
  
  struct Node;
  
  struct SplitHistory {
    struct SHNode {
      size_t value;
      bool leaf;
      std::shared_ptr<SHNode> left, right;
      SHNode();
    };
    
    SplitHistory(size_t mult);
    void insert(size_t axis, size_t a, size_t b);
    size_t get_and_split();

    std::shared_ptr<SHNode> root;
    std::vector<std::shared_ptr<SHNode> > node_ptr;
    size_t node_capacity;
  };
  
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
    
    std::shared_ptr<Node> topological_split(const Cell &new_entry, size_t &split_axis);
    std::shared_ptr<Node> overlap_min_split(const Cell &new_entry, size_t &split_axis);
    std::shared_ptr<Node> insert(const Cell &new_entry, size_t &split_axis);
    
    std::vector<Cell> entry;
    std::shared_ptr<SplitHistory> split_tree;
    size_t multiplier;
    size_t size;
  };
  
  Xtree();
  ~Xtree();
  
  size_t dimension() const;
  size_t size() const;
  bool empty() const;
  void print();

  void insert(const Point<N> point, const T data);
  
  std::shared_ptr<Node> choose_subtree(const std::shared_ptr<Node> current_node,
                                       const Rectangle<N> &box,
                                       const std::shared_ptr<T> data_ptr,
                                       size_t &split_axis);
  
  std::shared_ptr<Node> choose_leaf_node(const std::shared_ptr<Node> current_node,
                                         const Rectangle<N> &box,
                                         size_t &cell_pos);
  
  std::shared_ptr<Node> choose_dir_node(const std::shared_ptr<Node> current_node,
                                        const Rectangle<N> &box,
                                        size_t &cell_pos);
  
  std::shared_ptr<Node> adjust_tree(const std::shared_ptr<Node> &parent,
                                    const std::shared_ptr<Node> &left,
                                    const std::shared_ptr<Node> &right,
                                    size_t cell_pos,
                                    size_t split_axis,
                                    size_t &new_split_axis);

  void KNNsearch(const std::shared_ptr<Node> current_node, const Point<N> &query_point);
  std::vector<std::pair<Point<N>, T> > KNNquery(const Point<N> &query_point, const size_t k);

  size_t tree_entries;
  std::priority_queue<std::pair<float, std::shared_ptr<Cell> > > knn;
  std::shared_ptr<Node> root;
  
};

// Split History implementation!!

template<typename T, size_t N, size_t M, size_t m>
Xtree<T, N, M, m>::SplitHistory::SHNode::SHNode() : left(nullptr), right(nullptr), leaf(1) {}

template<typename T, size_t N, size_t M, size_t m>
Xtree<T, N, M, m>::SplitHistory::SplitHistory(size_t mult) : root(std::make_shared<SHNode>()) {
  node_capacity = M * mult;
  root->value = 0;
  node_ptr.resize(node_capacity);
  node_ptr[0] = root;
}

template<typename T, size_t N, size_t M, size_t m>
void Xtree<T, N, M, m>::SplitHistory::insert(size_t axis, size_t a, size_t b) {
  std::shared_ptr<SHNode> parent = node_ptr[a];
  parent->left = std::make_shared<SHNode>();
  parent->right = std::make_shared<SHNode>();
  parent->leaf = false;
  parent->value = axis;
  parent->left->value = a;
  parent->right->value = b;
  node_ptr[a] = parent->left;
  node_ptr[b] = parent->right;
}

// Node implementation!!

template<typename T, size_t N, size_t M, size_t m>
Xtree<T, N, M, m>::Node::Node(size_t mult) : size(0), multiplier(mult) {}


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
    // order {bounds, id}
    std::vector<std::pair<std::pair<float, float>, size_t> > order;
    for (size_t i = 0; i < size; ++i) {
      order.push_back({{entry[i].MBR[axis].first, entry[i].MBR[axis].second}, i});
    }
    order.push_back({{new_entry.MBR[axis].first, new_entry.MBR[axis].second}, size});
    std::sort(order.begin(), order.end());
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
  float overlap, area, overlap_percentage = 0;
  size_t M_val = max_size();
  size_t m_val = min_size();
  size_t index;
  
  // Distributions
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
      if (area)
        overlap_percentage = (overlap) / (area - overlap);
      minimum_area = area;
      index = i;
    }
  }
  //std::cout<<" # Overlap split: "<<minimum_overlap<<std::endl;
  return overlap_percentage > MAX_OVERLAP ? -1 : index;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::Node::topological_split(const Cell &new_entry, size_t &split_axis) {
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
    split_axis = axis;
    return new_node;
  }
  return nullptr;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::Node::overlap_min_split(const Cell &new_entry, size_t &split_axis) {
  return nullptr;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::Node::insert(const Cell &new_entry, size_t &split_axis) {
  split_axis = -1;
  if (size < max_size()) {
    entry.push_back(new_entry);
    size++;
    return nullptr;
  }
  std::shared_ptr<Node> new_node;
  // Try R* split
  new_node = topological_split(new_entry, split_axis);
  if(new_node)
    return new_node;
  // Try split-history
  new_node = overlap_min_split(new_entry, split_axis);
  if(new_node)
    return new_node;
  //Create supernode
  std::cout<<"SuperNode!!!!\n";
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
void Xtree<T, N, M, m>::insert(const Point<N> point, const T data) {
  Rectangle<N> box = point.get_rect();
  std::shared_ptr<T> data_ptr = std::make_shared<T>(data);
  size_t axis, tmp;
  std::shared_ptr<Node> splitted_node = choose_subtree(root, box, data_ptr, axis);
  if (!splitted_node)
    return;
  // Split the root !
  std::shared_ptr<Node> new_root = std::make_shared<Node>(1);
  new_root->split_tree = std::make_shared<SplitHistory>(1);
  Cell new_entry;
  new_entry.child = root;
  new_root->insert(new_entry, tmp);
  adjust_tree(new_root, root, splitted_node, 0, axis, tmp);
  root = new_root;
}

template<typename T, size_t N, size_t M, size_t m>
std::shared_ptr<typename Xtree<T, N, M, m>::Node>
Xtree<T, N, M, m>::choose_subtree(const std::shared_ptr<Node> current_node,
                                  const Rectangle<N> &box,
                                  const std::shared_ptr<T> data_ptr,
                                  size_t &split_axis) {
  // R* choose_subtree!!
  if (!current_node->is_leaf()) {
    size_t cell_pos, axis;
    std::shared_ptr<Node> next_node;
    if((*current_node)[0].child->is_leaf())
      next_node = choose_leaf_node(current_node, box, cell_pos);
    else
      next_node = choose_dir_node(current_node, box, cell_pos);
    std::shared_ptr<Node> splitted_node = choose_subtree(next_node, box, data_ptr, axis);
    return adjust_tree(current_node, next_node, splitted_node, cell_pos, axis, split_axis);
  }
  Cell new_entry;
  new_entry.MBR = box;
  new_entry.data = data_ptr;
  new_entry.child = nullptr;
  tree_entries++;
  return current_node->insert(new_entry, split_axis);
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
                               size_t cell_pos,
                               size_t split_axis, size_t &new_split_axis) {
  new_split_axis = -1;
  (*parent)[cell_pos].MBR.reset();
  for (Cell &current_entry : *left) {
    (*parent)[cell_pos].MBR.adjust(current_entry.MBR);
  }
  if (!right) {
    return nullptr;
  }
  //parent->split_tree->insert(split_axis, cell_pos, parent->size);
  Cell new_entry;
  new_entry.MBR.reset();
  for (Cell &current_entry : *right) {
    new_entry.MBR.adjust(current_entry.MBR);
  }
  new_entry.child = right;
  return parent->insert(new_entry, new_split_axis);
}

template<typename T, size_t N, size_t M, size_t m>
void Xtree<T, N, M, m>::KNNsearch(const std::shared_ptr<Node> current_node,
                                  const Point<N> &query_point) {
  if (current_node->is_leaf()) {
    for (Cell &cur_entry : current_node->entry) {
      float dist = cur_entry.MBR.MINDIST(query_point);
      if (dist < knn.top().first) {
        knn.pop();
        knn.push(std::make_pair(dist, std::make_shared<Cell>(cur_entry)));
      }
    }
  }
  else {
    std::vector<std::pair<float, int> > order;
    for (size_t i = 0; i < current_node->size; ++i) {
      order.emplace_back((*current_node)[i].MBR.MINDIST(query_point), i);
    }
    std::sort(order.begin(), order.end());
    for (std::pair<float, int> &elem : order) {
      if ((*current_node)[elem.second].MBR.MINDIST(query_point) < knn.top().first) {
        KNNsearch((*current_node)[elem.second].child, query_point);
      }
    }
  }
}

template<typename T, size_t N, size_t M, size_t m>
std::vector<std::pair<Point<N>, T> > Xtree<T, N, M, m>::KNNquery(const Point<N> &query_point,
                                                                 const size_t k) {
  for (size_t i = 0; i < k; ++i) {
    knn.push(std::make_pair(std::numeric_limits<float>::max(), nullptr));
  }
  KNNsearch(root, query_point);
  std::vector<std::pair<Point<N>, T> > result;
  while(!knn.empty()) {
    std::shared_ptr<Cell> cur = knn.top().second;
    knn.pop();
    Point<N> p;
    for (size_t i = 0; i < N; ++i) {
      p[i] = cur->MBR[i].first;
    }
    result.push_back(std::make_pair(p, *(cur->data)));
  }
  std::reverse(result.begin(), result.end());
  return result;
}

template<typename T, size_t N, size_t M, size_t m>
void Xtree<T, N, M, m>::print() {
  std::cout<<"Printing the tree ...";
  int lastLvl = -1;
  std::queue<std::pair<std::shared_ptr<Node>, int> > Q;
  Q.push(make_pair(root, 0));
  while (!Q.empty()) {
    std::shared_ptr<Node> current = Q.front().first;
    int lvl = Q.front().second;
    Q.pop();
    if (lvl != lastLvl) {
      std::cout<<std::endl;
      lastLvl = lvl;
    }
    std::cout<<"{ | ";
    for (Cell &obj : *current) {
      std::cout<<(obj.data ? *(obj.data) : -1);
      std::cout<<":("<<obj.MBR[0].first<<","<<obj.MBR[0].second<<") - ";
      std::cout<<"("<<obj.MBR[1].first<<","<<obj.MBR[1].second<<") | ";
      if (!current->is_leaf()) {
        Q.push(make_pair(obj.child, lvl+1));
      }
    }
    std::cout<<"} ";    
  }
  std::cout<<"\n";
}
