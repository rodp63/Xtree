#include "xtree.hpp"

using namespace std;

int main() {

  Rectangle<2> A;
  A[0].first = 0;
  A[0].second = 10;
  A[1].first = 0;
  A[1].second = 10;
  Rectangle<2> B;
  B[0].first = 8;
  B[1].first = 5;
  B[0].second = 15;
  B[1].second = 15;
  cout<<B.get_overlap(A)<<endl;
  
  Xtree<int, 2, 5> ttt;
  const int a = 11;
  int x = 1;
  ttt.insert(A, 69);
}
