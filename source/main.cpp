// Copyright 2020 Roger Peralta Aranibar Advanced Data Estructures
#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <fstream>

#include "rectangle.hpp"
#include "xtree.hpp"

#define DIM 14
#define POINTS 170653
#define PATH "../spotify_dataset.csv"

const float max_2 = 350000.0; // duration_ms
const float max_8 = 60.0; // loudness
const float max_12 = 244.0; // tempo


typedef std::pair<int, std::string> data_type;

template <typename>
class Timer;

template <typename R, typename... T>
class Timer<R(T...)> {
 public:
  typedef R (*function_type)(T...);
  function_type function;

  explicit Timer(function_type function, std::string process_name = "")
      : function_(function), process_name_(process_name) {}

  R operator()(T... args) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    R result = function_(std::forward<T>(args)...);

    end = std::chrono::high_resolution_clock::now();
    int64_t duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
            .count();

    std::cout << std::setw(10) << process_name_ << std::setw(30)
              << "Duration: " + std::to_string(duration) + " ns\n";
    return result;
  }

 private:
  function_type function_;
  std::string process_name_;
};

Xtree<data_type, DIM, 10> Cake;

int build_data_structure() {
  std::ifstream points(PATH);
  float value;
  int year;
  std::string info;
  getline(points, info); // ignore 1st line

  for (size_t p = 0; p < POINTS; ++p) {
    Point<DIM> pt;
    // Get point data
    for (size_t i = 0; i < DIM; ++i) {
      points >> pt[i];
    }
    pt[2] /= max_2;
    pt[8] /= max_8;
    pt[12] /= max_12;
    // Get song data
    points >> year;
    getline(points, info);
    Cake.insert(pt, std::make_pair(year, info));
    //std::cout<<Cake.size()<<std::endl;
  }
  return 0;
}

std::vector<std::vector<float>> query_knn(std::vector<float> query, int k) {
  Point<DIM> query_point;
  for (size_t i = 0; i < DIM; ++i) {
    query_point[i] = query[i];
  }
  std::vector<std::pair<Point<DIM>, data_type> > knnData = Cake.KNNquery(query_point, k);
  for (size_t i = 0; i < knnData.size(); ++i) {
    std::cout<<"   #"<<i+1<<" -> "<<knnData[i].second.second<<std::endl;
  }
  std::vector<std::vector<float> > result;
  for (size_t i = 0; i < knnData.size(); ++i) {
    result.push_back(knnData[i].first.get_vector());
  }
  return result;
}

int main() {
  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"---------- X Tree by Joaquin Palma ----------"<<std::endl;
  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"|> Inserting Points ... "<<std::endl;
  
  Timer<int()> timed_built(build_data_structure, "Index");
  timed_built();
  
  std::cout<<"\n-------- K Nearest Neighbors Queries --------\n"<<std::endl;
  while (true) {
    Timer<std::vector<std::vector<float>>(std::vector<float>, int)> timed_query(
        query_knn, "Query KNN");
    std::vector<float> query(DIM);
    int k;
    std::cout<<"|> Enter k: "; std::cin>>k;
    std::cout<<"|> Enter point coordinates: ";
    for (size_t i = 0; i < DIM; ++i) {
      std::cin>>query[i];
    }
    query[2] /= max_2;
    query[8] /= max_8;
    query[12] /= max_12;
    
    std::vector<std::vector<float>> result = timed_query(query, k);
    
    // Print Points
    for (int i = 0; i < k; ++i) {
      result[i][2] *= max_2;
      result[i][8] *= max_8;
      result[i][12] *= max_12;
      std::cout<<"   # { ";
      for (int j = 0; j < DIM; ++j)
        std::cout<<result[i][j]<<", ";
      std::cout<<"}"<<std::endl;
    }
    std::cout<<'\n';
  }
}
