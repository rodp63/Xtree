// Copyright 2020 Roger Peralta Aranibar Advanced Data Estructures
#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <chrono>

#include "rectangle.hpp"
#include "xtree.hpp"

#define DIM 14
#define POINTS 170653
#define PATH "../spotify_dataset.csv"

typedef std::pair<int, std::string> data_type; // {year, name}

std::map<size_t, float> normalizer;

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
    for (std::pair<size_t,float> norm : normalizer) {
      pt[norm.first] /= norm.second;
    }
    // Get song data
    points >> year;
    getline(points, info);
    Cake.insert(pt, std::make_pair(year, info));
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
  // Normalize values
  normalizer[3] = 350000.0; // duration_ms
  normalizer[7] = 11.0; // key
  normalizer[9] = 60.0; // loudness
  normalizer[11] = 100.0; // popularity
  normalizer[13] = 245.0; // tempo
  
  std::cout<<"*---------------------------------------------*"<<std::endl;
  std::cout<<"*---------- X Tree by Joaquin Palma ----------*"<<std::endl;
  std::cout<<"*---------------------------------------------*"<<std::endl;
  std::cout<<"|> Inserting Points ... "<<std::endl;
  
  Timer<int()> timed_built(build_data_structure, "Index");
  timed_built();

  std::cout<<"\n*-------- K Nearest Neighbors Queries --------*\n"<<std::endl;
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
    for (std::pair<size_t,float> norm : normalizer) {
      query[norm.first] /= norm.second;
    }
    
    std::vector<std::vector<float>> result = timed_query(query, k);

    for (int i = 0; i < k; ++i) {
      // get normalized dist
      float dist = 0;
      for (size_t j = 0; j < DIM; ++j) {
        dist += (result[i][j] - query[j]) * (result[i][j] - query[j]);
      }
      dist = sqrt(dist);
      std::cout<<"   #"<<i+1<<" [";
      std::cout<<std::fixed<<std::setprecision(5)<<dist;
      std::cout<<"] -> { ";
      for (std::pair<size_t,float> norm : normalizer) {
        result[i][norm.first] *= norm.second;
      }
      for (int j = 0; j < DIM; ++j)
        std::cout<<std::fixed<<std::setprecision(4)<<result[i][j]<<" ";
      std::cout<<"}\n";
    }
    std::cout<<'\n';
  }
}
