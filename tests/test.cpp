#include <Eigen/CXX11/Tensor>


#include <string>
#include <iostream>
struct Person {
  std::string name;
  int age;
};

int main() {
  Eigen::Tensor<Person, 3> t(2, 3, 4);
  t(0,0,0) = {"yjp", 20};
  std::cout << t(0,0,0).name << "\n";
  std::cout << t(0,0,0).age << "\n";
  return 0;
}