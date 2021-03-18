#include <iostream>

#include "core/MyFft.h"

int main(int argc, char** argv)
{
  std::vector<std::complex<double>> v{{1,-1},{4,4},{-3,3},{-2,-2}};
  //std::vector<std::complex<double>> v{{1,0},{4,0},{3,0},{2,0}};
  std::cout << "orig\n";
  for(auto e : v)
    std::cout << e << "\n";
  std::cout << "dft\n";
  auto r = dft(v);
  std::vector<double> rr;
  for(auto e : r) {
    rr.push_back(std::abs(e));
    std::cout << e << ":\n";
  }

  std::cout << "inverse\n";
  //auto r2 = idft(rr);
  auto r2 = idft(r);
  for(auto e : r2)
    std::cout << e << ":\n";
  exit(0);
}
