#include <iostream>

#include "core/MyFft.h"

int main(int argc, char** argv)
{
  std::vector<std::complex<double>> v{{1,-1},{4,4},{-3,3},{-2,-2}};
  //std::vector<std::complex<double>> v{{1,0},{4,0},{3,0},{2,0}};
  std::cout << "orig: ";
  for(auto e : v)
    std::cout << e << ", ";

  std::cout << "\ndft: ";
  auto r1 = dft(v);
  for(auto e : r1)
    std::cout << e << ", ";
  std::cout << "\nidft: ";
  auto inv1 = idft(r1);
  for(auto e : inv1)
    std::cout << e << ", ";

  auto r2 = FFT(v);
  std::cout << "\nFFT: ";
  for(auto e : r2)
    std::cout << e << ", ";
  std::cout << "\nIFFT: ";
  auto inv2 = IFFT(r2);
  for(auto e : inv2)
    std::cout << e << ", ";
  std::cout << "\n";

  exit(0);
}
