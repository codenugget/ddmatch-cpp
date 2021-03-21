#include <cmath>
#include <iostream>

#include "core/MyFft.h"

int main(int argc, char** argv)
{
  //std::vector<std::complex<double>> orig{{1,-1},{4,4},{-3,3},{-2,-2}};
  std::vector<std::complex<double>> orig{{1,0},{4,0},{3,0},{2,0}};
  std::cout << "orig: ";
  for(auto e : orig)
    std::cout << e << ", ";

  std::cout << "\n\ndft: ";
  auto r1 = dft(orig);
  for(auto e : r1)
    std::cout << e << ", ";
  std::cout << "\nidft: ";
  auto inv1 = idft(r1);
  for(auto e : inv1)
    std::cout << e << ", ";

  auto r2 = FFT(orig);
  std::cout << "\n\nFFT: ";
  for(auto e : r2)
    std::cout << e << ", ";
  std::cout << "\nIFFT: ";
  auto inv2 = IFFT(r2);
  for(auto e : inv2)
    std::cout << e << ", ";



  int M = orig.size();
  std::complex<double> W = std::exp(std::complex(0.0, -2.0*M_PI/M));
  std::complex<double> A = std::complex(1.0, 0.0);
  auto r3 = CZT(orig, M, W, A);
  std::cout << "\n\nCZT: ";
  for(auto e : r3)
    std::cout << e << ", ";
  std::cout << "\nICZT: ";
  auto inv3 = ICZT(r3, M, W, A);
  for(auto e : inv3)
    std::cout << e << ", ";
  std::cout << "\n";


  exit(0);
}
