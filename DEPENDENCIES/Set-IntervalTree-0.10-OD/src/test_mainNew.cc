#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>

#include "interval_tree.h"

int main(int argc, char *argv[])
{
  int low;
  int high;

  IntervalTree<int> intervalTree;
  int start = clock();
  intervalTree.LoadTree("test.txt");
  int end = clock();
  //std::cout << intervalTree.str();
  //std::cout << std::endl;
  std::cout << "it took " << end - start << "ticks, or " << ((float)end - start)/CLOCKS_PER_SEC << "seconds." << std::endl;  
  
  return 0;
}

