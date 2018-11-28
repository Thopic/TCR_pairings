// basic file operations
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <smmintrin.h>
using namespace std;


inline static int popcount(unsigned long long int x, unsigned long long int y)
{
  return __builtin_popcountll(x) + __builtin_popcountll(y);
}


int main (int argc, char ** argv) {

  if(argc < 5)
    return 0;

  std::string dfile = argv[1];
  std::string bf1 = argv[2];
  std::string bf2 = argv[3];
  std::string res = argv[4];
  
  std::string line;
  ifstream filedat (dfile);
  ifstream file1 (bf1);
  ifstream file2 (bf2);
  std::vector< std::vector < int > > d;
  unsigned int W = 0;
  int same = false;
  std::bitset<128> divider(0b00000000000000000000000000000000000000000000000000000000000000001111111111111111111111111111111111111111111111111111111111111111);
  
  if(filedat.is_open())
    {
      filedat >> W >> same;
      for(unsigned int i = 0; i <= W; ++i)
	d.push_back(std::vector< int >(W+1, 0));
      unsigned int n1, n2;
      int dv;
      while(filedat >> n1 >> n2 >> dv)
	{
	  d[n1][n2] = dv;
	  //	  std::cout << n1 << " " << n2 << " " << dv << std::endl;
	}
    }
  std::cout << "data read." << std::endl;
  std::cout << W << " " << (same ? "Same type" : "Different type") << std::endl;

  std::vector< std::vector<unsigned long long int > > pair1l(W+1, std::vector<unsigned long long int>()), pair2l(W+1, std::vector<unsigned long long int>()), pair1h(W+1, std::vector<unsigned long long int>()), pair2h(W+1, std::vector<unsigned long long int>());
  std::vector< std::vector< int > > pos1(W+1, std::vector< int>()), pos2(W+1, std::vector< int>());
  
  int i = 0;
  if(file1.is_open())
    {
      while(getline(file1, line))
	{
	  auto rKey = (std::bitset<128>(line) & divider).to_ullong();
	  auto lKey = ((std::bitset<128>(line) >> 64) & divider).to_ullong();
	  pair1h[popcount(rKey, lKey)].push_back(rKey);
	  pair1l[popcount(rKey, lKey)].push_back(lKey);
	  pos1[popcount(rKey, lKey)].push_back(i);
	  i += 1;
	}
    }
  file1.close();
  
  std::cout << "first file read." << std::endl;
  i = 0;
  if(file2.is_open())
    {
      while(getline(file2, line))
	{
	  auto rKey = (std::bitset<128>(line) & divider).to_ullong();
	  auto lKey = ((std::bitset<128>(line) >> 64) & divider).to_ullong();
	  pair2h[popcount(rKey, lKey)].push_back(rKey);
	  pair2l[popcount(rKey, lKey)].push_back(lKey);
	  pos2[popcount(rKey, lKey)].push_back(i);

	  i += 1;
	}
    }
  file2.close();
  std::cout << "second file read." << std::endl;

  // std::cout << "n: # of 1 with n bits set | # of 2 with n bits set" << std::endl;
  // for(unsigned int n = 0; n < W+1; ++n)
  //   {
  //     std::cout << n << ": " << pair1h[n].size() << " | " << pair2h[n].size() << std::endl;
  //   }
  
  ofstream fileres (res);
  std::cout << "Start computing the pairs" << std::endl;
  int dn1n2 = 0; 
  unsigned int startj = 0;
  for(unsigned int n1 = 0; n1 < W+1; ++n1)
    {
      std::cout << "\r" << std::setw(2) << n1 << std::flush;
      for(unsigned int n2 = 0; n2 < W+1; ++n2)
	{
	  std::vector<int> pairs_a, pairs_b;
	  dn1n2 = d[n1][n2];

	  auto p1h = pair1h[n1];
	  auto p1l = pair1l[n1];
	  auto p2h= pair2h[n2];
	  auto p2l = pair2l[n2];
	  
	  unsigned int s2 = pair2h[n2].size();
	  if(dn1n2 >= 0)
	    {
	      for(unsigned int i = 0; i < p1h.size(); ++i)
		{
		  startj = 0;
		  if(same and n1 == n2)
		    startj = i+1;
		  for(unsigned int j = startj; j < s2; ++j)
		    {
		      if(popcount(p1h[i]^p2h[j],p1l[i]^p2l[j]) <= dn1n2)
			{
			  pairs_a.push_back(i);
			  pairs_b.push_back(j);
			}
		    }
		}
	    
	      for(unsigned int i = 0; i < pairs_a.size(); ++i)
		fileres << pos1[n1][pairs_a[i]] << " " << pos2[n2][pairs_b[i]] << std::endl;
	      pairs_a.clear();
	      pairs_b.clear();
	    }
	}
    }
  std::cout << std::endl;
  return 0;
}


