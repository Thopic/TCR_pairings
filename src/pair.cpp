/**
    pair.cpp
    Purpose: Pairs sequences depending on the wells they were found in
    
    Take four arguments:
    - filename_p_values a file whose first line contain the number of well W.
    The following W lines stores a WxW matrix, d. d[n1, n2] contains the minimum 
    number of wells two sequences contained respectively in n1 and n2 wells need
    to share to be paired.
    - file_1 : the first list of sequences, N_1 lines of the form "i 00011100..01",
    where i is the index of the sequence, and the binary number represents the
    occupied wells.
    - file_2 : the second list of sequences.
    - output : the file were the results will be written

    @author Thomas Dupic
*/

#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <bitset>
#include <iomanip>





struct st_int128
{
  uint64_t hi,lo;
};


inline int popcount(st_int128 a)
{
  return __builtin_popcountll(a.hi) + __builtin_popcountll(a.lo);
}


inline int popcountxor(st_int128 a, st_int128 b)
{
  return __builtin_popcountll(a.hi^b.hi) + __builtin_popcountll(a.lo^b.lo);
}


st_int128 read_wells(std::string wells)
{
  st_int128 a;
  std::bitset<128> bset = std::bitset<128>(wells);
  a.hi = 0;
  a.lo = 0;
  uint64_t mask = 1;
  for(std::size_t i=0; i<64; ++i)
    {
      if(bset[i])
	a.hi |= mask;
      if(i < 63)
	a.hi <<= 1;
    }
  for(std::size_t i=64; i<128; ++i)
    {
      if(bset[i])
	a.lo |= mask;
      if(i < 127)
	a.lo <<= 1;
    }
  return a;
}


std::vector< std::pair<unsigned long, unsigned long> > apply_fdr
(const std::vector< std::pair<unsigned long, unsigned long> >& pairs_found,
 const std::vector<st_int128>& p1, const std::vector<st_int128>& p2,
 const std::vector<double>& pvalues, double fdr)
{
  int n1 = popcount(p1[pairs_found[0].first]);
  int n2 = popcount(p2[pairs_found[0].second]);
  int nbdiff;
  
  std::vector<double> pvalue_pairs;
  for(const auto& pair: pairs_found)
    {
      nbdiff = popcountxor(p1[pair.first], p2[pair.second]); 
      pvalue_pairs.push_back(pvalues[(n1 + n2 - nbdiff)/2]);
    }

  // sort
  std::sort(pvalue_pairs.begin(), pvalue_pairs.end());
  double plimit = 0.;
  int nb_pairs_null = p1.size()*p2.size();
  for(unsigned long ii=0; ii < pvalue_pairs.size(); ++ii)
    {
      if(ii > 0 and nb_pairs_null*pvalue_pairs[ii] > fdr*(ii+1))
	{
	  plimit = pvalue_pairs[ii-1];
	  break;
	}
    }

  // add the valid pvalues to a new vector
  std::vector< std::pair<unsigned long, unsigned long> > pairs_vetted
    = std::vector< std::pair<unsigned long, unsigned long> >();
  for(const auto& pair: pairs_found)
    {
      nbdiff = popcountxor(p1[pair.first], p2[pair.second]); 
      if(pvalues[(n1 + n2 - nbdiff)/2] < plimit)
	pairs_vetted.push_back(pair);
    }
  return pairs_vetted;
}


int main (int argc, char ** argv)
{

  //Open the files
  if(argc < 7)
    {
      std::cout << "This program needs 6 arguments: \n"
		<< "dvalues, pvalues, file_1, file_2, output, fdr\n";
      return 0;
    }

  std::string filename_dvalues = argv[1];
  std::string filename_pvalues = argv[2];
  std::string filename_1 = argv[3];
  std::string filename_2 = argv[4];
  std::string filename_output = argv[5];
  double fdr = std::stod(std::string(argv[6]));
  bool same_sequences = (filename_1 == filename_2); 

  std::ifstream file_dvalues(filename_dvalues);
  std::ifstream file_pvalues(filename_pvalues);
  std::ifstream file_1(filename_1);
  std::ifstream file_2(filename_2);

  if(!file_dvalues.is_open() or !file_pvalues.is_open()
     or !file_1.is_open() or !file_2.is_open())
    {
      std::cout << "Error: Cannot find the files.\n";
      return 0;
    }

  std::cout << "Read p-value file\n";
  unsigned int W; // number of wells
  file_dvalues >> W;
  
  std::vector<int> dvalues((W+1)*(W+1), -1);
  unsigned int n1, n2;
  int dv;
  while(file_dvalues >> n1 >> n2 >> dv)
      dvalues[n1 + (W+1)*n2] = dv;

  std::vector< std::vector < double > > pvalues;
  for(unsigned int i = 0; i < (W+1)*(W+1); ++i)
    pvalues.push_back(std::vector< double >(W+1, 1.));
  unsigned int n12;
  double pv;
  while(file_pvalues >> n1 >> n2 >> n12 >> pv)
    {
      pvalues[n1 + (W+1)*n2][n12] = pv;
    }

  

  // Read sequences file
  std::cout << "Read sequence files\n";
  std::vector< std::vector<st_int128> > sequences_1(W+1),
    sequences_2(W+1);
  std::vector< std::vector<long> > indexes_1(W+1), indexes_2(W+1);

  long number_sequence;
  std::string wells;
  file_1 >> wells >> wells; // skip header line
  while(file_1 >> number_sequence >> wells)
    {
      st_int128 int_well = read_wells(wells);
      sequences_1[popcount(int_well)].push_back(int_well);
      indexes_1[popcount(int_well)].push_back(number_sequence);
    }
  file_1.close();
  file_2 >> wells >> wells; // skip header line
  while(file_2 >> number_sequence >> wells)
    {
      st_int128 int_well = read_wells(wells);
      sequences_2[popcount(int_well)].push_back(int_well);
      indexes_2[popcount(int_well)].push_back(number_sequence);
    }
  file_2.close();

  // Compute the pairs
  std::ofstream file_output (filename_output);
  std::cout << "Start computing the pairs:\n";
  std::vector< std::pair<unsigned long, unsigned long> > pairs_found;
  std::vector< std::pair<unsigned long, unsigned long> > pairs_vetted;

  for(unsigned int n1=0; n1 < W+1; ++n1)
    {
      for(unsigned int n2=0; n2 < W+1; ++n2)
	{
	  std::cout << "\r" << std::setw(2) << n1 << "/" << W
		    << "\t" << std::setw(2) << n2 << "/" << W
		    << "                               " << std::flush;
	  int dn1n2 = dvalues[n1 + (W+1)*n2];
	  if(dn1n2 < 0)
	    continue;

	  auto p1 = sequences_1[n1];
	  auto p2 = sequences_2[n2];
	  for(ulong i=0 ; i < p1.size(); ++i)
	    {
	      for(ulong j=(same_sequences and n1==n2 ? i+1:0) ; j < p2.size() ; ++j)
		{
		  if(popcountxor(p1[i], p2[j]) <= dn1n2)
		    {
		      pairs_found.push_back(std::make_pair(i, j));
		    }
		}
	    }

	  std::cout << "\r" << std::setw(2) << n1 << "/" << W
		    << "\t" << std::setw(2) << n2 << "/" << W 
		    << "\t" << "vetting (" << pairs_found.size()
		    << " sequences)" << std::flush;
	  
	  pairs_vetted.clear();
	  if(pairs_found.size() > 0)
	    pairs_vetted = apply_fdr(pairs_found, p1, p2,
				     pvalues[n1 + (W+1)*n2], fdr);

	  for(auto p: pairs_vetted)
	      file_output << indexes_1[n1][p.first] <<  "\t" << indexes_2[n2][p.second] << "\n";
	  pairs_found.clear();
	}
    }
  std::cout << "\n";

  return 0;
}





