/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "split.h"
#include "cdflib.hpp"
#include "pdflib.hpp"
#include "var.hpp"
#include "index.hpp"
#include "phase.hpp"

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>
#include <limits.h>
#include <memory>

using namespace std;
using namespace vcflib;
void printVersion(void){

	    exit(1);
}

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "     HapLRT is a likelihood ratio test for haplotype lengths.  The lengths are modeled with an exponential distribution.  " << endl;
  cerr << "     The sign denotes if the target has longer haplotypes (1) or the background (-1).                                     " << endl << endl;

  cerr << "Output : 4 columns :                             " << endl;
  cerr << "     1. seqid                                    " << endl;
  cerr << "     2. position                                 " << endl;
  cerr << "     3. mean target haplotype length             " << endl;
  cerr << "     4. mean background haplotype length         " << endl;
  cerr << "     5. p-value from LRT                         " << endl;
  cerr << "     6. sign                                     " << endl << endl;

  cerr << "INFO: Usage: hapLRT  --target 0,1,2,3,4,5,6,7 --background 11,12,13,16,17,19,22 --type GP --file my.vcf                                     " << endl;
  cerr << endl;

  cerr << "INFO: required: t,target     -- argument: a zero base comma separated list of target individuals corresponding to VCF columns        " << endl;
  cerr << "INFO: required: b,background -- argument: a zero base comma separated list of background individuals corresponding to VCF columns    " << endl;
  cerr << "INFO: required: f,file       -- argument: a properly formatted phased VCF file                                                       " << endl;
  cerr << "INFO: required: y,type       -- argument: type of genotype likelihood: PL, GL, GT or GP                                                  " << endl;
  cerr << "INFO: optional: r,region     -- argument: a genomic range to calculate hapLrt on in the format : \"seqid:start-end\" or \"seqid\" " << endl;
  cerr << endl;
  cerr << endl << "Type: genotype" << endl << endl;

  printVersion();

  exit(1);
}

void clearHaplotypes(std::vector<std::pair<std::string,std::string>>& haplotypes) {
    for (int i = 0; i < haplotypes.size(); i++) {
        haplotypes[i].first.clear();
        haplotypes[i].second.clear();
    }
}
//changed (carson) --> void findLengths(string haplotypes[][2], vector<int> group, int core, int lengths[]), int maxI){
void findLengths(std::vector<std::pair<std::string, std::string>>& haplotypes, vector<int> group, int core, int lengths[]){
  int gmax = group.size();
  int smax = haplotypes[0].first.length();
  int matrix[gmax*2][gmax*2];

  for(int i = 0 ; i < gmax*2; i++){
    lengths[i] = 0;
  }

  for(int i = 0; i < gmax*2; i++){ //get group member 1
    int g = (i < gmax) ? group[i] : group[i-gmax];

    string& currentHap = (i < gmax) ? haplotypes[g].first : haplotypes[g].second;
    for(int j = i+1; j < gmax*2; j++){ //get group member 2 for comparison
      int ag = (j < gmax) ? group[j] : group[j-gmax];

      string& altHap = (j < gmax) ? haplotypes[ag].first : haplotypes[ag].second;

      typedef string::const_iterator iter;
      iter cBeg = currentHap.begin(), cEnd = currentHap.end();
      iter aBeg = altHap.begin(), aEnd = altHap.end();

      iter citB = cBeg + core;
      iter citE = cBeg + core;
      iter aitB = aBeg + core;
      iter aitE = aBeg + core;

      int len = 0;
      if(*citB != *aitB){ //block length is 0 on no match
	matrix[i][j] = len;
	continue;
      }
      len++; //length is at least 1 on match

      //block longer than 1
      while(len < smax){
	int to_add = 0;
	if(citB > cBeg){
	  citB--;
	  aitB--;
	  if(*citB != *aitB) break;
	  to_add++;
	}

	if(citE < cEnd-1){
	  citE++;
          aitE++;
          if(*citE != *aitE) break;
	  to_add++;
	}

	len += to_add;
      }

      //set length
      matrix[i][j] = len;
      if(lengths[i] < len) lengths[i] = len;
      if(lengths[j] < len) lengths[j] = len;
    }
  }

  //check for bugs
  //for(int i = 0 ; i < gmax*2; i++){
  //  if(lengths[i] == -1){
  //  cerr << endl;
  //  cerr << core;
  //  for(int i = 0; i < gmax*2; i++){ //get group member 1
  //    int g = (i < gmax) ? group[i] : group[i-gmax];
  //    int c = (i < gmax) ? 0 : 1;
  //
  //    string& currentHap = haplotypes[g][c];
  //    typedef string::const_iterator iter;
  //    iter cBeg = currentHap.begin(), cEnd = currentHap.end();
  //    iter citB = cBeg + core;
  //    iter citE = cBeg + core;
  //
  //    cerr << *citB << endl;
  //  }
  //  cerr << endl;
  //  }
  //}
}

double mean(int data[], int n){

  if(!n)
    return -std::numeric_limits<double>::quiet_NaN();

  int sum = 0;

  for(int i = 0; i < n; i++){
    sum += data[i];
    //cerr << data[i] << endl;
  }
  return (double(sum) / double(n));
}

double lfactorial(int n) {

  double fact=0;
  double i;

  for (i=1.0; i<=n; i++){
    fact += log(i);
  }

  return fact;
}

// using the wikipedia's def

double lnbinomial(double k, double r, double m ){

  //parameterized by scale and mu
  // R example dnbinom(x=6, size=0.195089, mu=7.375, log=TRUE)

  double ans = lgamma( r+k ) - ( lfactorial(k)  + lgamma(r)  ) ;
  ans += log(pow((m/(r+m)),k))   ;
  ans += log(pow((r/(r+m)),r)) ;

  //cerr << "k: " << k << "\t" << "m: " << m << "\t" << "r: " << r << "\t" << "ans: " << ans << endl;

  return ans;

}

double lexp(double x, double lambda){

  double ans = lambda * pow ( exp(1) , (-lambda * x));

  return log(ans);

}

double totalLL(int dat[], int n, double m){

  double ll = 0;

  for(int j = 0; j < n; j++){
    ll += lexp(dat[j], 1/m);
  }

  return ll;
}

double var(int dat[], int n, double mean){

  double sum = 0;

  for(int i = 0; i < n; i++){
    sum += pow( (double(dat[i]) - mean), 2);
  }

  double var = sum / (n - 1);

  return var;

}

void calc(std::vector<std::pair<std::string, std::string>>& haplotypes, int nhaps, const vector<long int>& pos, const vector<double>& afs, const vector<int> & target, const vector<int> & background,  const vector<int>& total,  const string& seqid){

  //moved (carson)
  int tl = 2*target.size();
  int bl = 2*background.size();
  int al = 2*total.size();

  //treat as buffers to seed findLengths (carson)


  for(int snp = 0; snp < haplotypes[0].first.length(); snp++){

    int targetLengths[tl];
    int backgroundLengths[bl];
    int totalLengths[al];

    //changed (carson) --> findLengths(haplotypes, target,     snp, targetLengths, tl);
    //changed (carson) --> findLengths(haplotypes, background, snp, backgroundLengths, bl);
    findLengths(haplotypes, target,     snp, targetLengths);
    findLengths(haplotypes, background, snp, backgroundLengths);

    copy(targetLengths, targetLengths + tl, totalLengths);
    copy(backgroundLengths, backgroundLengths +bl, totalLengths + tl);


    double tm = mean(targetLengths, tl);
    double bm = mean(backgroundLengths, bl);
    double am = mean(totalLengths, al);

    double dir = 1;

    if(tm < bm){
      dir = -1;
    }


    double Alt = totalLL(targetLengths, 2*target.size(), tm)
      + totalLL(backgroundLengths, 2*background.size(), bm);


    double Null = totalLL(targetLengths, 2*target.size(), am)
      + totalLL(backgroundLengths, 2*background.size(), am);

    double l = 2 * (Alt - Null);

    if(l < 0){
      continue;
    }

    int     which = 1;
    double  p ;
    double  q ;
    double  x  = l;
    double  df = 2;
    int     status;
    double  bound ;

    cdfchi(&which, &p, &q, &x, &df, &status, &bound );

    cout << seqid << "\t" << pos[snp] << "\t" << tm << "\t" << bm  <<  "\t" << 1-p <<  "\t" << dir <<  endl;

  }
}

int main(int argc, char** argv) {

  // set the random seed for MCMC

  srand((unsigned)time(NULL));

  // the filename

  string filename = "NA";

  // set region to scaffold

  string region = "NA";

  // using vcflib; thanks to Erik Garrison

  VariantCallFile variantFile;

  // zero based index for the target and background individuals

  map<int, int> targetIndex, backgroundIndex;

  // deltaaf is the difference of allele frequency we bother to look at

  // ancestral state is set to zero by default

  // phased

  int phased = 0;

  string type = "NA";

    const struct option longopts[] =
    {
      {"version"   , 0, 0, 'v'},
      {"help"      , 0, 0, 'h'},
      {"file"      , 1, 0, 'f'},
      {"target"    , 1, 0, 't'},
      {"background", 1, 0, 'b'},
      {"region"    , 1, 0, 'r'},
      {"type"      , 1, 0, 'y'},

      {0,0,0,0}
    };

    int findex;
    int iarg=0;

    while(iarg != -1)
    {
      iarg = getopt_long(argc, argv, "y:r:t:b:f:hv", longopts, &findex);

      switch (iarg)
      {
        case 'h':
          printHelp();
        case 'v':
          printVersion();
        case 'y':
          type = optarg;
          break;
        case 't':
          loadIndices(targetIndex, optarg);
          cerr << "INFO: there are " << targetIndex.size() << " individuals in the target" << endl;
          cerr << "INFO: target ids: " << optarg << endl;
          break;
        case 'b':
          loadIndices(backgroundIndex, optarg);
          cerr << "INFO: there are " << backgroundIndex.size() << " individuals in the background" << endl;
          cerr << "INFO: background ids: " << optarg << endl;
          break;
        case 'f':
          cerr << "INFO: file: " << optarg  <<  endl;
          filename = optarg;
          break;
        case 'r':
          cerr << "INFO: set seqid region to : " << optarg << endl;
          region = optarg;
          break;
        default:
          break;
      }
    }

    map<string, int> okayGenotypeLikelihoods;
    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;
    okayGenotypeLikelihoods["GT"] = 1;

    if(type == "NA"){
      cerr << "FATAL: failed to specify genotype likelihood format : PL, GL, GT or GP" << endl;
      printHelp();
      return 1;
    }
    if(okayGenotypeLikelihoods.find(type) == okayGenotypeLikelihoods.end()){
      cerr << "FATAL: genotype likelihood is incorrectly formatted, only use: PL, GL GT or GP" << endl;
      printHelp();
      return 1;
    }

    if(filename == "NA"){
      cerr << "FATAL: did not specify a file" << endl;
      printHelp();
      return(1);
    }


    variantFile.open(filename);

    if(region != "NA"){
      if(!variantFile.setRegion(region)){ //check if region is even specified in header
        bool region_exists = false;
        vector<string> headerLines = split (variantFile.header, "\n");
        for(const auto& headerLine : headerLines){
          if(headerLine.substr(0,8) == "##contig"){
            string contigInfo = headerLine.substr(10, headerLine.length() -11);
            vector<string> info = split(contigInfo, ",");
            for(const auto& sub : info){
              vector<string> subfield = split(sub, "=");
              if (subfield[0] != "ID") break;
              if(subfield[1] == region){
                region_exists = true;
                break;
              }

              vector<string> subregion = split(region, ":");
              if(subfield[1] == subregion[0]){
                region_exists = true;
                break;
              }
            }
            if(region_exists){
              break;
            }
          }
        }

        if(region_exists){
          cerr << "WARNING: There are no variants for the specified region" << endl;
          return 0;
        }
        else{
          cerr << "FATAL: You specified an invalid region" << endl;
          return 1;
        }
      }
    }

    if (!variantFile.is_open()) {
        return 1;
    }


    Variant var(variantFile);

    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    vector<int> ibi, iti, itot;

    int index = 0, indexi = 0;

    // TODO: fix loop
    for(const auto& _ : samples){
      if(targetIndex.find(index) != targetIndex.end()){
        iti.push_back(indexi);
        //	itot.push_back(indexi);
        indexi++;
      }
      if(backgroundIndex.find(index) != backgroundIndex.end()){
        ibi.push_back(indexi);
        //	itot.push_back(indexi);
        indexi++;
      }
      index++;
    }

    //   itot.insert(itot.end(), iti.begin(), iti.end());

    itot = iti;
    itot.insert(itot.end(), ibi.begin(), ibi.end());

    vector<long int> positions;
    vector<double>   afs;

    std::vector<std::pair<std::string, std::string>> haplotypes(nsamples);

    string currentSeqid = "NA";

    int count = 0;
    while (variantFile.getNextVariant(var)) {
      count++;

      if(!var.isPhased()){
        cerr <<"FATAL: Found an unphased variant. All genotypes must be phased!" << endl;
        printHelp();
        return(1);
      }

      if (var.alt.size() > 1){
        continue;
      }

      if(currentSeqid != var.sequenceName){
        if(haplotypes[0].first.length() > 10){
          calc(haplotypes, nsamples, positions, afs, iti, ibi, itot, currentSeqid);
        }
        clearHaplotypes(haplotypes);
        positions.clear();
        currentSeqid = var.sequenceName;
        afs.clear();
      }

      vector < map< string, vector<string> > > target, background, total;

      int sindex = 0;

      for(int nsamp = 0; nsamp < nsamples; nsamp++){
        map<string, vector<string> > sample = var.samples[samples[nsamp]];

        if(!sample[type].size())
        {
          cerr << "Bad file format: genotype field " << type << " is not present for: " << var.sequenceName << " " << var.position << endl;
          exit(1);
        }

        if((type == "GL" || type == "GP" || type == "PL") && sample[type].size() != 3)
        {
          cerr << "Bad file format: genotype field " << type << " should have 3 values but has only ";
          cerr << sample[type].size() << " for: " << var.sequenceName << " " << var.position;
          cerr << " in sample " << nsamp << endl;
          exit(1);
        }

        if(targetIndex.find(sindex) != targetIndex.end() ){
          target.push_back(sample);
          total.push_back(sample);
        }
        if(backgroundIndex.find(sindex) != backgroundIndex.end()){
          background.push_back(sample);
          total.push_back(sample);
        }

        sindex += 1;
      }


      std::unique_ptr<genotype> populationTarget;
      std::unique_ptr<genotype> populationBackground;
      std::unique_ptr<genotype> populationTotal;

      if (type == "PL"){
        populationTarget     = std::make_unique<pl>();
        populationBackground = std::make_unique<pl>();
        populationTotal      = std::make_unique<pl>();
      } else if (type == "GL"){
        populationTarget     = std::make_unique<gl>();
        populationBackground = std::make_unique<gl>();
        populationTotal      = std::make_unique<gl>();
      } else if (type == "GP"){
        populationTarget     = std::make_unique<gp>();
        populationBackground = std::make_unique<gp>();
        populationTotal      = std::make_unique<gp>();
      } else if (type == "GT"){
        populationTarget     = std::make_unique<gt>();
        populationBackground = std::make_unique<gt>();
        populationTotal      = std::make_unique<gt>();
      }

      populationTarget->loadPop(target,         var.position);

      populationBackground->loadPop(background, var.position);

      populationTotal->loadPop(total,           var.position);


      if(populationTotal->af > 0.95 || populationTotal->af < 0.05){
        continue;
      }

      afs.push_back(populationTotal->af);
      positions.push_back(var.position);
      loadPhased(haplotypes, populationTotal.get());
    }  

    calc(haplotypes, nsamples, positions, afs, iti, ibi, itot, currentSeqid);

    return 0;
}
