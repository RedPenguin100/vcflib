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
#include "gpatInfo.hpp"
#include "phase.hpp"

#include <memory>
#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <getopt.h>

// maaas speed

#if defined HAS_OPENMP
#include <omp.h>
// print lock
omp_lock_t lock;
#endif


struct opts{
  int         threads             ;
  int         pos                 ;
  std::string filename            ;
  std::string mapFile             ;
  std::string seqid               ;
  std::string geneticMapFile      ;
  std::string type                ;
  std::string region              ;
  std::map<int, double> geneticMap;
  double      af                  ;

}globalOpts;


using namespace std;
using namespace vcflib;

void printHelp(void){
  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      " << endl << endl;

  cerr << R"( meltEHH provides the data to plot extended haplotype homozygosity
(EHH) curves and produces the data to generate the following plot:
<img src="https://github.com/vcflib/vcflib/blob/master/examples/example-ehh.png?raw=true" alt="" width=400>

INFO: help
INFO: description:
     meltEHH provides the data to plot EHH curves.
Output : 4 columns :
     1. seqid
     2. position
     3. EHH
     4. ref or alt [0 == ref]
Usage:
      meltEHH --target 0,1,2,3,4,5,6,7 --pos 10 --file my.phased.vcf  \
           --region chr1:1-1000 > STDOUT 2> STDERR

Params:
       required: t,target   <STRING>  A zero base comma separated list of target
                                     individuals corresponding to VCF columns
       required: r,region   <STRING>  A tabix compliant genomic range
                                     format: "seqid:start-end" or "seqid"
       required: f,file     <STRING>  Proper formatted and phased VCF.
       required: y,type     <STRING>  Genotype likelihood format: GT,PL,GL,GP
       required: p,position <INT>     Variant position to melt.
       optional: a,af       <DOUBLE>  Alternative alleles with frequencies less
                                     than [0.05] are skipped.

)" << endl;

  cerr << endl << "Type: statistics" << endl << endl;
  cerr << endl;

  printVersion();

  exit(1);
}


bool gDist(int start, int end, double * gd){

  if(globalOpts.geneticMap.find(start) == globalOpts.geneticMap.end()){
    return false;
  }
  if(globalOpts.geneticMap.find(end) == globalOpts.geneticMap.end()){
    return false;
  }
  *gd = abs(globalOpts.geneticMap[start] - globalOpts.geneticMap[end]);
  return true;
}

void loadGeneticMap(int start, int end){

  if(globalOpts.geneticMapFile.empty()){
    std::cerr << "WARNING: No genetic map." << std::endl;
    std::cerr << "WARNING: A constant genetic distance is being used: 0.001." << std::endl;
    return;
  }

  ifstream featureFile (globalOpts.geneticMapFile.c_str());

  string line;

  int lastpos      = 0;
  double lastvalue = 0;

  if(featureFile.is_open()){

    while(getline(featureFile, line)){

      vector<string> region = split(line, "\t");

      if(region.front() != globalOpts.seqid){
	std::cerr << "WARNING: seqid MisMatch: " << region.front() << " " << globalOpts.seqid << std::endl;
	continue;
      }

      int   pos = atoi(region[3].c_str()) ;
      double cm = atof(region[2].c_str()) ;

      if(lastpos == 0 && start > pos){
	lastpos = pos;
	continue;
      }

      int diff     = abs(pos - lastpos);
      double vdiff = abs(lastvalue - cm );
      double chunk = vdiff/double(diff);

      double running = lastvalue;

      for(int i = lastpos; i < pos; i++){
	globalOpts.geneticMap[i] = running;
	running += chunk;
      }

      if(pos > end){
	break;
      }


      lastpos = pos;
      lastvalue = cm;
    }
  }

  featureFile.close();

  if(globalOpts.geneticMap.size() < 1){
    std::cerr << "FATAL: Problem loading genetic map" << std::endl;
    exit(1);
  }
}

void clearHaplotypes(std::vector<std::pair<std::string, std::string>>& haplotypes) {
    for (int i = 0; i < haplotypes.size(); i++) {
        haplotypes[i].first.clear();
        haplotypes[i].second.clear();
    }
}
void countHaps(int nhaps, map<string, int> & targetH,
	       const std::vector<std::pair<std::string, std::string>>& haplotypes, int start, int end){

  for(int i = 0; i < nhaps; i++){

    std::string h1 =  haplotypes[i].first.substr(start, (end - start)) ;
    std::string h2 =  haplotypes[i].second.substr(start, (end - start)) ;

    if(targetH.find(h1)  == targetH.end()){
      targetH[h1] = 1;
    }
    else{
      targetH[h1]++;
    }
    if(targetH.find(h2)  == targetH.end()){
      targetH[h2] = 1;
    }
    else{
      targetH[h2]++;
    }
  }
}

void computeNs(map<string, int> & targetH, int start,
	       int end, double * sumT, char ref, bool dir){

  for(const auto& th : targetH){

    if(th.second < 2){
      continue;
    }


    // end is extending ; check first base
    if(dir){
      if( th.first[0] == ref){

	//	std::cerr << "count dat: " << th->first << " " << th->second << " " << ref << " " << dir << endl;


	*sumT += r8_choose(th.second, 2);
      }
    }

    // start is extending ; check last base
    else{

      int last = th.first.size() -1;
      if( th.first[last] == ref ){
	//	std::cerr << "count dat:" << th->first << " " << th->second << " " << ref << " " << dir << endl;


      	*sumT += r8_choose(th.second, 2);
      }
    }
  }
}

bool calcEhh(const std::vector<std::pair<std::string, std::string>>& haplotypes, int start,
	     int end, char ref, int nhaps,
	     double * ehh, double  div, bool dir){

  double sum = 0 ;
  map<string , int> refH;

  countHaps(nhaps, refH, haplotypes, start, end);
  computeNs(refH, start, end, &sum, ref, dir   );

  double internalEHH = sum / (r8_choose(div, 2));

  if(internalEHH > 1){
    std::cerr << "FATAL: internal error." << std::endl;
    exit(1);
  }

  *ehh = internalEHH;

  return true;
}

int integrate(std::vector<std::pair<std::string, std::string>>& haplotypes,
	      vector<long int> & pos,
	      bool         direction,
	      int               maxl,
	      int                snp,
	      char               ref,
	      int              nhaps,
	      double *           iHH,
	      double           denom ){

  double ehh = 1;

  int start = snp;
  int end   = snp;

  // controls the substring madness
  if(!direction){
    start += 1;
    end += 1;
  }

  while(ehh > 0.01){
    if(direction){
      end += 1;
    }
    else{
      start -= 1;
    }
    if(start < 0){
      return 1;
    }
    if(end > maxl){
      return 1;
    }
    double ehhRT = 0;
    if(!calcEhh(haplotypes,
		start, end,
		ref, nhaps,
		&ehhRT, denom,
		direction)){
      return 1;
    }

    if(ehhRT <= 0.01){
      return 0;
    }

    double delta_gDist = 0.001;

    if(direction){
      gDist(pos[end-1], pos[end], &delta_gDist);
    }
    else{
      gDist(pos[start + 1], pos[start], &delta_gDist);
    }
    *iHH += ((ehh + ehhRT)/2)*delta_gDist;

    if(direction){
    std::cout << pos[end] << "\t" << ehh << "\t" << ref << "\t" << direction << std::endl;
    }
    else{
    std::cout << pos[start] << "\t" << ehh << "\t" << ref << "\t" << direction << std::endl;
    }


    ehh = ehhRT;

  }

  return 10;
}

void calc(std::vector<std::pair<std::string, std::string>>& haplotypes, int nhaps,
	  const vector<double> &, vector<long int> & pos,
	  const vector<int> &, const vector<int> &, const string&){

  int maxl = haplotypes[0].first.length();


  for(int snp = 0; snp < maxl; snp++){

    if(pos[snp] != globalOpts.pos ){
      continue;
    }


    double ihhR     = 0;
    double ihhA     = 0;

    map<string , int> refH;

    countHaps(nhaps, refH, haplotypes, snp, snp+1);


    double denomP1 = double(refH["0"]);
    double denomP2 = double(refH["1"]);

    int refFail = 0;
    int altFail = 0;

    std::cout << pos[snp] << "\t" << "1" << "\t" << "0" << "\t" << "0" << std::endl;



    refFail += integrate(haplotypes, pos, true,  maxl, snp, '0', nhaps, &ihhR, denomP1);

    refFail += integrate(haplotypes, pos, false, maxl, snp, '0', nhaps, &ihhR,  denomP1);

    altFail += integrate(haplotypes, pos, true, maxl, snp,  '1', nhaps, &ihhA, denomP2);

    altFail += integrate(haplotypes, pos, false, maxl, snp,  '1', nhaps, &ihhA, denomP2);

  }
}

int main(int argc, char** argv) {

  globalOpts.threads = 1   ;
  globalOpts.af      = 0.05;

  // zero based index for the target and background individuals

  map<int, int> it, ib;

    const struct option longopts[] =
      {
	{"version"   , 0, 0, 'v'},
	{"help"      , 0, 0, 'h'},
        {"file"      , 1, 0, 'f'},
	{"target"    , 1, 0, 't'},
	{"region"    , 1, 0, 'r'},
	{"gen"       , 1, 0, 'g'},
	{"type"      , 1, 0, 'y'},
	{"threads"   , 1, 0, 'x'},
	{"af"        , 1, 0, 'a'},
	{"pos"       , 1, 0, 'p'},
	{0,0,0,0}
      };

    int findex;
    int iarg=0;

    while(iarg != -1)
      {
	iarg = getopt_long(argc, argv, "a:x:g:y:r:d:t:b:f:p:hv", longopts, &findex);

	switch (iarg)
	  {
	  case 'p':
	    {
	      globalOpts.pos = atoi(optarg);
	      break;
	    }

	  case 'a':
	    {
	      globalOpts.af = atof(optarg);
	      break;
	    }
	  case 'x':
	    {
	      globalOpts.threads = atoi(optarg);
	      break;
	    }
	  case 'g':
	    {
	      globalOpts.geneticMapFile = optarg;
	      break;
	    }
	  case 'h':
	    {
	      printHelp();
	      break;
	    }
	  case 'v':
	    {
	      printVersion();
	      break;
	    }
	  case 'y':
	    {
	      globalOpts.type = optarg;
	      break;
	    }
	  case 't':
	    {
	      loadIndices(it, optarg);
	      cerr << "INFO: there are " << it.size() << " individuals in the target" << endl;
	      cerr << "INFO: target ids: " << optarg << endl;
	      break;
	    }
	  case 'f':
	    {
	      cerr << "INFO: file: " << optarg  <<  endl;
	      globalOpts.filename = optarg;
	      break;
	    }
	  case 'r':
	    {
	      cerr << "INFO: set seqid region to : " << optarg << endl;
	      globalOpts.region = optarg;
	      break;
	    default:
	      break;
	    }
	  }
      }

#if defined HAS_OPENMP
  omp_set_num_threads(globalOpts.threads);
#endif

    map<string, int> okayGenotypeLikelihoods;
    okayGenotypeLikelihoods["PL"] = 1;
    okayGenotypeLikelihoods["GL"] = 1;
    okayGenotypeLikelihoods["GP"] = 1;
    okayGenotypeLikelihoods["GT"] = 1;


    // add an option for dumping

//    for(std::map<int, double>::iterator gm = geneticMap.begin(); gm != geneticMap.end(); gm++){
//      cerr << "pos: " << gm->first << " cm: " << gm->second << endl;
//    }

    if(globalOpts.type.empty()){
      cerr << "FATAL: failed to specify genotype likelihood format : PL or GL" << endl;
      printHelp();
      exit(1);
    }
    if(okayGenotypeLikelihoods.find(globalOpts.type) == okayGenotypeLikelihoods.end()){
      cerr << "FATAL: genotype likelihood is incorrectly formatted, only use: PL or GL" << endl;
      printHelp();
      exit(1);
    }

    if(globalOpts.filename.empty()){
      cerr << "FATAL: did not specify a file" << endl;
      printHelp();
      exit(1);
    }

    if(it.size() < 2){
      cerr << "FATAL: target option is required -- or -- less than two individuals in target\n";
      printHelp();
      exit(1);
    }

    // using vcflib; thanksErik

    VariantCallFile variantFile;

    variantFile.open(globalOpts.filename);

    if(globalOpts.region.empty()){
      cerr << "FATAL: region required" << endl;
      exit(1);
    }
    if(! variantFile.setRegion(globalOpts.region)){
      cerr <<"FATAL: unable to set region" << endl;
      exit(1);
    }

    if (!variantFile.is_open()) {
      exit(1);
    }

    Variant var( variantFile );
    vector<int> target_h, background_h;

    int index   = 0;
    int indexi  = 0;


    vector<string> samples = variantFile.sampleNames;
    int nsamples = samples.size();

    // TODO: fix loop
    for(const auto& _ : samples){

      if(it.find(index) != it.end() ){
		target_h.push_back(indexi);
		indexi++;
      }
      index++;
    }


    vector<long int> positions;

    vector<double> afs;

    std::vector<std::pair<std::string, std::string>> haplotypes(target_h.size());


    while (variantFile.getNextVariant(var)) {

      globalOpts.seqid = var.sequenceName;

      if(!var.isPhased()){
	cerr << "FATAL: Found an unphased variant. All genotypes must be phased!" << endl;
	exit(1);
      }

      if(var.alleles.size() > 2){
	continue;
      }

      vector < map< string, vector<string> > > target, background, total;

      int sindex = 0;

      for(int nsamp = 0; nsamp < nsamples; nsamp++){

	map<string, vector<string> > sample = var.samples[ samples[nsamp]];

	if(it.find(sindex) != it.end() ){
	  target.push_back(sample);
	}
	sindex += 1;
      }

      std::unique_ptr<genotype> populationTarget;

      if(globalOpts.type == "PL"){
	populationTarget     = std::make_unique<pl>();
      }
      if(globalOpts.type == "GL"){
	populationTarget     = std::make_unique<gl>();
      }
      if(globalOpts.type == "GP"){
	populationTarget     = std::make_unique<gp>();
      }
      if(globalOpts.type == "GT"){
	populationTarget     = std::make_unique<gt>();
      }

      populationTarget->loadPop(target, var.position);

      if(populationTarget->af <= globalOpts.af
	 || populationTarget->nref < 2
	 || populationTarget->nalt < 2){
	;
	continue;
      }
      positions.push_back(var.position);
      afs.push_back(populationTarget->af);
      loadPhased(haplotypes, populationTarget.get());
    }

    if(!globalOpts.geneticMapFile.empty()){
      cerr << "INFO: loading genetics map" << endl;
      loadGeneticMap(positions.front(), positions.back());
      cerr << "INFO: finished loading genetics map" << endl;
    }

    calc(haplotypes, target_h.size(), afs, positions,
	 target_h, background_h, globalOpts.seqid);
    clearHaplotypes(haplotypes);

    exit(0);

}
