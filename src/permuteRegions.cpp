/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

/*

This program was created at:  Fri Apr 17 14:59:53 2015
This program was created by:  Zev N. Kronenberg


Contact: zev.kronenber@gmail.com

Organization: Unviersity of Utah
    School of Medicine
    Salt Lake City, Utah


The MIT License (MIT)

Copyright (c) <2015> <Zev N. Kronenberg>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/
#include "split.h"
#include "stats.hpp"
#include "gpatInfo.hpp"

#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <omp.h>
// print lock
omp_lock_t lock;

struct options{
  std::string file    ;
  std::string smoothed;
  std::string format  ;
  double npermutation ;
  double nsuc         ;
  int threads         ;
  int chrIndex        ;
  int nIndex          ; 
  int valueIndex      ; 
}globalOpts;

struct score{
  std::string seqid;
  long int pos     ;
  double score     ;
};

struct smoothedInputData{
  std::string line;
  double score ; 
  double n     ;
  double nPer  ;
  double nSuc  ;
  double ePv   ;
};

static const char *optString = "x:y:u:f:n:s:";

using namespace std;

std::map<std::string, int> FORMATS;

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : NA

 Function does   : prints help

 Function returns: NA

*/
void printHelp()
{

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "     permuteSmoothFst is a method for adding empirical p-values  smoothed wcFst scores." << endl ;
  cerr << endl;
  cerr << "INFO: usage:  permuteSmoothFst -s wcFst.smooth.txt -f wcFst.txt -n 5 -s 1 "<< endl;
  cerr << endl;
  cerr << "Required:" << endl;
  cerr << "      file:     f   -- argument: original wcFst data     "<< endl;
  cerr << "      smoothed: s   -- argument: smoothed wcFst data     "<< endl;
  cerr << "      format:   y   -- argument: [swcFst, segwcFst]      "<< endl;
  cerr << "Optional:" << endl;
  cerr << "      number:   n   -- argument: the number of permutations to run for each value [1000]" << endl;
  cerr << "      success:  u   -- argument: stop permutations after \'s\' successes [1]"             << endl;
  cerr << "      success:  x   -- argument: number of threads [1]"             << endl;
  cerr << endl;
  cerr << "OUTPUT: permuteSmoothFst will append three additional columns:" << endl;
  cerr << "        1. The number of successes                            " << endl;
  cerr << "        2. The number of trials                               " << endl;
  cerr << "        3. The empirical p-value                              " << endl;
  cerr << endl;
  printVersion();

}


//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
{
  int opt = 0;
  globalOpts.file = "NA";
  
  globalOpts.nsuc         = 1;
  globalOpts.npermutation = 1000;
  
  opt = getopt(argc, argv, optString);
  while(opt != -1){
    switch(opt){
    case 'x':
      {
	globalOpts.threads = atoi(optarg);
      }
    case 'f':
      {
	globalOpts.file =  optarg;
	break;
      }
    case 'y':
      {
	globalOpts.format = optarg;
	if(FORMATS.find(globalOpts.format) == FORMATS.end()){
	  std::cerr << "FATAL: format not supported: " << globalOpts.format;
	  std::cerr << endl;
	  printHelp();
	  exit(1);
	}
	if(globalOpts.format == "swcFst"){
	  globalOpts.chrIndex   = 0;
	  globalOpts.nIndex     = 3;
	  globalOpts.valueIndex = 4;
	}
	if(globalOpts.format == "segwcFst"){
	  globalOpts.chrIndex   = 0;
	  globalOpts.valueIndex = 3;
	  globalOpts.nIndex     = 5;
	}	
	break;
      }
    case 'n':
      {
	globalOpts.npermutation = atof(((string)optarg).c_str());
	cerr << "INFO: permuteGPAT++ will do N permutations: " << globalOpts.npermutation << endl;
	break;
      }
    case 's':
      {
	globalOpts.smoothed = optarg;
	cerr << "INFO: smoothed file: " << globalOpts.smoothed << endl;
	break;
      }	
    case 'u':
      {
	globalOpts.nsuc = atof(((string)optarg).c_str());
	cerr << "INFO: permuteGPAT++ will stop permutations after N successes: " << globalOpts.nsuc << endl;
	break;
	}
    case '?':
      {
	break;
      }
    }
    
    opt = getopt( argc, argv, optString ); 
  }
  return 1;
}


//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : data, vector to load

 Function does   : returns a contguous window

 Function returns: bool

*/

bool getContiguousWindow(const vector<score> & data, vector<double> & load, int n){
  int r = rand() % data.size();    

  if(r+n >= data.size()){
    return false;
  }

  if(data[r].seqid != data[r+n].seqid){
    return false;
  }
  
  for(int i = r; i < r+n; i++){
    load.push_back(data[i].score);
  }
  return true;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : score, n, data

 Function does   : permutes the score.  requires that the window is contiguous

 Function returns: NA

*/

void permute(double s, int n, vector<score *> & data, 
	     double * nRep, double *nSuc, double * ePv){
 
  
  *ePv = 1 / globalOpts.npermutation;
  
  while(*nSuc < globalOpts.nsuc && *nRep < globalOpts.npermutation ){
    *nRep += 1;    


    std::vector<double> scores;

    bool getWindow = false;
    while(!getWindow){
      getWindow = getContiguousWindow(data, scores, n);
    }

    double ns = mean(scores);
    
    if(ns > s){
      *nSuc += 1;
    }
  }
  if(*nSuc > 0){
    *ePv = *nSuc / *nRep;
  }
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
srand (time(NULL));

 FORMATS["swcFst"]   = 1;
 FORMATS["segwcFst"] = 1;

 globalOpts.threads = 1;

int parse = parseOpts(argc, argv);

 omp_set_num_threads(globalOpts.threads);


 if(globalOpts.file.compare("NA") == 0){
   cerr << "FATAL: no file was provided" << endl;
   printHelp();
   exit(1);
 }
 if(globalOpts.format.empty()){
   std::cerr << "FATAL: no format specified." << std::endl;
   exit(1);
 }

 vector< score > data;

 ifstream wcDat (globalOpts.file.c_str());

 string line;

 if(wcDat.is_open()){

   while(getline(wcDat, line)){
     vector<string> region = split(line, "\t");
     // will change for other output
     double value = atof(region[4].c_str());

     if(globalOpts.format == "swcFst" || globalOpts.format == "segwcFst"){
       
       if(region.size() != 5){
	 cerr << "FATAL: wrong number of columns in wcFst input" << endl;
	 exit(1);
       } 
       if(value < 0 ){
	 value = 0;
       }
     }

     data.emplace_back();
     auto& sp = data.back();

     sp.seqid = region[0]              ;
     sp.pos   = atoi(region[1].c_str());
     sp.score = value                  ;
   }
 }
 else{
   cerr << "FATAL: could not open file: " << globalOpts.file << endl;
   exit(1);
 }


 wcDat.close();
 line.clear();

 cerr << "INFO: raw values to permute against: " << data.size() << endl;

 ifstream smoothedFile (globalOpts.smoothed.c_str());

 vector<smoothedInputData> sData;

 if(smoothedFile.is_open()){
   while(getline(smoothedFile, line)){
   
     vector<string> region = split(line, "\t");
     sData.emplace_back();
     auto& sp = sData.back();

     sp.line = line;
     sp.score = atof(region[globalOpts.valueIndex].c_str());
     sp.n     = atoi(region[globalOpts.nIndex].c_str());
     sp.nPer = 0;
     sp.nSuc = 0;
     sp.ePv  = 0;
     
   }
 }
 smoothedFile.close();

 cerr << "INFO: Number of smoothed windows to permute : " << sData.size() << endl;

#pragma omp parallel for schedule(dynamic, 20)
 for(int i = 0; i < sData.size(); i++){
   permute(sData[i].score, sData[i].n, 
	   data, &sData[i].nPer, 
	   &sData[i].nSuc, &sData[i].ePv);
   
   omp_set_lock(&lock);
   cout << sData[i].line 
	<< "\t" << sData[i].nSuc
	<< "\t" << sData[i].nPer 
	<< "\t" << sData[i].ePv << endl;
   omp_unset_lock(&lock);
 }

 
 return 0;
}
