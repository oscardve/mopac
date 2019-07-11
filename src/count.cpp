/* MoPAC FASTQ internal function
version: 1.0
date: May 2018
description: internal function used by the read.FATSQ module
author:  Oscar Villarreal
affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
contact: oscardvillarreal AT gmail.com
*/

#include<Rcpp.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<algorithm>
#include<vector>
#include<map>
#include<set>
using namespace std;
using namespace Rcpp;
using map_of_vectors = map<string,vector<string>>;
using map_of_maps = map<string,map<int,int>>;

// [[Rcpp::plugins("cpp11")]]

// Obtain sgRNA annotation:
map_of_vectors get_annotation(const string library){
  ifstream lib{library};
  string line, spacer, gene;
  map_of_vectors anno;
  while(getline(lib,line)){ //loop over file
    stringstream ss(line);
    ss>>spacer>>gene;
  // while(!lib.eof()){ //loop over file
  //   lib>>spacer>>gene;
    anno[spacer].push_back(gene);
  }
  lib.close();
  return(anno);
}

// Obtain table of counts from all fastq files:
map_of_maps get_counts(map_of_vectors &anno, const vector<string> &fastq, const int line_start, const int line_interval, const int spacer_start, const int spacer_length, const bool rev, const bool comp, ofstream &total){
  ifstream reads;
  string line, spacer;
  map_of_maps counts;
  // Loop over fastq files:
  total<<"File"<<'\t'<<"Total_reads"<<endl;
  for(unsigned i{};i<fastq.size();i++){
    int N{};//number of reads
    reads.open(fastq[i]);
    Rcpp::Rcout<<"Processing file: "<<fastq[i]<<std::endl;
    for(int j{};j<line_start;j++)
      getline(reads,line); //skip starting lines
    while(!reads.eof()){ //loop over lines
      N++;
      spacer = line.substr(spacer_start,spacer_length); //trim
      /*try{spacer = line.substr(line.find("CAATTC")-8,spacer_length);} //trim
      catch(out_of_range){
      for(int j{};j<line_interval;j++) getline(reads,line); //skip lines in-between spacers
      continue;
      }*/
      if(rev==true)
        reverse(spacer.begin(),spacer.end()); //reverse
      if(comp==true){
        replace(spacer.begin(),spacer.end(),'A','X');replace(spacer.begin(),spacer.end(),'T','A');replace(spacer.begin(),spacer.end(),'X','T'); //complement A<->T
        replace(spacer.begin(),spacer.end(),'C','X');replace(spacer.begin(),spacer.end(),'G','C');replace(spacer.begin(),spacer.end(),'X','G'); //complement C<->G
      }
      if(anno.find(spacer) != anno.end())
        ++counts[spacer][i]; //if the spacer is in the annotation library, increase thefrequency of this spacer in this file
      for(int j{};j<line_interval;j++)
        getline(reads,line); //skip lines in-between spacers
    }
    reads.close();
    total<<fastq[i].substr(fastq[i].find_last_of("/")+1)<<'\t'<<N<<endl;
}
  total.close();
  // Add zeroes where no sequence was found (i.first=spacer; i.second=<file,count> pairs):
  for(const auto &i:counts)
    for(unsigned j{};j<fastq.size();j++)
      if(counts[i.first].find(j) == counts[i.first].end())
        counts[i.first][j] = 0;
  return counts;
}

// [[Rcpp::export]]
void readFASTQC(const std::string library, const std::string id, const int spacer_start, const int spacer_length, const bool rev, const bool comp, const StringVector fastqIN){
  // Parse user input:
  int line_start{2}, line_interval{4};
  string file;
  vector<string>fastq;
  for(int i=0;i<fastqIN.size();i++){
    file = fastqIN(i);
    fastq.push_back(file); //fastq files
  }
  // Obtain library annotations:
  map_of_vectors anno = get_annotation(library);
  // Obtain table of counts:
  ofstream total;//number of total reads per file
  total.open(id+"/fastq_reads.txt");
  map_of_maps counts = get_counts(anno, fastq, line_start, line_interval, spacer_start, spacer_length, rev, comp, total);
  // Save to output:
  ofstream out{};
  out.open(id+"/fastq_counts.txt");
  out<<"sgRNA"<<'\t'<<"Gene"; //start header
  for(unsigned i{};i<fastq.size();i++)
    out<<'\t'<<fastq[i].substr(fastq[i].find_last_of("/")+1); //file names without path
  out<<endl;
  for(const auto& i:counts){
    out<<i.first<<'\t'<<anno[i.first][0]<<'\t'; //sgRNA index, spacer, gene  out<<i.first; //sgRNA index
    for(auto const &j:i.second)
      out<<'\t'<<j.second; //counts
    out<<endl;
  }
  out.close();
}

// Obtain table of counts from all fastq files:
map_of_maps get_counts_shiny(map_of_vectors &anno, const vector<string> &fastq, const int line_start, const int line_interval, const int spacer_start, const int spacer_length, const bool rev, const bool comp, ofstream &total, Rcpp::Function shinyF){
  ifstream reads;
  string line, spacer;
  map_of_maps counts;
  // Loop over fastq files:
  total<<"File"<<'\t'<<"Total_reads"<<endl;
  for(unsigned i{};i<fastq.size();i++){
    int N{};//number of reads
    reads.open(fastq[i]);
    Rcpp::Rcout<<"Processing file: "<<fastq[i]<<std::endl;
    shinyF(1.0/fastq.size(),"Processing File:",fastq[i].substr(fastq[i].find_last_of("/")+1));
    for(int j{};j<line_start;j++)
      getline(reads,line); //skip starting lines
    while(!reads.eof()){ //loop over lines
      N++;
      spacer = line.substr(spacer_start,spacer_length); //trim
      /*try{spacer = line.substr(line.find("CAATTC")-8,spacer_length);} //trim
      catch(out_of_range){
      for(int j{};j<line_interval;j++) getline(reads,line); //skip lines in-between spacers
      continue;
      }*/
      if(rev==true)
        reverse(spacer.begin(),spacer.end()); //reverse
      if(comp==true){
        replace(spacer.begin(),spacer.end(),'A','X');replace(spacer.begin(),spacer.end(),'T','A');replace(spacer.begin(),spacer.end(),'X','T'); //complement A<->T
        replace(spacer.begin(),spacer.end(),'C','X');replace(spacer.begin(),spacer.end(),'G','C');replace(spacer.begin(),spacer.end(),'X','G'); //complement C<->G
      }
      if(anno.find(spacer) != anno.end())
        ++counts[spacer][i]; //if the spacer is in the annotation library, increase thefrequency of this spacer in this file
      for(int j{};j<line_interval;j++)
        getline(reads,line); //skip lines in-between spacers
    }
    reads.close();
    total<<fastq[i].substr(fastq[i].find_last_of("/")+1)<<'\t'<<N<<endl;
  }
  total.close();
  // Add zeroes where no sequence was found (i.first=spacer; i.second=<file,count> pairs):
  for(const auto &i:counts)
    for(unsigned j{};j<fastq.size();j++)
      if(counts[i.first].find(j) == counts[i.first].end())
        counts[i.first][j] = 0;
      return counts;
}

// [[Rcpp::export]]
void readFASTQC_shiny(const std::string library, const std::string id, const int spacer_start, const int spacer_length, const bool rev, const bool comp, const StringVector fastqIN, Rcpp::Function shinyF){
  // Parse user input:
  int line_start{2}, line_interval{4};
  string file;
  vector<string>fastq;
  for(int i=0;i<fastqIN.size();i++){
    file = fastqIN(i);
    fastq.push_back(file); //fastq files
  }
  // Obtain library annotations:
  map_of_vectors anno = get_annotation(library);
  // Obtain table of counts:
  ofstream total;//number of total reads per file
  total.open(id+"/fastq_reads.txt");
  map_of_maps counts = get_counts_shiny(anno, fastq, line_start, line_interval, spacer_start, spacer_length, rev, comp, total, shinyF);
  // Save to output:
  ofstream out{};
  out.open(id+"/fastq_counts.txt");
  out<<"sgRNA"<<'\t'<<"Gene"; //start header
  for(unsigned i{};i<fastq.size();i++)
    out<<'\t'<<fastq[i].substr(fastq[i].find_last_of("/")+1); //file names without path
  out<<endl;
  for(const auto& i:counts){
    out<<i.first<<'\t'<<anno[i.first][0]<<'\t'; //sgRNA index, spacer, gene  out<<i.first; //sgRNA index
    for(auto const &j:i.second)
      out<<'\t'<<j.second; //counts
    out<<endl;
  }
  out.close();
}
