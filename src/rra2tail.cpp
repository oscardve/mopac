/* MoPAC 2-tail RRA internal function
version: 1.0
date: May 2018
description: internal function used by the RRA.2tail module
author:  Oscar Villarreal
affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
contact: oscardvillarreal AT gmail.com
*/

#include<Rcpp.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<algorithm>
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

#include"rngs.h"
#include"math_api.h"
#include"rvgs.h"

// INTERNAL CPP FUNCTIONS --------------------------------------------------------------------------

// sgRNA structure:
struct item {
  item(const string& name, const string& value)
    : name{name},value{stod(value)}{}
  item(const string& name, const double& percentile) //input the percentile directly (for permutation)
    : name{name},value{percentile},percentile{percentile}{}
  string name; //name of the sgRNA
  double value{}; //score of the sgRNA
  mutable double percentile{}; //sgRNA percentile
};

// The sgRNAs are automatically sorted by score when inserted into a set:
struct item_order{
  bool operator()(const item& a, const item& b) const{return a.value-0.000000001<b.value+0.000000001;}
};

// Calculate lo-values for one gene from beta distribution when the sgrnas contain name, score and percentile:
double ComputeLoValue(set<item,item_order> &items, const double& alpha){
  vector<double>percentiles;
  for(const auto &j:items) //sgRNA iteration
    percentiles.push_back(j.percentile);
  sort(percentiles.begin(),percentiles.end());
  double loValue{1},l{};
  for(unsigned k{};k<percentiles.size();k++){
    if(percentiles[k]<=alpha || k==0){ //if percentile < alpha or if it's the minimum percentile
      double tmp = BetaNoncentralCdf(l+1,items.size()-l,0.0,percentiles[l],1E-10);
      if(tmp<loValue) loValue = tmp;
      l++;
    }
  }
  return(loValue);
}

// Calculate lo-values for one gene from beta distribution when the "sgrnas" consist of permuted percentiles without name:
double ComputeLoValue(vector<double> &items, const double& alpha){
  double loValue{1},k{};
  for(const auto &j:items){ //sgRNA iteration
    if(j<=alpha || k==0){ //if percentile < alpha or if it's the minimum percentile
      double tmp = BetaNoncentralCdf(k+1,items.size()-k,0.0,j,1E-10);
      if(tmp<loValue) loValue = tmp;
      k++;
    }
  }
  return(loValue);
}

// Calculate lo-values for all genes in a single direction:
map<string,double> getLoValuesC(set<item,item_order> &items, map<string,set<item,item_order>> &groups, const double& alpha){
  // The statistical significance, as expressed in a P-value, is calculated as the fraction of permutation values that are at least as extreme as the original statistic, which was derived from non-permuted data.
  // We performed a permutation test where the sgRNAs are randomly assigned to genes (the numbers of sgRNAs targeting each gene remain unchanged). By default, 100 × ng permutations are performed, where ng is the number of genes.
  // Get percentiles for each sgRNA in each gene:
  double R{};
  map<string,double> percentiles;
  map<string,double> loValues;
  for(const auto &i:items){
    percentiles[i.name]=(R+R+1)/(2*items.size());
    R++;
  }
  // Get loValues for each gene:
  for(const auto &i:groups){ //gene iteration
    for(const auto &j:groups[i.first]){ //sgRNA iteration
      j.percentile = percentiles[j.name]; //normalized percentile
    }
    loValues[i.first] = ComputeLoValue(groups[i.first], alpha);
  }
  return(loValues);
}

// Calculate lo-values for all genes in both directions:
map<string,double> getLoValues2C(set<item,item_order> &items, map<string,set<item,item_order>> &groups, const double& alpha){
  // The statistical significance, as expressed in a P-value, is calculated as the fraction of permutation values that are at least as extreme as the original statistic, which was derived from non-permuted data.
  // We performed a permutation test where the sgRNAs are randomly assigned to genes (the numbers of sgRNAs targeting each gene remain unchanged). By default, 100 × ng permutations are performed, where ng is the number of genes.
  // Get percentiles for each sgRNA in each gene:
  double R{};
  map<string,double> percentiles1,percentiles2;
  map<string,double> loValues1,loValues2,loValues;
  for(const auto &i:items){
    percentiles1[i.name]=(R+R+1)/(2*items.size());
    percentiles2[i.name]=1-percentiles1[i.name];
    R++;
  }
  // Get loValues for each gene in one direction:
  for(const auto &i:groups){ //gene iteration
    for(const auto &j:groups[i.first]) //sgRNA iteration
      j.percentile = percentiles1[j.name]; //normalized percentile
    loValues1[i.first] = ComputeLoValue(groups[i.first], alpha);
  }
  // Get loValues for each gene in the opposite direction:
  for(const auto &i:groups){ //gene iteration
    for(const auto &j:groups[i.first]) //sgRNA iteration
      j.percentile = percentiles2[j.name]; //normalized percentile
    loValues2[i.first] = ComputeLoValue(groups[i.first], alpha);
  }
  for(const auto &i:groups) //gene iteration
    loValues[i.first] = min(loValues1[i.first],loValues2[i.first]);
  return(loValues);
}

// Calculate p values for all genes (in this function, 'groups' is used only to get gene names and sizes):
map<string,double> getPValuesC(map<string,set<item,item_order>> &groups, map<string,double> &loValues, const double& alpha){
  // Set up randomization:
  PlantSeeds(123456);
  int rand_passnum{100};
  if(100*groups.size()<100000) //if there are less than 1000 genes, increase permutations
    rand_passnum = 100000/groups.size() + 1;
  int scanPass = rand_passnum + 1; // ~ 100
  // Map gene size to number of genes:
  map<unsigned,unsigned> itemNumMap;
  for(const auto &i:groups)
    ++itemNumMap[groups[i.first].size()];
  // Compute p values:
  map<string,double> pValues;
  vector<double> randLoValue(groups.size()*scanPass);
  //ofstream out1{"testing1.txt"};
  for(const auto &i:itemNumMap){ //for each gene size
    if(i.first<=100){ //avoid permuting genes with >100 sgRNAs
      vector<double>perm1(i.first),perm2(i.first);
      Rcpp::Rcout<<"Permuting genes with "<<i.first<<" sgrnas..."<<std::endl;
      int randLoValueNum = 0;
      for(int l=0;l<scanPass;l++){
        for(const auto &j:groups){ //gene iteration
          for(unsigned m=0;m<i.first;m++){ //sgRNA iteration
            perm1[m] = Uniform(0.0,1.0);
            perm2[m] = 1-perm1[m];
          }
          sort(perm1.begin(),perm1.end());//all percentiles for one gene sorted
          sort(perm2.begin(),perm2.end());//all percentiles for one gene sorted
          double loValue1 = ComputeLoValue(perm1, alpha);
          double loValue2 = ComputeLoValue(perm2, alpha);
          randLoValue[randLoValueNum] = min(loValue1,loValue2);
          randLoValueNum++;
        }
      }
    }
    sort(randLoValue.begin(),randLoValue.end());
    for(const auto &j:groups) //gene iteration
      if(groups[j.first].size()==i.first) {
        double index = lower_bound(randLoValue.begin(),randLoValue.end(),loValues[j.first])-randLoValue.begin();
        pValues[j.first] = (2*index+1.0)/(2*randLoValue.size());
        //out1<<j.first<<'\t'<<index<<'\t'<<loValues[j.first]<<'\t'<<randLoValue[0]<<'\t'<<randLoValue[1]<<'\t'<<randLoValue[2]<<'\t'<<randLoValue[3]<<'\t'<<pValues[j.first]<<endl;
      }
  }
  return(pValues);
}

// R FUNCTIONS --------------------------------------------------------------------------------------

// [[Rcpp::export]]
DataFrame getPValues(StringVector genes, StringVector sgrnas, NumericVector values, double alpha){
  // Set up items and groups:
  set<item,item_order> items; //ordered sgrnas
  map<string,set<item,item_order>> groups; //ordered sgrnas grouped into genes
  string gene,sgrna;
  double value{};
  for(int i=0;i<values.size();i++){
    gene = genes(i);
    sgrna = sgrnas(i);
    value = values(i);
    item item1 = item{sgrna,value};
    items.insert(item1); //order sgrnas in real-time
    groups[gene].insert(item1);
  }
  // Get lo values:
  map<string,double> loValues = getLoValues2C(items, groups, alpha);
  // Get p values:
  map<string,double> pValues = getPValuesC(groups, loValues, alpha);
  // Output data frame:
  StringVector names(pValues.size());
  NumericVector sizes(pValues.size()),result(pValues.size()),pvalues1(pValues.size());
  int iter = 0;
  for(const auto &i:pValues){
    names(iter) = i.first;
    sizes(iter) = groups[i.first].size();
    result(iter) = loValues[i.first];
    pvalues1(iter) = pValues[i.first];
    iter++;
  }
  return DataFrame::create(_["Gene"]=names,_["Size"]=sizes,_["rho"]=result,_["p"]=pvalues1,_["stringsAsFactors"]=false);
}

// ##################################################################################################
/* MoPAC 1-tail RRA internal C++ function ###########################################################
 version: 1.0
date: May 2018
description: translation to C++ of the RRA function used by MAGeCK-RRA
translated by: Oscar Villarreal
affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
contact: oscardvillarreal AT gmail.com
*/
// ##################################################################################################

// INTERNAL CPP FUNCTIONS --------------------------------------------------------------------------

// Calculate lo-values for one gene from beta distribution when the sgrnas contain name, score and percentile:
double ComputeLoValue_1tail (set<item,item_order> &items){
  double loValue{1},k{};
  for(const auto &j:items){ //sgRNA iteration
    if(j.percentile<=0.1 || k==0){ //if percentile < alpha or if it's the minimum percentile
      double tmp = BetaNoncentralCdf(k+1,items.size()-k,0.0,j.percentile,1E-10);
      if(tmp<loValue) loValue = tmp;
      k++;
    }
  }
  return(loValue);
}

// Calculate lo-values for one gene from beta distribution when the "sgrnas" consist of permuted percentiles without name:
double ComputeLoValue_1tail (vector<double> &items){
  double loValue{1},k{};
  for(const auto &j:items){ //sgRNA iteration
    if(j<=0.1 || k==0){ //if percentile < alpha or if it's the minimum percentile
      double tmp = BetaNoncentralCdf(k+1,items.size()-k,0.0,j,1E-10);
      if(tmp<loValue) loValue = tmp;
      k++;
    }
  }
  return(loValue);
}

// Calculate lo-values for all genes:
map<string,double> getLoValuesC_1tail (set<item,item_order> &items, map<string,set<item,item_order>> &groups){
  // The statistical significance, as expressed in a P-value, is calculated as the fraction of permutation values that are at least as extreme as the original statistic, which was derived from non-permuted data.
  // We performed a permutation test where the sgRNAs are randomly assigned to genes (the numbers of sgRNAs targeting each gene remain unchanged). By default, 100 × ng permutations are performed, where ng is the number of genes.
  // Get percentiles for each sgRNA in each gene:
  double R{};
  map<string,double> percentiles;
  map<string,double> loValues;
  for(const auto &i:items){
    percentiles[i.name]=(R+R+1)/(2*items.size());
    R++;
  }
  // Get loValues for each gene:
  for(const auto &i:groups){ //gene iteration
    for(const auto &j:groups[i.first]){ //sgRNA iteration
      j.percentile = percentiles[j.name]; //normalized percentile
    }
    loValues[i.first] = ComputeLoValue_1tail(groups[i.first]);
  }
  return(loValues);
}

// Calculate p values for all genes:
map<string,double> getPValuesC_1tail (map<string,set<item,item_order>> &groups, map<string,double> &loValues){
  // Set up randomization:
  PlantSeeds(123456);
  int rand_passnum{100};
  if(100*groups.size()<100000) //if there are less than 1000 genes, increase permutations
    rand_passnum = 100000/groups.size() + 1;
  int scanPass = rand_passnum + 1; // ~ 100
  // Map gene size to number of genes:
  map<unsigned,unsigned> itemNumMap;
  for(const auto &i:groups)
    ++itemNumMap[groups[i.first].size()];
  // Compute p values:
  map<string,double> pValues;
  vector<double> randLoValue(groups.size()*scanPass);
  for(const auto &i:itemNumMap){ //for each gene size
    if(i.first<=100){ //avoid permuting genes with >100 sgRNAs
      vector<double>perm(i.first);
      Rcpp::Rcout<<"Permuting genes with "<<i.first<<" sgrnas..."<<std::endl;
      int randLoValueNum = 0;
      for(int l=0;l<scanPass;l++){
        for(const auto &j:groups){ //gene iteration
          for(unsigned m=0;m<i.first;m++) //sgRNA iteration
            perm[m] = Uniform(0.0,1.0);
          sort(perm.begin(),perm.end());//all percentiles for one gene sorted
          randLoValue[randLoValueNum] = ComputeLoValue_1tail(perm);
          randLoValueNum++;
        }
      }
    }
    sort(randLoValue.begin(),randLoValue.end());
    for(const auto &j:groups) //gene iteration
      if(groups[j.first].size()==i.first) {
        double index = lower_bound(randLoValue.begin(),randLoValue.end(),loValues[j.first])-randLoValue.begin();
        pValues[j.first] = (2*index+1.0)/(2*randLoValue.size());
      }
  }
  return(pValues);
}

// R FUNCTIONS --------------------------------------------------------------------------------------

// [[Rcpp::export]]
DataFrame RRA_1tail(StringVector genes, StringVector sgrnas, NumericVector values){
  // Set up items and groups:
  set<item,item_order> items; //ordered sgrnas
  map<string,set<item,item_order>> groups; //ordered sgrnas grouped into genes
  string gene,sgrna;
  double value;
  for(int i=0;i<values.size();i++){
    gene = genes(i);
    sgrna = sgrnas(i);
    value = values(i);
    item item1 = item{sgrna,value};
    items.insert(item1); //order sgrnas in real-time
    groups[gene].insert(item1);
  }
  // Get lo values:
  map<string,double> loValues = getLoValuesC_1tail(items, groups);
  // Get p values:
  map<string,double> pValues = getPValuesC_1tail(groups, loValues);
  // Output data frame:
  StringVector names(pValues.size());
  NumericVector sizes(pValues.size()),result(pValues.size()),pvalues1(pValues.size());
  int iter = 0;
  for(const auto &i:pValues){
    names(iter) = i.first;
    sizes(iter) = groups[i.first].size();
    result(iter) = loValues[i.first];
    pvalues1(iter) = pValues[i.first];
    iter++;
  }
  return DataFrame::create(_["Gene"]=names,_["Size"]=sizes,_["loValue"]=result,_["p"]=pvalues1,_["stringsAsFactors"]=false);
}
