#include "genomeProcessing.h"
#include <fstream>
#include "misc.h"
#include "globalVariables.h"

string twoBitsToString(bool b1, bool b2) {
  if (b1) {
    if (b2) {
      return "T";
    } else {
      return "G";
    }
  } else {
    if (b2) {
      return "C";
    } else {
      return "A";
    }
  }
}

std::vector<bool> fullstr2num(const string& str) {
  std::vector<bool> res;
  for(uint i(0);i<str.size();i++){
    switch (str[i]){
      case 'A':res.push_back(false);res.push_back(false);break;
      case 'C':res.push_back(false);res.push_back(true);break;
      case 'G':res.push_back(true);res.push_back(false);break;
      default:res.push_back(true);res.push_back(true);break;
    }
  }
  return res;
}

std::string fullnum2str(vector<bool> num) {
  string str(num.size()/2, 'N');
  uint j = 0;
  for(uint i(0);i<num.size();i+=2){
    if(num[i]){
      if(num[i+1]){
      	str[j] = 'T';
      }else{
        str[j] = 'G';
      }
    }else{
      if(num[i+1]){
        str[j] = 'C';
      }else{
        str[j] = 'A';
      }
    }
    j++;
  }
  return str;
}

robin_hood::unordered_map<std::string, std::vector<bool>> indexGenome(std::string readsFile) {
  robin_hood::unordered_map<std::string, std::vector<bool>> index;

  std::ifstream f(readsFile);
  std::string header, sequence, seq;
  int nbLines = 0;

  getline(f, header);
  while (header.length() > 0) {
    // Get header
    header.erase(0, 1);
    header = splitString(header, " ")[0];
    
    // Get sequence, watching out for multiline FASTA/FASTQ
    getline(f, seq);
    sequence = seq;
    nbLines = 1;
    getline(f, seq);
    while (seq.length() > 0 and seq[0] != '>' and seq[0] != '+') {
      sequence += seq;
      nbLines++;
      getline(f, seq);
    }

    // Index header/sequence pair
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
    index[header] = fullstr2num(sequence);

    // Skip remaining lines if FASTQ
    if (seq[0] == '+') {
      getline(f, seq);
      for (int i = 1; i < nbLines; i++) {
        getline(f, seq);
      }
      getline(f, seq);
    }

    // Next header has been read, update and loop
    header = seq;
  }

  return index;
}

string extractGenomeRegion(string region) {
  vector<string> t = splitString(region, ":");
  vector<string> tt = splitString(t[1], "-");

  return fullnum2str(genomeIndex[t[0]]).substr(static_cast<uint32_t>(stoul(tt[0])), windowSize);

}