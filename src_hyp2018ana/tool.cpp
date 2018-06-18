#include <iostream>
#include <fstream>

#include "kn_El.h"
#include "TTree.h"
#include "KnuclRootData.h"

int main(int argc, char** argv)
{
  std::cout<<"===== tool START  ====="<<std::endl;
  if( argc!=2 ){
    std::cout<<"Please input only 1 file"<<std::endl;
    return 0;
  }

  std::ifstream ifs(argv[1]);
  std::string str;
  std::string confname;
  std::string outfilename;
  std::string indir;
  int min, max;
  while( getline(ifs, str) ){
    char c_str[500];
    if( sscanf(str.c_str(), "ConfFile: %s", c_str)==1 ){
      confname=c_str;
    }
    if( sscanf(str.c_str(), "OutFile: %s", c_str)==1 ){
      outfilename=c_str;
    }
    if( sscanf(str.c_str(), "Input: %s %d %d", c_str, &min, &max)==3 ){
      indir=c_str;
    }
  }
  std::cout<<"> ConfFile: "<<confname<<std::endl;
  std::cout<<"> OutFile: "<<outfilename<<std::endl;
  std::cout<<"> Input: "<<indir<<" "<<min<<" "<<max<<std::endl;

  ConfMan *conf = new ConfMan(confname);
  conf-> Initialize();

  TFile *of = new TFile(outfilename.c_str(), "recreate");
  initHist();
  for( int i=min; i<max; i++ ){
     std::cout<<"> File : "<<Form("%s/%03d/test_mc.root", indir.c_str(), i)<<std::endl;
     TFile *infile = new TFile(Form("%s/%03d/test_mc.root", indir.c_str(), i));
     //     TFile *infile = new TFile("~/test/geant/knucl6/build/test_mc.root");
     if( !infile->IsOpen() || infile->IsZombie() ) continue;
     TTree *tree2 = (TTree*)infile-> Get("tree2");
     TTree *tree = (TTree*)infile-> Get("tree");
     if( !tree ) continue;
     DetectorData *detData =0;
     MCData *mcData =0;
     ReactionData *reacData =0;
     tree-> SetBranchAddress("DetectorData", &detData);
     tree-> SetBranchAddress("MCData", &mcData);
     tree-> SetBranchAddress("ReactionData", &reacData);
    
     for( int ev=0; ev<tree->GetEntries(); ev++ ){
       tree-> GetEntry(ev);
       of-> cd();
       fillHist(detData, mcData, reacData);
     }
     delete infile;
  }
  of-> cd();
  writeHist();

  std::cout<<"===== tool FINISH ====="<<std::endl;
  return 0;
}
