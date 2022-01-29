//H. Asano
//This is a code to read a text file "run62.list" 
//and make a good run list for RUN62

#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;


int main(){
  
  ifstream ifs("run62.list");
  std::string str;
  
  const int maxline=1024;
  int runnum[maxline];
  
  int iline=0;
  while(getline(ifs,str)){
    
    //for(int i=0;i<arr.size();++i){
    //  //cout << arr[i] << "\t";
    //}
    //std::cout << std::endl;
    runnum[iline] = atoi(str.c_str());
    //std::cout << arr.size() << std::endl;
    
    iline++;
  }

  
  int pro_nruns=0;
  //generate shell script for hadd
  ofstream ofs_pro_haddcsh("hadd_h2.csh");
  ofs_pro_haddcsh << "#!/bin/tcsh -f" << endl;

  ofs_pro_haddcsh << "hadd -To evanaIMsigma " ;
  for(int i=0; i<iline;i++){
    char fname[256];
    cout << runnum[i] << endl; ;
    sprintf(fname,"evanaIMsigma_h2_0%03d.root",runnum[i]);
    ofs_pro_haddcsh << fname << " "  ;
  }

  cout << endl;
}


