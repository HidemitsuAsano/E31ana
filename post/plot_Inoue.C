
#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

vector <string> split(const string& src, const char* delimiter = "\t"){
  vector <string> wordvec;
  string::size_type len = src.length();

  for(string::size_type i=0,n; i< len;i=n+1){
    n =src.find_first_of(delimiter,i);
    if( n== string::npos){
      n = len;
    }
    wordvec.push_back(src.substr(i,n-i));
  }

  return wordvec;
}



void plot_Inoue()
{
  
  ifstream ifs("pimL_CS.txt");
  std::string str;
  
  const int maxline=1024;
  double mass[maxline];
  double binhalf[maxline];
  double CS[maxline];
  double error[maxline];
  
  int iline=0;
  while(getline(ifs,str)){
    vector<string> arr = split(str);

    for(int i=0;i<arr.size();++i){
      cout << arr[i] << "\t" ;
    }
    cout << endl;

    if(iline!=0 && arr.size()==5){
      mass[iline] = atof(arr.at(0).c_str());
      cout << mass[iline] << "\t";
      CS[iline] = atof(arr.at(1).c_str());
      cout << CS[iline] << "\t";
      binhalf[iline] = atof(arr.at(2).c_str());
      error[iline] = atof(arr.at(4).c_str());
      cout << error[iline] <<endl;
    }

    iline++;

  }

  TGraphErrors *gr = new TGraphErrors(iline,mass,CS,binhalf,error);
  gr->SetTitle("Forward Proton");
  gr->SetMarkerColor(4);
  gr->SetLineColor(4);
  gr->SetMarkerStyle(20);
  gr->Print();
  gr->RemovePoint(26);
  gr->RemovePoint(0);
  gr->Draw("AP");
  
  TFile *f = new TFile("pimL_CS.root","RECREATE");
  gr->Write();
  f->Close();
}


