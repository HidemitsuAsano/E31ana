
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



void plot_InouepS()
{
  
  ifstream ifs("../pimSp_CS.txt");
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
   
    if(iline!=0 && arr.size()==4){
      mass[iline] = atof(arr.at(0).c_str());
      cout << mass[iline] << "\t";
      CS[iline] = atof(arr.at(1).c_str());
      cout << CS[iline] << "\t";
      binhalf[iline] = atof(arr.at(2).c_str());
      error[iline] = atof(arr.at(3).c_str());
      cout << error[iline] <<endl;
    }

    iline++;

  }

  TGraphErrors *grSp = new TGraphErrors(iline,mass,CS,binhalf,error);
  grSp->SetTitle("Forward Proton");
  grSp->SetMarkerColor(2);
  grSp->SetLineColor(2);
  grSp->SetMarkerStyle(20);
  grSp->Print();
  grSp->GetXaxis()->SetRangeUser(1.35,1.5);
  //grSp->RemovePoint(26);
  //grSp->RemovePoint(0);
  grSp->Draw("AP");
   

  ifstream ifsSm("../pipSm_CS.txt");
  std::string str2;
  
  double massSm[maxline];
  double binhalfSm[maxline];
  double CSSm[maxline];
  double errorSm[maxline];
  int iline2=0;
  while(getline(ifsSm,str2)){
    vector<string> arr = split(str2);

    for(int i=0;i<arr.size();++i){
      cout << arr[i] << "\t" ;
    }
    cout << endl;

    if(iline2!=0 && arr.size()==4){
      massSm[iline2] = atof(arr.at(0).c_str());
      cout << massSm[iline2] << "\t";
      CSSm[iline2] = atof(arr.at(1).c_str());
      cout << CSSm[iline2] << "\t";
      binhalfSm[iline2] = atof(arr.at(2).c_str());
      errorSm[iline2] = atof(arr.at(3).c_str());
      cout << errorSm[iline2] <<endl;
    }

    iline2++;
  }

  TGraphErrors *grSm = new TGraphErrors(iline2,massSm,CSSm,binhalfSm,errorSm);
  grSm->SetTitle("pi+Sigma-");
  grSm->SetMarkerColor(4);
  grSm->SetLineColor(4);
  grSm->SetMarkerStyle(20);
  grSm->Print();
  //grSp->RemovePoint(26);
  //grSp->RemovePoint(0);
  grSm->Draw("P");
   


  TFile *f = new TFile("pSinoue_CS.root","RECREATE");
  grSp->Write();
  grSm->Write();
  f->Close();
}


