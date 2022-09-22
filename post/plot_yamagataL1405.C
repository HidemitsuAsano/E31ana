
#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

vector <string> split(const string& src, const char* delimiter = " "){
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



void plot_yamagataL1405()
{
  
  ifstream ifs("yamagata_plot.csv");
  std::string str;
  
  const int maxline=1024;
  double cosval[maxline];
  double CS[maxline];
  
  int iline=0;
  while(getline(ifs,str)){
    //std::cout << str << std::endl;
    vector<string> arr = split(str);

    for(int i=0;i<arr.size();++i){
      cout << arr[i] << "\t" ;
    }
    cout << endl;
   
    if( arr.size()==2){
      cosval[iline] = atof(arr.at(0).c_str());
      cout << cosval[iline] << "\t";
      CS[iline] = atof(arr.at(1).c_str());
      cout << CS[iline] << "\t";
      //error[iline] = sqrt(staterr[iline]*staterr[iline]);//+shapeerr[iline]*shapeerr[iline]);//+cserr[iline]*cserr[iline]);
    }

    iline++;

  }
  
  TGraph *grL1405 = new TGraph(iline,cosval,CS);
  grL1405->SetTitle("cos");
  grL1405->SetName("gr_yamagata");
  grL1405->SetMarkerColor(2);
  grL1405->SetLineColor(2);
  grL1405->SetMarkerStyle(20);
  grL1405->Print();
  //grL1405->GetXaxis()->SetRangeUser(1.35,1.5);
  //grSpcs->SetFillColor(0);
  //grSpcs->SetFillStyle(0);
  //grSp->Draw("AP");
  grL1405->Draw("ap");
  
  TFile *f = new TFile("yamagataL1405.root","RECREATE");
  grL1405->Write();
  f->Close();
}

