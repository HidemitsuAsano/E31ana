
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
  
  ifstream ifs("../CS_pimSp_after.txt");
  std::string str;
  
  const int maxline=1024;
  double mass[maxline];
  double binhalf[maxline];
  double CS[maxline];
  double staterr[maxline];
  double shapeerr[maxline];
  double cserr[maxline];//may be correrated ?
  double error[maxline];//
  
  int iline=0;
  while(getline(ifs,str)){
    vector<string> arr = split(str);

    for(int i=0;i<arr.size();++i){
      cout << arr[i] << "\t" ;
    }
    cout << endl;
   
    if(iline!=0 && arr.size()==6){
      mass[iline] = atof(arr.at(0).c_str());
      cout << mass[iline] << "\t";
      CS[iline] = atof(arr.at(1).c_str());
      cout << CS[iline] << "\t";
      binhalf[iline] = atof(arr.at(2).c_str());
      staterr[iline] = atof(arr.at(3).c_str());
      cout << staterr[iline] <<endl;
      shapeerr[iline] = atof(arr.at(4).c_str());
      cout << shapeerr[iline] <<endl;
      cserr[iline] = atof(arr.at(5).c_str());
      cout << cserr[iline] <<endl;
      //error[iline] = sqrt(staterr[iline]*staterr[iline]);//+shapeerr[iline]*shapeerr[iline]);//+cserr[iline]*cserr[iline]);
      error[iline] = sqrt(staterr[iline]*staterr[iline]+shapeerr[iline]*shapeerr[iline]);
    }

    iline++;

  }

  TGraphErrors *grSp = new TGraphErrors(iline,mass,CS,binhalf,error);
  TGraphErrors *grSpstat = new TGraphErrors(iline,mass,CS,binhalf,staterr);
  TGraphErrors *grSpshape = new TGraphErrors(iline,mass,CS,binhalf,shapeerr);
  TGraphErrors *grSpcs = new TGraphErrors(iline,mass,CS,binhalf,cserr);
  grSp->SetTitle("pi-Sigma+");
  grSp->SetName("gr_inoueSp");
  grSpstat->SetName("gr_inoueSpstat");
  grSpshape->SetName("gr_inoueSpshape");
  grSpcs->SetName("gr_inoueSpcs");
  grSp->SetMarkerColor(2);
  grSp->SetLineColor(2);
  grSp->SetMarkerStyle(20);
  grSp->Print();
  grSp->GetXaxis()->SetRangeUser(1.35,1.5);
  //grSp->RemovePoint(26);
  //grSp->RemovePoint(0);
  grSpcs->SetFillColor(0);
  grSpcs->SetFillStyle(0);
  grSp->Draw("AP");
  grSpcs->Draw("5");
   

  ifstream ifsSm("../CS_pipSm_after.txt");
  std::string str2;
  
  double massSm[maxline];
  double binhalfSm[maxline];
  double CSSm[maxline];
  double staterrSm[maxline];
  double shapeerrSm[maxline];
  double cserrSm[maxline];
  double errorSm[maxline];
  int iline2=0;
  while(getline(ifsSm,str2)){
    vector<string> arr = split(str2);

    for(int i=0;i<arr.size();++i){
      cout << arr[i] << "\t" ;
    }
    cout << endl;

    if(iline2!=0 && arr.size()==6){
      massSm[iline2] = atof(arr.at(0).c_str());
      cout << massSm[iline2] << "\t";
      CSSm[iline2] = atof(arr.at(1).c_str());
      cout << CSSm[iline2] << "\t";
      binhalfSm[iline2] = atof(arr.at(2).c_str());
      staterrSm[iline2] = atof(arr.at(3).c_str());
      cout << staterrSm[iline2] <<endl;
      shapeerrSm[iline2] = atof(arr.at(4).c_str());
      cout << shapeerrSm[iline2] <<endl;
      cserrSm[iline2] = atof(arr.at(5).c_str());
      cout << cserrSm[iline2] <<endl;
      //errorSm[iline2] = sqrt(staterrSm[iline2]*staterrSm[iline2]+shapeerrSm[iline2]*shapeerrSm[iline2]+cserrSm[iline2]*cserrSm[iline2]);
      errorSm[iline2] = sqrt(staterrSm[iline2]*staterrSm[iline2]+shapeerrSm[iline2]*shapeerrSm[iline2]);
    }

    iline2++;
  }

  TGraphErrors *grSm = new TGraphErrors(iline2,massSm,CSSm,binhalfSm,errorSm);
  TGraphErrors *grSmstat = new TGraphErrors(iline2,massSm,CSSm,binhalfSm,staterrSm);
  TGraphErrors *grSmshape = new TGraphErrors(iline2,massSm,CSSm,binhalfSm,shapeerrSm);
  TGraphErrors *grSmcs = new TGraphErrors(iline2,massSm,CSSm,binhalfSm,cserrSm);
  grSm->SetTitle("pi+Sigma-");
  grSm->SetName("gr_inoueSm");
  grSmstat->SetName("gr_inoueSmstat");
  grSmshape->SetName("gr_inoueSmshape");
  grSmcs->SetName("gr_inoueSmcs");
  grSm->SetMarkerColor(4);
  grSm->SetLineColor(4);
  grSm->SetMarkerStyle(20);
  grSm->Print();
  //grSp->RemovePoint(26);
  //grSp->RemovePoint(0);
  grSm->Draw("P");
  grSmcs->SetFillColor(0);
  grSmcs->SetFillStyle(0);
  grSmcs->Draw("5");
   

  TFile *f = new TFile("pSinoue_CS.root","RECREATE");
  grSp->Write();
  grSpstat->Write();
  grSpshape->Write();
  grSpcs->Write();
  grSm->Write();
  grSmstat->Write();
  grSmshape->Write();
  grSmcs->Write();
  f->Close();
}


