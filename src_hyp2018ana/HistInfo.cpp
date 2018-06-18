#include "HistInfo.h"

using namespace std;

HistInfo::HistInfo()
{
}

HistInfo::HistInfo(std::vector<double> bin, int n)
  : fBin(bin), fUnder(n, 0.0), fOver(n, 0.0), fName(n), fVal(bin.size()-1, std::vector<double>(n))
{
}

HistInfo::HistInfo(std::vector<double> bin, std::vector<TString> names)
  : fBin(bin), fUnder(names.size(), 0.0), fOver(names.size(), 0.0), fName(names), fVal(bin.size()-1, std::vector<double>(names.size()))
{
}

int HistInfo::index(double val)
{
  if( val<fBin[0] ) return -1;

  for( int i=0; i<fBin.size();  i++ ){
    if( val<fBin[i+1] ) return i;
  }
  return -999;
}

double HistInfo::val(int i1, TString name) const
{
  int j=-1;
  for( int i=0; i<fName.size(); i++ ){
    if( name==fName[i]) j=i;
  }
  if( j<0 ){
    cout<<"  !!! HistInfo::val  name="<<name<<" not found !!!"<<endl;
    return -1e100;
  }

  return fVal.at(i1).at(j);
}

bool HistInfo::setVal(int i1, TString name, double val)
{
  int j=-1;
  for( int i=0; i<fName.size(); i++ ){
    if( name==fName[i]) j=i;
  }
  if( j<0 ){
    cout<<"  !!! HistInfo::setVal  name="<<name<<" not found !!!"<<endl;
    return false;
  }
  fVal.at(i1).at(j)=val;

  return true;
}

void HistInfo::dump()
{
  string str="min   max   ";
  for( int i=0; i<fName.size(); i++ ){
    str+=Form("%8s  ", fName[i].Data());
  }
  cout<<str<<endl;

  for( int i=0; i<fVal.size(); i++ ){
    string str2=Form("%5.3f %5.3f ", fBin[i], fBin[i+1]);
    for( int j=0; j<fVal[i].size(); j++ ) str2+=Form("%8g  ", fVal[i][j]);
    cout<<str2<<endl;
  }
}
