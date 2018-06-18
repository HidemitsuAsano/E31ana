#include "AnaInfo.h"

using namespace std;

AnaInfo::AnaInfo()
{
}

CDSInfo *AnaInfo::minDCA()
{
  CDSInfo *minDCA=0;
  double min_dca=DBL_MAX;
  for( int i=0; i<nCDS(); i++ ){
    if( CDS(i)->flag() && CDS(i)->dca()<min_dca ){
      min_dca=CDS(i)->dca();
      minDCA=CDS(i);
    }
  }
  return minDCA;
}

int AnaInfo::nCDS(const int &pid)
{
  int num=0;
  for( int i=0; i<(int)fCDSInfoContainer.size(); i++ ){
    if( fCDSInfoContainer[i].pid()==pid ) num++;
  }
  return num;
}

CDSInfo* AnaInfo::CDSbyID(const int &id)
{
  for( int i=0; i<(int)fCDSInfoContainer.size(); i++ ){
    if( fCDSInfoContainer[i].trackID()==id ) return &fCDSInfoContainer[i];
  }
  return 0;
}

CDSInfo* AnaInfo::CDS(const int &pid, const int &id)
{
  int num=0;
  for( int i=0; i<(int)fCDSInfoContainer.size(); i++ ){
    if( fCDSInfoContainer[i].pid()==pid ){
      if( num==id ) return &fCDSInfoContainer[i];
      num++;
    }
  }
  return 0;
}

int AnaInfo::nCDS2(int pid1, int pid2)
{
  int num=0;
  if( pid1>pid2 ) std::swap(pid1, pid2);
  for( int i=0; i<(int)fCDS2InfoContainer.size(); i++ ){
    int pid1_2=fCDS2InfoContainer[i].pid1(), pid2_2=fCDS2InfoContainer[i].pid2();
    if( pid1_2>pid2_2 ) std::swap(pid1_2, pid2_2);
    if( pid1_2==pid1 && pid2_2==pid2 ) num++;
  }
  return num;
}

CDS2Info* AnaInfo::CDS2(int pid1, int pid2, const int &id)
{
  int num=0;
  if( pid1>pid2 ) std::swap(pid1, pid2);
  for( int i=0; i<(int)fCDS2InfoContainer.size(); i++ ){
    int pid1_2=fCDS2InfoContainer[i].pid1(), pid2_2=fCDS2InfoContainer[i].pid2();
    if( pid1_2>pid2_2 ) std::swap(pid1_2, pid2_2);
    if( pid1_2==pid1 && pid2_2==pid2 ){
      if( num==id ) return &fCDS2InfoContainer[i];
      num++;
    }
  }
  return 0;
}

void AnaInfo::Clear()
{
  fBeamInfoContainer.clear();
  fCDSInfoContainer.clear();
  fCDS2InfoContainer.clear();
  fFCInfoContainer.clear();
  fFNInfoContainer.clear();
}

void AnaInfo::dump()
{
  cout<<"===== AnaInfo::dump ====="<<endl;
  if( fBeamInfoContainer.size()==1 ){
    cout<<"> Beam : ";
    if( beam(0)->pid()==Beam_Kaon ) cout<<"Kaon"<<endl;
    else if( beam(0)->pid()==Beam_Pion     ) cout<<"Pion"<<endl;
    else if( beam(0)->pid()==Beam_Proton   ) cout<<"Proton"<<endl;
    else if( beam(0)->pid()==Beam_Deuteron ) cout<<"Deuteron"<<endl;
    TVector3 pos=beam(0)->T0pos();
    cout<<"D5 mom  : "<<beam(0)->D5mom()<<" [GeV/c2]    vtx : "<<Form("(%lf, %lf, %lf)", pos.X(), pos.Y(), pos.Z())<<endl;
    pos=beam(0)->vertex();
    cout<<"Vtx mom : "<<beam(0)->vertex_mom()<<" [GeV/c2]    vtx : "<<Form("(%lf, %lf, %lf)", pos.X(), pos.Y(), pos.Z())<<endl;
  }
  else{
    cout<<">n Beam : "<<fBeamInfoContainer.size()<<endl;
  }
  cout<<"========================="<<endl;
  cout<<"> CDS pi+ : "<<nCDS(CDS_PiPlus)<<endl;
  cout<<"> CDS p   : "<<nCDS(CDS_Proton)<<endl;
  cout<<"> CDS d   : "<<nCDS(CDS_Deuteron)<<endl;
  cout<<"> CDS He3 : "<<nCDS(CDS_Triton)<<endl;
  cout<<"> CDS t   : "<<nCDS(CDS_Helium3)<<endl;
  cout<<"> CDS pi- : "<<nCDS(CDS_PiMinus)<<endl;
  cout<<"> CDS K-  : "<<nCDS(CDS_Kaon)<<endl;
  cout<<"> CDS Other : "<<nCDS(CDS_Other)<<endl;
  cout<<"> CDS Default : "<<nCDS(CDS_DEFAULT)<<endl;

  if( nFNeutral()==1 ){
    if( forwardNeutral(0)->pid()==F_Gamma ){
      cout<<"> Gamma "<<endl;
    }
    else if( forwardNeutral(0)->pid()==F_Neutron ){
      cout<<"> Neutron "<<endl;
    }
  }
  if( nFCharge()==1 ){
    if( forwardCharge(0)->pid()==F_Pion ){
      cout<<"> Forward Pion "<<endl;
    }
    else if( forwardCharge(0)->pid()==F_Proton ){
      cout<<"> Forward Proton "<<endl;
    }
    else if( forwardCharge(0)->pid()==F_Deuteron ){
      cout<<"> Forward Deuteron "<<endl;
    }
  }
}
