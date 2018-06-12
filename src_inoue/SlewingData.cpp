#include "SlewingData.h"

SlewingData::SlewingData()
{
}

SlewingData::SlewingData(const int cid1, const int cid2) : CID1(cid1), CID2(cid2)
{
  clear();
}

void SlewingData::Set(HodoscopeLikeHit *hit1,HodoscopeLikeHit *hit2)
{
  Seg1=hit1->seg(), Eu1=hit1->eu(), Ed1=hit1->ed(), Emean1=hit1->emean(), Tu1=hit1->tu(), Td1=hit1->td(), Ctu1=hit1->ctu(), Ctd1=hit1->ctd(), Ctmean1=hit1->ctmean();
  Seg2=hit2->seg(), Eu2=hit2->eu(), Ed2=hit2->ed(), Emean2=hit2->emean(), Tu2=hit2->tu(), Td2=hit2->td(), Ctu2=hit2->ctu(), Ctd2=hit2->ctd(), Ctmean2=hit2->ctmean();
}

void SlewingData::Calc(const double par[4], const double offset[4])
{
  Ctu1 = Tu1+par[0]/sqrt(Eu1)-offset[0];
  Ctd1 = Td1+par[1]/sqrt(Ed1)-offset[1];
  Ctmean1 =0.5*(Ctu1+Ctd1);

  Ctu2 = Tu2+par[2]/sqrt(Eu2)-offset[2];
  Ctd2 = Td2+par[3]/sqrt(Ed2)-offset[3];
  Ctmean2 =0.5*(Ctu2+Ctd2);
}

void SlewingData::Calc(const double par[4], const double par1[4], const double offset[4])
{
  Ctu1 = Tu1+par[0]/sqrt(Eu1)+par1[0]*Eu1-offset[0];
  Ctd1 = Td1+par[1]/sqrt(Ed1)+par1[1]*Ed1-offset[1];
  Ctmean1 =0.5*(Ctu1+Ctd1);

  Ctu2 = Tu2+par[2]/sqrt(Eu2)+par1[2]*Eu2-offset[2];
  Ctd2 = Td2+par[3]/sqrt(Ed2)+par1[3]*Ed2-offset[3];
  Ctmean2 =0.5*(Ctu2+Ctd2);
}

void SlewingData::clear()
{
  BeamPID = Beam_Other;
  Seg1=-1, Eu1=0, Ed1=0, Emean1=0, Tu1=0, Td1=0, Ctu1=0, Ctd1=0, Ctmean1=0;
  Seg2=-1, Eu2=0, Ed2=0, Emean2=0, Tu2=0, Td2=0, Ctu2=0, Ctd2=0, Ctmean2=0;
}

void SlewingData::SetOffset(const int &id, const double &offset)
{
  if( id==1 ){
    Tu1 -= offset, Td1 -= offset, Ctu1 -= offset, Ctd1 -= offset, Ctmean1 -= offset;
  }
  else if( id==2 ){
    Tu2 -= offset, Td2 -= offset, Ctu2 -= offset, Ctd2 -= offset, Ctmean2 -= offset;
  }
  else{
    std::cout<<"  !!! SlewingData::SetOffset id="<<id<<" !!!"<<std::endl;
    exit(-1);
  }
}
