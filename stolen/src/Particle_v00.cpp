#include "Particle.h"
#include "ELossTools.h"
#include "GeomTools.h"

ClassImp(Particle);
ClassImp(pBeam);
ClassImp(pCDS);
ClassImp(pNC);
ClassImp(pPC);

Particle::Particle() : TObject()
{
  Clear();
}

void Particle::Clear(){
  ProductContainer.clear();
  BeamContainer.clear();
  CDSContainer.clear();
  NCContainer.clear();
  PCContainer.clear();
  CVCContainer.clear();

  NCrawContainer.clear();
  PCrawContainer.clear();
  IHrawContainer.clear();
  CDHrawContainer.clear();
  CVCrawContainer.clear();
  BVCrawContainer.clear();
  
  CDSTrackIDContainer.clear();
  KaonTrackIDContainer.clear();
  PiplusTrackIDContainer.clear();
  PiminusTrackIDContainer.clear();
  ProtonTrackIDContainer.clear();
  DeuteronTrackIDContainer.clear();
  Helium3TrackIDContainer.clear();
  TritonTrackIDContainer.clear();
  OtherTrackIDContainer.clear();
}

void Particle::AddBeam( const pBeam &beam )
{
  BeamContainer.push_back(beam);
}

void Particle::AddNC( const pNC &nc )
{
  NCContainer.push_back(nc);
}

void Particle::AddPC( const pPC &pc )
{
  PCContainer.push_back(pc);
}

void Particle::AddCVC( const pPC &pc )
{
  CVCContainer.push_back(pc);
}

//void Particle::AddProduct( const pProduct &pro )
void Particle::AddProduct( const pCDS &pro )
{
  ProductContainer.push_back(pro);
}

void Particle::AddCDS( const pCDS &track )
{
  int pid=track.pid();  
  int id=track.id();  
  CDSContainer.push_back(track);
  CDSTrackIDContainer[ id ] = CDSTrackIDContainer.size()-1;
  switch( pid ){
  case CDS_Kaon:
    KaonTrackIDContainer.push_back(id);
    break;
  case CDS_PiPlus:
    PiplusTrackIDContainer.push_back(id);
    break;
  case CDS_PiMinus:
    PiminusTrackIDContainer.push_back(id);
    break;
  case CDS_Proton:
    ProtonTrackIDContainer.push_back(id);
    break;
  case CDS_Deuteron:
    DeuteronTrackIDContainer.push_back(id);
    break;
  case CDS_Helium3:
    Helium3TrackIDContainer.push_back(id);
    break;
  case CDS_Triton:
    TritonTrackIDContainer.push_back(id);
    break;
  case CDS_Other:
    OtherTrackIDContainer.push_back(id);
    break;
  default:
    break;
  }
  return;
}

void Particle::CalcAngleCM(const double &targetmass){
  // if(nBeam()!=1){
  //   std::cout<<"!!! error in Particle::CalcAngleCM !!!"<<std::endl;
  //   return;
  // }
  // for(int i=0;i<nProduct();i++) 
  //   product(i)->CalcAngleCM(beam(0)->GetLorentzVector(),targetmass);
  // for(int i=0;i<nCDS();i++) 
  //   cdsi(i)->CalcAngleCM(beam(0)->GetLorentzVector(),targetmass);
  return;
}

pCDS* Particle::cds(const int &id){
  unsigned int i;
  std::map<int,int>::iterator is = CDSTrackIDContainer.find(id);
  if( is != CDSTrackIDContainer.end()){
    i=is->second;
  }else{
    return 0;
  }
  return (0<=i&&i<CDSContainer.size()) ? &CDSContainer[i] : 0;
}

pCDS* Particle::kaon(const int &id){
  return (0<=id&&id<(int)KaonTrackIDContainer.size()) ? cds(KaonTrackIDContainer[id]) : 0;
}

pCDS* Particle::pip(const int &id){
  return (0<=id&&id<(int)PiplusTrackIDContainer.size()) ? cds(PiplusTrackIDContainer[id]) : 0;
}

pCDS* Particle::pim(const int &id){
  return (0<=id&&id<(int)PiminusTrackIDContainer.size()) ? cds(PiminusTrackIDContainer[id]) : 0;
}

pCDS* Particle::proton(const int &id){
  return (0<=id&&id<(int)ProtonTrackIDContainer.size()) ? cds(ProtonTrackIDContainer[id]) : 0;
}

pCDS* Particle::deuteron(const int &id){
  return (0<=id&&id<(int)DeuteronTrackIDContainer.size()) ? cds(DeuteronTrackIDContainer[id]) : 0;
}

pCDS* Particle::triton(const int &id){
  return (0<=id&&id<(int)TritonTrackIDContainer.size()) ? cds(TritonTrackIDContainer[id]) : 0;
}

pCDS* Particle::helium3(const int &id){
  return (0<=id&&id<(int)Helium3TrackIDContainer.size()) ? cds(Helium3TrackIDContainer[id]) : 0;
}

pCDS* Particle::other(const int &id){
  return (0<=id&&id<(int)OtherTrackIDContainer.size()) ? cds(OtherTrackIDContainer[id]) : 0;
}

// void pCDS::CalcAngleCM(TLorentzVector beam, const double &targetmass){
//   TLorentzVector target(0.,0.,0.,targetmass);
//   TLorentzVector tmp=GetLorentzVector();

//   TVector3 boost=-1*(beam+target).BoostVector();
//   target.Boost(boost);
//   tmp.Boost(boost);
//   beam.Boost(boost);

//   //  AngleCM=beam.Vect().Angle(tmp.Vect());
//   return;
// }


pBeam::pBeam() : TObject(){
  BLCTRACK=false;
  BPCTRACK=false;
  MOM=false;
  PID=-1;
  T0seg=-1;
  BHDseg=-1;
  tofbhdt0=-999;
  Momentum=-999;
  BHDX=-999;
}

double pBeam::CalcVertexTime(const TVector3 &vertex)
{
  double momout,tof;
  ELossTools::CalcElossBeamTGeo(t0pos(),vertex,mom(),mass(),momout,tof);
  return t0time()+tof;
}
double pBeam::CalcVertexMom(const TVector3 &vertex)
{
  double momout,tof;
  ELossTools::CalcElossBeamTGeo(t0pos(),vertex,mom(),mass(),momout,tof);
  return momout;
}
pNC::pNC() : TObject(), Segment(DEFAULTI), PID(DEFAULTI)
	   , HitPos(DEFVECT), Mom(DEFVECT)
	   , Time(DEFAULTD), Energy(DEFAULTD), Beta(DEFAULTD), Mass(DEFAULTD), TOF(DEFAULTD), FLength(DEFAULTD)
{
}

pNC::pNC(int seg, double time, TVector3 pos) :TObject(), Segment(seg),PID(DEFAULTI),
					      HitPos(pos), Mom(DEFVECT)
					     ,Time(time), Energy(DEFAULTD), Beta(DEFAULTD), Mass(DEFAULTD), TOF(DEFAULTD), FLength(DEFAULTD)
{
}

void pNC::CalcMom(pBeam *beam,const  TVector3 &vertex,double offs){
  double beamtime=beam->CalcVertexTime(vertex); 
  FLength=(vertex-hitpos()).Mag();
  TOF=time()-beamtime-offs;
  double betan=fl()/tof()/(Const*100);
  double momn= nMass/sqrt(1/betan/betan-1); 
  TVector3 mom=(hitpos()-vertex).Unit()*momn;  
  SetBeta(betan);
  SetMomentum(mom);
  if(1./betan>1.05){
    PID=F_Neutron;
    Mass=nMass;
  }
  else{
    PID=F_Gamma;
    Mass=0;
  }
}

pPC::pPC() : pNC(), FDC1Pos(DEFVECT), Angle(DEFAULTD), Radius(DEFAULTD)
{
}

pPC::pPC(int seg, double time, TVector3 pos) :pNC(seg,time,pos), FDC1Pos(DEFVECT), Angle(DEFAULTD), Radius(DEFAULTD)
{
}

double pPC::momUSWK(const TVector3 &vertex)
{
  double momout,tmptof;
  ELossTools::CalcElossBeamTGeo(vertex,fdc1pos(),mom().Mag(),mass(),momout,tmptof);
  return momout;
}

void pPC::CalcMomBVC(const  TVector3 &vertex, int pid, bool SIM){
  TVector3 pdir=(fdc1pos()-vertex).Unit();
  TVector3 mom(DEFVECT);
  SetMomentum(mom);
  if(nbvchit()<1) return ;
  TVector3 uswkpos=fdc1pos()+(250-fdc1pos().Z())/pdir.Z()*pdir;
  Angle=(uswkpos-vertex).Angle(hitpos()-uswkpos);

  TVector3 bvcpos=fdc1pos()+(bvchit(0).pos().Z()-fdc1pos().Z())/pdir.Z()*pdir;
  TOF=time()-bvchit(0).ctmean();
  
  //consider theve difference between the lenght of 2 lines and arc trajectory in USWK
  double leff=100; //cm
  Radius=leff/sqrt(2*(1-TMath::Cos(angle())));
  double larc=r()*angle();

  double line=2*leff/sqrt(2*(1+TMath::Cos(angle())));
  double diff=line-larc;
  
  double flength=(uswkpos-bvcpos).Mag()+(hitpos()-uswkpos).Mag()-diff;
  if(!SIM)  flength-=1.*cm;
  FLength=flength;
  Beta=fl()/tof()/(Const*100);
  double tmpmomr=momr();
  Mass2=tmpmomr*tmpmomr*(1/beta()/beta()-1);
  if(pid<0){
    if(Mass2>(0.02-0.007*5)&&Mass2<(0.02+0.007*5)) PID=F_Pion;
    else if(Mass2>(0.885-0.046*5)&&Mass2<(0.885+0.046*5)) PID=F_Proton;
    else if(Mass2>(3.51-0.3*2)&&Mass2<(3.51+0.3*2))       PID=F_Deuteron;
    else PID=F_Other;
    //    std::cout<<Mass2<<"  "<<PID<<std::endl;
  }else   PID=pid;
  Mass=parMass[PID];
  double mom1= Mass/sqrt(1/beta()/beta()-1);    
  if(beta()>1){ SetMomentum(mom); return; }
  double momraw=mom1;
  double mom2;      
  double tmpmom,tmptof;
  double tmpdiff;
  ELossTools::CalcElossForwardTGeo(bvcpos, fdc1pos(), fl(),mom1, Mass, tmpmom,tmptof);
  tmpdiff=tmptof - tof();
  mom2=mom1;
  int count=0;
  // if(PID==F_Pion)
  //   std::cout<<Mass2<<"  "<<mom1<<std::endl;
  do{
    mom2+=0.01*tmpdiff/fabs(tmpdiff);
    if(!ELossTools::CalcElossForwardTGeo(bvcpos, fdc1pos(), fl(),mom2, Mass, tmpmom,tmptof)) return;
    diff=tmptof - tof();
    count++;
    if(count>50){
      std::cout<<" too many loops !!!  "<<tmpdiff<<"  "<<mom1<<"  "<<mom2<<"  "<<tmpdiff<<"  "<<diff<<std::endl;
      return; 
    }
  } while (tmpdiff*diff>0);
  if(mom1>mom2){
    double tmp=mom1;
    mom1=mom2;
    mom2=tmp;
  }
  do{
    double mommid=(mom1+mom2)*0.5;
    if(!ELossTools::CalcElossForwardTGeo(bvcpos, fdc1pos(), fl(),mommid, Mass, tmpmom,tmptof))return;
    diff=tmptof - tof();
    if(diff>0) mom1=mommid;
    else mom2=mommid;
    count++;
    if(count>50) return; 
  } while (TMath::Abs(mom2-mom1)>0.0001);
  mom=(fdc1pos()-vertex).Unit()*mom1;
  SetMomentum(mom);
}
void pPC::CalcMom(pBeam *beam,const  TVector3 &vertex, int pid, bool SIM){
  double beamtime=beam->CalcVertexTime(vertex); 
  //  TVector3 fdc1  
  TVector3 pdir=(fdc1pos()-vertex).Unit();
  TVector3 uswkpos=fdc1pos()+(250-fdc1pos().Z())/pdir.Z()*pdir;
  Angle=(uswkpos-vertex).Angle(hitpos()-uswkpos);

  TOF=time()-beamtime;
  
  //consider theve difference between the lenght of 2 lines and arc trajectory in USWK
  double leff=100; //cm
  Radius=leff/sqrt(2*(1-TMath::Cos(angle())));
  double larc=r()*angle();

  double line=2*leff/sqrt(2*(1+TMath::Cos(angle())));
  double diff=line-larc;
  
  double flength=(uswkpos-vertex).Mag()+(hitpos()-uswkpos).Mag()-diff;
  if(!SIM)  flength-=1.5*cm;
  //  flength+=5.0*cm;
  FLength=flength;
  Beta=fl()/tof()/(Const*100);
  double tmpmomr=momr();
  Mass2=tmpmomr*tmpmomr*(1/beta()/beta()-1);
  if(pid<0){
    if(Mass2>(0.02-0.007*5)&&Mass2<(0.02+0.007*5)) PID=F_Pion;
    else if(Mass2>(0.885-0.046*5)&&Mass2<(0.885+0.046*5)) PID=F_Proton;
    else if(Mass2>(3.51-0.3*2)&&Mass2<(3.51+0.3*2))       PID=F_Deuteron;
    else PID=F_Other;
    //    std::cout<<Mass2<<"  "<<PID<<std::endl;
  }else   PID=pid;
  Mass=parMass[PID];
  
  double mom1= Mass/sqrt(1/beta()/beta()-1);    
  TVector3 mom(DEFVECT);
  if(beta()>1){ SetMomentum(mom); return; }
  if(1){
    double momraw=mom1;
    double mom2;      
    double tmpmom,tmptof;
    double tmpdiff;
    ELossTools::CalcElossForwardTGeo(vertex, fdc1pos(), fl(),mom1, Mass, tmpmom,tmptof);
    tmpdiff=tmptof - tof();
    mom2=mom1;
    int count=0;
    // if(PID==F_Pion)
    //   std::cout<<Mass2<<"  "<<mom1<<std::endl;
    do{
      mom2+=0.01*tmpdiff/fabs(tmpdiff);
      if(!ELossTools::CalcElossForwardTGeo(vertex, fdc1pos(), fl(),mom2, Mass, tmpmom,tmptof)) return;
      diff=tmptof - tof();
      count++;
      if(count>50){
	std::cout<<" too many loops !!!  "<<tmpdiff<<"  "<<mom1<<"  "<<mom2<<"  "<<tmpdiff<<"  "<<diff<<std::endl;
	return; 
      }
    } while (tmpdiff*diff>0);
    if(mom1>mom2){
      double tmp=mom1;
      mom1=mom2;
      mom2=tmp;
    }
    do{
      double mommid=(mom1+mom2)*0.5;
      if(!ELossTools::CalcElossForwardTGeo(vertex, fdc1pos(), fl(),mommid, Mass, tmpmom,tmptof))return;
      diff=tmptof - tof();
      if(diff>0) mom1=mommid;
      else mom2=mommid;
      count++;
      if(count>50) return; 
    } while (TMath::Abs(mom2-mom1)>0.0001);
    // if(PID==F_Pion)
    //   std::cout<<"fl,tof,beta,mass2,momraw,momcorr   "<<fl()<<"  "<<tof()<<"  "<<1./beta()<<"  "<<Mass2<<"  "<<momraw<<"  "<<mom1<<std::endl;
    mom=(fdc1pos()-vertex).Unit()*mom1;
    SetMomentum(mom);
  }
}

pCDS::pCDS() : TObject(){
  trackID=-1;
  daughterID1=-1;
  daughterID2=-1;
  CombID=-1;
  PID=-1;
  Momentum=-1;
  Beta=-1;
  TOF=-999;
  Mass=-1;
  PDGMass=-1;
  VertexDistance=-1;
  VertexBeamDistance=-999;
  ProductBeamDCA=-999;
  //  for(int i=0;i<5;i++) Param[i]=-999.;
  CDHseg.clear();
  IHseg.clear();
  AngleLab=-999.;
  OpenAngle=-999;
  //  AngleCM=-999.;
  FlightLength=-999;
}


double pCDS::mom(ConfMan *conf){
  CDSFittingParamMan *CDSParam=conf->GetCDSFittingParamManager();
  double magneticfield=CDSParam->GetMagneticField();
  double tmppt=magneticfield*Const/100./Param[2];
  return tmppt*sqrt(1+Param[4]*Param[4]);
}

// double pCDS::cmom2(const TVector3 &vertex){
//   double cmom2;
//   double tof;
//   ELossTools::CalcHelixElossToVertexTGeo(Param,vertex,mom(),pdgmass(),cmom2,tof);
//   // for(int i=0;i<5;i++)
//   //   std::cout<<Param[i]<<"  ";
//   // std::cout<<std::endl;
//   //  std::cout<<pdgmass()<<"  "<<mom()<<" -> "<<cmom2<<"  "<<cmom()<<std::endl;
//   return cmom2;
// }

// double pCDS::cmom2(){
//   return cmom2(vcdc());
// }

pProduct::pProduct() : TObject(){
  trackID1=-1;
  trackID2=-1;
  CombID=-1;
  PID=-1;
  Momentum=-999;
  Mass=-1;
  PDGMass=-1;
  VertexDistance=-999;
  VertexBeamDistance=-999;
  ProductBeamDCA=-999;
  Beta=-1;
  Gamma=-1;
  TOF=-999;
  OpenAngle=-999;
  AngleLab=-999;
  AngleCM=-999;
}

void pProduct::SetPDGMass(const double &mass)
{
  
}

void pProduct::CalcAngleCM(TLorentzVector beam, const double &targetmass){
  TLorentzVector target(0.,0.,0.,targetmass);
  TVector3 boost=-1*(beam+target).BoostVector();
  TLorentzVector tmp=GetLorentzVector();
  TLorentzVector ctmp=GetCLorentzVector();

  tmp.Boost(boost);
  ctmp.Boost(boost);
  beam.Boost(boost);
  target.Boost(boost);

  AngleCM=beam.Vect().Angle(tmp.Vect());
  CAngleCM=beam.Vect().Angle(ctmp.Vect());
}

void Particle::CalcNCMom()
{
  if(nBeam()!=1) return;
  TVector3 vertex=beam(0)->vertex();
  CalcNCMom(vertex);
}
void Particle::CalcNCMom(const TVector3 &vertex)
{
  //  if(nNC()!=1) return;
  if(nBeam()!=1) return;
  for(int i=0;i<nNC();i++){
    nc(i)->CalcMom(beam(0),vertex);
  }
}

HodoscopeLikeHit* Particle::ncrawseg(const int &seg)
{
  for( hodoContainer::iterator it=NCrawContainer.begin();
       it!=NCrawContainer.end(); it++ ){
    if(it->seg()==seg) return &(*it);
  }
  return 0;
}
