#include "ProtonArm.h"

using namespace std;
static const double DELTA_DT=0.1;

namespace ProtonArm{
  static TMinuit* minuit=0;
  static double SCALE;
  static TVector3 VERTEX_POS;
  static double START_Z;
  static double MASS;
  static TVector3 FDC1_POS;
  static HodoscopeLikeHit *COUNTER_HIT=0;
  void SetCounterHit(HodoscopeLikeHit* hit);
  static ConfMan *CONFMAN=0;

  //  static double DT;

  static int USWK_NX;
  static int USWK_NY;
  static int USWK_NZ;
  static double USWK_DX;
  static double USWK_DY;
  static double USWK_DZ;
  static double USWK_XMAX;
  static double USWK_YMAX;
  static double USWK_ZMAX;
  static TVector3 USWK_GPOS;
  static std::map<unsigned long, TVector3> USWK_FIELD_MAP;

  inline double fieldRange_minX(){ return USWK_GPOS.X()-USWK_XMAX; }
  inline double fieldRange_maxX(){ return USWK_GPOS.X()+USWK_XMAX; }
  inline double fieldRange_minY(){ return USWK_GPOS.Y()-USWK_YMAX; }
  inline double fieldRange_maxY(){ return USWK_GPOS.Y()+USWK_YMAX; }
  inline double fieldRange_minZ(){ return USWK_GPOS.Z()-USWK_ZMAX; }
  inline double fieldRange_maxZ(){ return USWK_GPOS.Z()+USWK_ZMAX; }
}

static void fcnFC(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //  cout<<"===== Call fcnBT("<<npar<<")"<<endl;
  //  TVector3 pos(par[0], par[1], ProtonArm::VERTEX_POS.Z());
  //  double diffVtx=(pos-ProtonArm::VERTEX_POS).Mag();
  TVector3 pos=ProtonArm::VERTEX_POS;
  TVector3 mom(par[0], par[1], par[2]);

  ChargeParticle particle(ProtonArm::MASS, pos, mom);

  if( particle.pos().Z()<ProtonArm::fieldRange_minZ() ){
    //    if( !particle.calcELossBeamTGeo(ProtonArm::fieldRange_minZ()) ){
    if( !particle.calcELossBeamTGeo(ProtonArm::fieldRange_minZ(), ProtonArm::CONFMAN->GetCDSFittingParamManager()->GetMagneticField()) ){
      //      cout<<"  !!! calcELossBeamTGeo return false 0 !!!"<<endl;
      f=10e10;
      return;
    }
  }
  double diffFDC1=-1.0;
  while( ProtonArm::isUSWK_Field(particle.pos()) ){
    TVector3 pos=particle.pos();
    if( !particle.nextStepRK4(DELTA_DT) ){
      //      cout<<"  !!! nextStep Runge Kutta return false  !!!"<<endl;
      f=10e10;
      return;
    }
    if( pos.Z()<ProtonArm::FDC1_POS.Z() && ProtonArm::FDC1_POS.Z()<particle.pos().Z() ){
      double rate=(ProtonArm::FDC1_POS.Z()-pos.Z())/(particle.pos().Z()-pos.Z());
      TVector3 FDC1pos=(1.0-rate)*pos+rate*particle.pos();
      particle.SetFDC1pos(FDC1pos);
      diffFDC1=(ProtonArm::FDC1_POS-FDC1pos).Mag();
      // cout<<Form("(%lf, %lf, %lf)", pos.x(), pos.y(), pos.z())<<endl;
      // cout<<Form("(%lf, %lf, %lf)", particle.pos().x(), particle.pos().y(), particle.pos().z())<<endl;
      // cout<<Form("(%lf, %lf, %lf)", FDC1pos.x(), FDC1pos.y(), FDC1pos.z())<<endl;
      // cout<<Form("(%lf, %lf, %lf)", ProtonArm::FDC1_POS.x(), ProtonArm::FDC1_POS.y(), ProtonArm::FDC1_POS.z())<<endl;
    }
  }
  if( diffFDC1<0. ){
    //    cout<<"  !!! cann't find FDC1 pos !!!"<<endl;
    TVector3 pos=particle.pos();
    // cout<<Form("(%lf, %lf, %lf)", pos.x(), pos.y(), pos.z())<<endl;
    // cout<<Form("(%lf, %lf, %lf)", ProtonArm::FDC1_POS.x(), ProtonArm::FDC1_POS.y(), ProtonArm::FDC1_POS.z())<<endl;
    f=10e10;
    return;
  }

  TVector3 counter_hit_pos=ProtonArm::crossPointPCCVC(particle.pos(), particle.mom());
  double diffCounter=ProtonArm::disPCCVC(particle.pos(), particle.mom());
  particle.SetDiffCounter(diffCounter);

  if( !particle.calcELossBeamTGeo(counter_hit_pos.Z()) ){
    //    cout<<"  !!! calcELossBeamTGeo return false 2 !!!"<<endl;
    f=10e10;
    return;
  }
  //  f=diffVtx+diffFDC1+diffCounter;
  f=diffFDC1+fabs(diffCounter);
  // cout<<"diff Vtx : "<<diffVtx<<"  diffFDC1 : "<<diffFDC1<<"  diffCounter : "<<diffCounter<<endl;
  // cout<<"Factor : "<<f<<endl;
}

static void fcnBT(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //  cout<<"===== Call fcnBT("<<npar<<")"<<endl;
  TVector3 pos(par[0], par[1], ProtonArm::START_Z);
  TVector3 mom(par[2], par[3], par[4]);
  TVector3 BPCpos=pos+((ProtonArm::VERTEX_POS.Z()-pos.Z())/mom.Z())*mom;
  double diffBPC=(ProtonArm::VERTEX_POS-BPCpos).Mag();

  ChargeParticle particle(ProtonArm::MASS, pos, mom);
  if( particle.pos().Z()<ProtonArm::fieldRange_minZ() ){
    if( !particle.calcELossBeamTGeo(ProtonArm::fieldRange_minZ()) ){
      //      cout<<"  !!! calcELossBeamTGeo return false !!!"<<endl;
      f=10e10;
      return;
    }
  }
  double diffFDC1=-1.0;
  while( ProtonArm::isUSWK_Field(particle.pos()) ){
    TVector3 pos=particle.pos();
    if( !particle.nextStepRK4(DELTA_DT) ){
      //      cout<<"  !!! nextStep Runge Kutta return false  !!!"<<endl;
      f=10e10;
      return;
    }
    if( pos.Z()<ProtonArm::FDC1_POS.Z() && ProtonArm::FDC1_POS.Z()<particle.pos().Z() ){
      double rate=(ProtonArm::FDC1_POS.Z()-pos.Z())/(particle.pos().Z()-pos.Z());
      TVector3 FDC1pos=(1.0-rate)*pos+rate*particle.pos();
      diffFDC1=(ProtonArm::FDC1_POS-FDC1pos).Mag();
    }
  }
  if( diffFDC1<0. ){
    //    cout<<"  !!! cann't find FDC1 pos !!!"<<endl;
    f=10e10;
    return;
  }

  TVector3 counter_hit_pos=ProtonArm::crossPointPCCVC(particle.pos(), particle.mom());
  if( !particle.calcELossBeamTGeo(counter_hit_pos.Z()) ){
    //    cout<<"  !!! calcELossBeamTGeo return false !!!"<<endl;
    f=10e10;
    return;
  }
  double diffCounter=ProtonArm::disPCCVC(particle.pos(), particle.mom());
  f=diffBPC+diffFDC1+fabs(diffCounter);
  // cout<<"diff BPC : "<<diffBPC<<"  diffFDC1 : "<<diffFDC1<<"  diffCounter : "<<diffCounter<<endl;
  // cout<<"Factor : "<<f<<endl;
}

bool ProtonArm::fit(const double &mass, TVector3 &vtx, TVector3 &mom, TVector3 &FDC1pos, HodoscopeLikeHit *fc_hit)
{
  VERTEX_POS=vtx;
  MASS=mass;
  FDC1_POS=FDC1pos;
  COUNTER_HIT=fc_hit;

  minuit-> SetPrintLevel(-1);
  minuit-> SetFCN(fcnFC);
  TString par_name[3] = { "mom_{x}", "mom_{y}", "mom_{z}" };
  double init_par[3] = { mom.X(), mom.Y(), mom.Z() };
  double step[3] = { 0.001, 0.001, 0.005 };
  double min[3] = { -0.1, -0.1, 0.5 };
  double max[3] = { 0.1,  0.1, 2.0 };

  int flag=0;
  for( int i=0; i<3; i++ ) minuit->mnparm(i, par_name[i], init_par[i], step[i], min[i], max[i], flag);

  double arglist[20]; arglist[0]=3000, arglist[1]=1;
  int ierflag=0;
  //  int time0=time(0);
  minuit->Command("SET STRategy 0");
  //  cout<<"SIMPLEX"<<endl;
  minuit->mnexcm("SIMPLEX", arglist, 1, ierflag);
  //  cout<<"MIGRAD"<<endl;
  minuit->mnexcm("MIGRAD", arglist, 1, ierflag);
  // int time2=time(0);
  // cout<<"MIGRAD fnish    Time : "<<time2-time0<<endl;

  double par[3];
  double err[3];
  for( int i=0; i<minuit->GetNumPars(); i++ ){
    TString name; double bnd1, bnd2; int flag;
    minuit-> mnpout(i, name, par[i], err[i], bnd1, bnd2, flag);
  }
  //  vtx.SetXYZ(par[0], par[1], vtx.Z());
  mom.SetXYZ(par[0], par[1], par[2]);

  return true;
}

bool ProtonArm::fine_fit(const double &mass, TVector3 &vtx, TVector3 &mom, TVector3 &FDC1pos, HodoscopeLikeHit *fc_hit)
{
  VERTEX_POS=vtx;
  MASS=mass;
  FDC1_POS=FDC1pos;
  COUNTER_HIT=fc_hit;

  minuit-> SetPrintLevel(-1);
  minuit-> SetFCN(fcnFC);
  TString par_name[3] = {"mom_{x}", "mom_{y}", "mom_{z}" };
  double init_par[3] = { mom.X(), mom.Y(), mom.Z() };
  double step[3] = {0.01, 0.01, 0.05 };
  double min[3] = { -0.1, -0.1, mom.Z()-0.1 };
  double max[3] = { 0.1,  0.1, mom.Z()+0.1 };

  int flag=0;
  for( int i=0; i<3; i++ ) minuit->mnparm(i, par_name[i], init_par[i], step[i], min[i], max[i], flag);

  double arglist[20]; arglist[0]=3000, arglist[1]=1;
  int ierflag=0;
  minuit->mnexcm("SIMPLEX", arglist, 1, ierflag);
  minuit->mnexcm("MIGRAD", arglist, 1, ierflag);

  double par[5];
  double err[5];
  for( int i=0; i<minuit->GetNumPars(); i++ ){
    TString name; double bnd1, bnd2; int flag;
    minuit-> mnpout(i, name, par[i], err[i], bnd1, bnd2, flag);
  }
  vtx.SetXYZ(par[0], par[1], vtx.Z());
  mom.SetXYZ(par[2], par[3], par[4]);

  return true;
}

bool ProtonArm::fitBT(const double &mass, const double &mom, const TVector3 &BPCpos, const TVector3 &FDC1pos, ConfMan *conf)
{
  VERTEX_POS=BPCpos;
  MASS=mass;
  FDC1_POS=FDC1pos;

  TVector3 momentum=FDC1pos-BPCpos;
  TVector3 T0gpos;
  conf-> GetGeomMapManager()-> GetGPos(CID_T0, 1, T0gpos);
  START_Z=T0gpos.Z()-0.5;
  momentum.SetMag(mom);
  TVector3 T0pos=BPCpos+((START_Z-BPCpos.Z())/momentum.Z())*momentum;

  minuit-> SetPrintLevel(-1);
  minuit-> SetFCN(fcnBT);
  TString par_name[5] = { "pos_{x}", "pos_{y}", "mom_{x}", "mom_{y}", "mom_{z}" };
  double init_par[5] = { T0pos.X(), T0pos.Y(), momentum.X(), momentum.Y(), momentum.Z() };
  double step[5] = { 0.5, 0.5, 0.01, 0.01, 0.01 };
  double min[5] = { T0pos.X()-5.0, T0pos.Y()-5.0, momentum.X()-0.1, momentum.Y()-0.1, momentum.Z()-0.25 };
  double max[5] = { T0pos.X()+5.0, T0pos.Y()+5.0, momentum.X()+0.1, momentum.Y()+0.1, momentum.Z()+0.25 };

  int flag=0;
  for( int i=0; i<5; i++ ) minuit->mnparm(i, par_name[i], init_par[i], step[i], min[i], max[i], flag);

  double arglist[20]; arglist[0]=3000, arglist[1]=1;
  int ierflag=0;
  minuit->mnexcm("SIMPLEX", arglist, 1, ierflag);
  minuit->mnexcm("MIGRAD", arglist, 1, ierflag);

  return true;
}

ChargeParticle ProtonArm::shootChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom, ConfMan *conf)
{
  // cout<<"===== ProtonArm::shootChargeParticle ====="<<endl;
  // cout<<"Init mom ("<<mom.X()<<", "<<mom.Y()<<", "<<mom.Z()<<")"<<endl;
  ChargeParticle particle(mass, pos, mom);
  if( !checkPreAna() ){
    particle.SetFlag(false);
    return particle;
  }

  if( particle.pos().Z()<ProtonArm::fieldRange_minZ() ){
    //    if( !particle.calcELossBeamTGeo(ProtonArm::fieldRange_minZ()) ){
    if( !particle.calcELossBeamTGeo(ProtonArm::fieldRange_minZ(), ProtonArm::CONFMAN->GetCDSFittingParamManager()->GetMagneticField()) ){
      //      cout<<"  !!! calcELossBeamTGeo return false 0 !!!"<<endl;
      particle.SetFlag(false);
      return particle;
    }
  }
  //  cout<<"field in mom ("<<particle.mom().X()<<", "<<particle.mom().Y()<<", "<<particle.mom().Z()<<")"<<endl;

  //  cout<<"FL before USWK : "<<particle.fl()<<endl;
  double diffFDC1=-1.0;
  while( isUSWK_Field(particle.pos()) ){
    TVector3 pos=particle.pos();
    if( !particle.nextStepRK4(DELTA_DT) ){
      //      cout<<"  !!! nextStep Runge Kutta return false  !!!"<<endl;
      particle.SetFlag(false);
      return particle;
    }
    if( pos.Z()<ProtonArm::FDC1_POS.Z() && ProtonArm::FDC1_POS.Z()<particle.pos().Z() ){
      double rate=(ProtonArm::FDC1_POS.Z()-pos.Z())/(particle.pos().Z()-pos.Z());
      TVector3 FDC1pos=(1.0-rate)*pos+rate*particle.pos();
      diffFDC1=(ProtonArm::FDC1_POS-FDC1pos).Mag();
      particle.SetFDC1pos(FDC1pos);
      // cout<<Form("(%lf, %lf, %lf)", pos.x(), pos.y(), pos.z())<<endl;
      // cout<<Form("(%lf, %lf, %lf)", particle.pos().x(), particle.pos().y(), particle.pos().z())<<endl;
      // cout<<Form("(%lf, %lf, %lf)", FDC1pos.x(), FDC1pos.y(), FDC1pos.z())<<endl;
      // cout<<Form("(%lf, %lf, %lf)", ProtonArm::FDC1_POS.x(), ProtonArm::FDC1_POS.y(), ProtonArm::FDC1_POS.z())<<endl;
    }
  }
  //  cout<<"FL after  USWK : "<<particle.fl()<<endl;
  if( diffFDC1<0. ){
    //    cout<<"  !!! cann't find FDC1 pos !!!"<<endl;
    particle.SetFlag(false);
    return particle;
  }
  //  cout<<"field out mom ("<<particle.mom().X()<<", "<<particle.mom().Y()<<", "<<particle.mom().Z()<<")"<<endl;

  TVector3 counter_hit_pos=crossPointPCCVC(particle.pos(), particle.mom());
  if( !particle.calcELossBeamTGeo(counter_hit_pos.Z()) ){
    particle.SetFlag(false);
    return particle;
  }
  //  particle.dump();
  //  cout<<"FL Counter pos : "<<particle.fl()<<endl;
  particle.SetDiffCounter(ProtonArm::disPCCVC(particle.pos(), particle.mom()));
  //  cout<<"Couter hit mom ("<<particle.mom().X()<<", "<<particle.mom().Y()<<", "<<particle.mom().Z()<<")"<<endl;

  return particle;
}

ChargeParticle ProtonArm::shootChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom, ConfMan *conf, vector<vector<double> > &trajectory)
{
  cout<<"===== ProtonArm::shootChargeParticle ====="<<endl;
  ChargeParticle particle(mass, pos, mom);
  if( !checkPreAna() ){
    particle.SetFlag(false);
    return particle;
  }

  vector<double> point;
  point.push_back(pos.X()), point.push_back(pos.Y()), point.push_back(pos.Z());
  trajectory.push_back(point);
  //  particle.dump();
  if( particle.pos().Z()<fieldRange_minZ() ){
    //    if( !particle.calcELossBeamTGeo(fieldRange_minZ() ) ){ return particle; }
    if( !particle.calcELossBeamTGeo(fieldRange_minZ(), CONFMAN->GetCDSFittingParamManager()->GetMagneticField()) ){ return particle; }
  }
  //  particle.dump();
  point[0]=particle.pos().x(), point[1]=particle.pos().y(), point[2]=particle.pos().z();
  trajectory.push_back(point);

  while( isUSWK_Field(particle.pos()) ){
    if( !particle.nextStepRK4(DELTA_DT) ) return particle;
    point[0]=particle.pos().x(), point[1]=particle.pos().y(), point[2]=particle.pos().z();
    trajectory.push_back(point);
  }
  //  particle.dump();

  TVector3 counter_hit_pos=crossPointPCCVC(particle.pos(), particle.mom());
  if( !particle.calcELossBeamTGeo(counter_hit_pos.Z()) ) return particle;
  //  particle.dump();
  point[0]=particle.pos().x(), point[1]=particle.pos().y(), point[2]=particle.pos().z();
  trajectory.push_back(point);

  return particle;
}

ChargeParticle ProtonArm::shootChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom, ConfMan *conf, 
					      HodoscopeLikeHit *hit, vector<vector<double> > &trajectory)
{
  cout<<"===== ProtonArm::shootChargeParticle ====="<<endl;
  COUNTER_HIT=hit;

  ChargeParticle particle(mass, pos, mom);
  if( !checkPreAna() ){
    particle.SetFlag(false);
    return particle;
  }

  vector<double> point;
  point.push_back(pos.X()), point.push_back(pos.Y()), point.push_back(pos.Z());
  trajectory.push_back(point);
  //  particle.dump();
  if( particle.pos().Z()<fieldRange_minZ() ){
    //    if( !particle.calcELossBeamTGeo(fieldRange_minZ()) ) return particle;
    if( !particle.calcELossBeamTGeo(fieldRange_minZ(), CONFMAN->GetCDSFittingParamManager()->GetMagneticField()) ){ return particle; }
  }
  //  particle.dump();
  point[0]=particle.pos().x(), point[1]=particle.pos().y(), point[2]=particle.pos().z();
  trajectory.push_back(point);

  while( isUSWK_Field(particle.pos()) ){
    if( !particle.nextStepRK4(DELTA_DT) ) return particle;
    point[0]=particle.pos().x(), point[1]=particle.pos().y(), point[2]=particle.pos().z();
    trajectory.push_back(point);
  }
  //  particle.dump();

  TVector3 counter_hit_pos=crossPointPCCVC(particle.pos(), particle.mom());
  if( !particle.calcELossBeamTGeo(counter_hit_pos.Z()) ) return particle;
  //  particle.dump();
  point[0]=particle.pos().x(), point[1]=particle.pos().y(), point[2]=particle.pos().z();
  trajectory.push_back(point);

  return particle;
}

TVector3 ProtonArm::crossPointPCCVC(const TVector3 &pos, const TVector3 &dir)
{
  TVector3 cross_pos;
  TVector3 gpos, grot;
  CONFMAN-> GetGeomMapManager()-> GetPos(COUNTER_HIT->cid(), 0, gpos);
  CONFMAN-> GetGeomMapManager()-> GetRot(COUNTER_HIT->cid(), 0, grot);

  TVector3 cpos, crot;
  double param[20];
  CONFMAN-> GetGeomMapManager()-> GetPos(COUNTER_HIT->cid(), COUNTER_HIT->seg(), cpos);
  CONFMAN-> GetGeomMapManager()-> GetRot(COUNTER_HIT->cid(), COUNTER_HIT->seg(), crot);
  CONFMAN-> GetGeomMapManager()-> GetParam(COUNTER_HIT->cid(), COUNTER_HIT->seg(), param);
  double dz=param[8];

  TVector3 cpos2(cpos.X(), cpos.Y(), cpos.Z()-0.5*dz);
  TVector3 cdir(1.0, 0.0, 0.0);
  cpos2.RotateY(crot.Y()*Deg2Rad);
  cdir.RotateY(crot.Y()*Deg2Rad);
  cpos2.RotateY(grot.Y()*Deg2Rad);
  cdir.RotateY(grot.Y()*Deg2Rad);
  cpos2+=gpos;

  // double dist;
  // TVector3 vtx;
  //  MathTools::LineToLine(pos, dir, cpos2, cdir, 0.0, dist, cross_pos, vtx);
  double A00=cpos2.X();
  double A10=cpos2.Z();
  double A01=cdir.X();
  double A11=cdir.Z();

  double B00=pos.X();
  double B10=pos.Z();
  double B01=dir.X();
  double B11=dir.Z();

  double s=(A00/A01-A10/A11-B00/A01+B10/A11)/(B01/A01-B11/A11);
  double t=(B00+s*B01-A00)/A01;

  if( A11==0 ){
    s=(A10-B10)/B11;
    t=(B00-A00+s*B01)/A01;
  }
  cross_pos=pos+s*dir;
  TVector3 pos2=cpos2+t*cdir;

  // cout<<"A0 : "<<A00<<"  "<<A01<<endl;
  // cout<<"A1 : "<<A10<<"  "<<A11<<endl;
  // cout<<"B0 : "<<B00<<"  "<<B01<<endl;
  // cout<<"B1 : "<<B10<<"  "<<B11<<endl;
  // cout<<Form("Cross Positon  : (%5.3lf, %5.3lf, %5.3lf)", cross_pos.X(), cross_pos.Y(), cross_pos.Z())<<endl;
  // cout<<Form("Cross Positon2 : (%5.3lf, %5.3lf, %5.3lf)", pos2.X(), pos2.Y(), pos2.Z())<<endl;
  return cross_pos;
}

double ProtonArm::momByAng(const double &angle)
{
  double radius=100./sqrt(2*(1+TMath::Cos(angle)));
  return radius*0.9782*Const/100.+0.023;
}

double ProtonArm::disPCCVC(const TVector3 &pos, const TVector3 &dir)
{
  TVector3 cross_pos;
  TVector3 gpos, grot;
  CONFMAN-> GetGeomMapManager()-> GetPos(COUNTER_HIT->cid(), 0, gpos);
  CONFMAN-> GetGeomMapManager()-> GetRot(COUNTER_HIT->cid(), 0, grot);

  TVector3 cpos, crot;
  double param[20];
  CONFMAN-> GetGeomMapManager()-> GetPos(COUNTER_HIT->cid(), COUNTER_HIT->seg(), cpos);
  CONFMAN-> GetGeomMapManager()-> GetRot(COUNTER_HIT->cid(), COUNTER_HIT->seg(), crot);
  CONFMAN-> GetGeomMapManager()-> GetParam(COUNTER_HIT->cid(), COUNTER_HIT->seg(), param);
  double dz=param[8];

  TVector3 cpos2(cpos.X(), cpos.Y(), cpos.Z()-0.5*dz);
  TVector3 cdir(1.0, 0.0, 0.0);
  cpos2.RotateY(crot.Y()*Deg2Rad);
  cdir.RotateY(crot.Y()*Deg2Rad);
  cpos2.RotateY(grot.Y()*Deg2Rad);
  cdir.RotateY(grot.Y()*Deg2Rad);
  cpos2+=gpos;

  // double dist;
  // TVector3 vtx;
  // MathTools::LineToLine(pos, dir, cpos2, cdir, 0.0, dist, cross_pos, vtx);
  //  cout<<cdir.X()<<", "<<cdir.Y()<<", "<<cdir.Z()<<endl;
  // cout<<cross_pos.X()<<", "<<cross_pos.Y()<<", "<<cross_pos.Z()<<endl;
  // cout<<vtx.X()<<", "<<vtx.Y()<<", "<<vtx.Z()<<endl;
  // cout<<cpos2.X()<<", "<<cpos2.Y()<<", "<<cpos2.Z()<<endl;
  double A00=cpos2.X();
  double A10=cpos2.Z();
  double A01=cdir.X();
  double A11=cdir.Z();

  double B00=pos.X();
  double B10=pos.Z();
  double B01=dir.X();
  double B11=dir.Z();

  double s=(A00/A01-A10/A11-B00/A01+B10/A11)/(B01/A01-B11/A11);
  double t=(B00+s*B01-A00)/A01;

  // cout<<"A0 : "<<A00<<"  "<<A01<<endl;
  // cout<<"A1 : "<<A10<<"  "<<A11<<endl;
  // cout<<"B0 : "<<B00<<"  "<<B01<<endl;
  // cout<<"B1 : "<<B10<<"  "<<B11<<endl;
  if( A11==0 ){
    s=(A10-B10)/B11;
    t=(B00-A00+s*B01)/A01;
  }
  cross_pos=pos+s*dir;
  TVector3 pos2=cpos2+t*cdir;

  double disX=cross_pos.X()-cpos2.X();
  double disZ=cross_pos.Z()-cpos2.Z();
  double dis=sqrt(disX*disX+disZ*disZ);

  // cout<<"s : "<<s<<endl;
  // cout<<"t : "<<t<<endl;

  // cout<<Form("Cross Positon  : (%5.3lf, %5.3lf, %5.3lf)", cross_pos.X(), cross_pos.Y(), cross_pos.Z())<<endl;
  // cout<<Form("Cross Positon2 : (%5.3lf, %5.3lf, %5.3lf)", pos2.X(), pos2.Y(), pos2.Z())<<endl;

  if( disX>0 ) return dis;
  else return -dis;
}

TVector3 ProtonArm::getField(const TVector3 &pos)
{
  double x=pos.x(), y=pos.y(), z=pos.z();
  TVector3 field(0.0, 0.0, 0.0);

  double int_x, int_y, int_z;
  double rate_x = modf((x-USWK_GPOS.X()+USWK_XMAX)/USWK_DX, &int_x);
  double rate_y = modf((y-USWK_GPOS.Y()+USWK_YMAX)/USWK_DY, &int_y);
  double rate_z = modf((z-USWK_GPOS.Z()+USWK_ZMAX)/USWK_DZ, &int_z);
  int ix = (int)int_x;
  int iy = (int)int_y;
  int iz = (int)int_z;

  if( ix<0 || USWK_NX<=ix ) return field;
  if( iy<0 || USWK_NY<=iy ) return field;
  if( iz<0 || USWK_NZ<=iz ) return field;

  unsigned long key1 = ix   + USWK_NX*iy     + USWK_NX*USWK_NY*iz;
  unsigned long key2 = ix+1 + USWK_NX*iy     + USWK_NX*USWK_NY*iz;
  unsigned long key3 = ix+1 + USWK_NX*(iy+1) + USWK_NX*USWK_NY*iz;
  unsigned long key4 = ix   + USWK_NX*(iy+1) + USWK_NX*USWK_NY*iz;
  unsigned long key5 = ix   + USWK_NX*iy     + USWK_NX*USWK_NY*(iz+1);
  unsigned long key6 = ix+1 + USWK_NX*iy     + USWK_NX*USWK_NY*(iz+1);
  unsigned long key7 = ix+1 + USWK_NX*(iy+1) + USWK_NX*USWK_NY*(iz+1);
  unsigned long key8 = ix   + USWK_NX*(iy+1) + USWK_NX*USWK_NY*(iz+1);

  TVector3 field1 = USWK_FIELD_MAP[key1];
  TVector3 field2 = USWK_FIELD_MAP[key2];
  TVector3 field3 = USWK_FIELD_MAP[key3];
  TVector3 field4 = USWK_FIELD_MAP[key4];
  TVector3 field5 = USWK_FIELD_MAP[key5];
  TVector3 field6 = USWK_FIELD_MAP[key6];
  TVector3 field7 = USWK_FIELD_MAP[key7];
  TVector3 field8 = USWK_FIELD_MAP[key8];

  TVector3 ave_x1 = rate_x*field1+(1-rate_x)*field2;
  TVector3 ave_x2 = rate_x*field4+(1-rate_x)*field3;
  TVector3 ave_x3 = rate_x*field5+(1-rate_x)*field6;
  TVector3 ave_x4 = rate_x*field8+(1-rate_x)*field7;

  TVector3 ave_y1 = rate_y*ave_x1+(1-rate_y)*ave_x2;
  TVector3 ave_y2 = rate_y*ave_x3+(1-rate_y)*ave_x4;

  field = rate_z*ave_y1+(1-rate_z)*ave_y2;

  return 0.0001*SCALE*field;
}

void ProtonArm::SetCounterHit(HodoscopeLikeHit *hit){ COUNTER_HIT=hit; }

void ProtonArm::Initialize(ConfMan *conf, double val)
{
  CONFMAN=conf;

  minuit = new TMinuit(3);
  TROOT minexam("ProtonArm","ProtonArm fit using TMinuit");

  SCALE=val;
  cout<<" USWK map scale: "<<SCALE<<endl;
  //  SCALE=0.8;
  double scale;

  double x, y, z;
  double Bx, By, Bz;
  int nx, ny, nz;
  FILE *fp;
  char str[1024];

  if( (fp=fopen("param/FieldMap/UshiwakaMap.param","r"))==0 ){
    cout<<"  !!! Ushiwaka Field Map file open error !!!"<<endl;
    exit(0);
  }
  cout << "param/FieldMap/UshiwakaMap.param Opened" << endl;
  while( fgets(str,1024,fp)!=0 ){
    if( str[0]=='#' ) continue;

    if( sscanf(str, "GPOS: %lf %lf %lf",&x,&y,&z)==3 ){
      USWK_GPOS.SetXYZ(x, y, z);
    }
    if( sscanf(str, "MAX: %lf %lf %lf",&x,&y,&z)==3 ){
      USWK_XMAX = x, USWK_YMAX = y, USWK_ZMAX = z;
    }
    if( sscanf(str, "SCALE: %lf",&scale)==3 ){
      SCALE=scale;
    }
    if( sscanf(str, "dl: %lf %lf %lf",&x,&y,&z)==3 ){
      USWK_DX = x, USWK_DY = y, USWK_DZ = z;
    }
    if( sscanf(str, "n: %d %d %d",&nx,&ny,&nz)==3 ){
      USWK_NX = nx, USWK_NY = ny, USWK_NZ = nz;
    }
    if( sscanf(str, " %lf %lf %lf %lf %lf %lf",&x, &y, &z, &Bx, &By, &Bz)==6 ){
      TVector3 field(Bx, By, Bz);
      int index_x = (int)((x+USWK_XMAX)/USWK_DX);
      int index_y = (int)((y+USWK_YMAX)/USWK_DY);
      int index_z = (int)((z+USWK_ZMAX)/USWK_DZ);
      unsigned long key = index_x + USWK_NX*index_y + USWK_NX*USWK_NY*index_z;
      USWK_FIELD_MAP[key] = field;
    }
  }
  //  cout<<"NX*NY*NZ="<<USWK_NX*USWK_NY*USWK_NZ<<"   field map size="<<USWK_FIELD_MAP.size()<<endl;
}

bool ProtonArm::isUSWK_Field(const TVector3 &pos)
{
  if( fabs(pos.X()-USWK_GPOS.X())>USWK_XMAX+0.001 ) return false;
  if( fabs(pos.Y()-USWK_GPOS.Y())>USWK_YMAX+0.001 ) return false;
  if( fabs(pos.Z()-USWK_GPOS.Z())>USWK_ZMAX+0.001 ) return false;

  return true;
}

bool ProtonArm::checkPreAna()
{
  if( !COUNTER_HIT ){
    cout<<"  !!! ProtonArm::COUNTER_HIT is NULL !!!"<<endl;
    return false;
  }
  if( COUNTER_HIT->cid()!=CID_CVC && COUNTER_HIT->cid()!=CID_PC ){
    cout<<"  !!! ProtonArm::COUNTER_HIT isn't CVC/PC !!!"<<endl;
    return false;
  }
  if( !CONFMAN ){
    cout<<"  !!! ProtonArm::CONFMAN is NULL !!!"<<endl;
    return false;
  }
  return true;
}
