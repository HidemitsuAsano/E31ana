#include "ChargeParticle.h"

using namespace std;

ChargeParticle::ChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom)
  : fFlag(true), fMass(mass), fInitMom(mom.Mag()), fFDC1pos(DEFVECT), fDiffCounter(DEFAULTD),
    fDeltaTime(0), fFlightLength(0), fTime(0), fPos(pos), fMom(mom)
{
}

bool ChargeParticle::calcELossBeamTGeo(const double &z)
{
  if( (z-fPos.Z())*fMom.Z()<0 ){
    fFlag=false;
    return false;
  }
  TVector3 next_pos=fPos+((z-fPos.Z())/fMom.Z())*fMom;

  if( fMass!=0 ){
    double mom_out, tof;
    if( !ELossTools::CalcElossBeamTGeo(fPos, next_pos, fMom.Mag(), fMass, mom_out, tof) ){
      fFlag=false;
      return false;
    }
    fFlightLength += (next_pos-fPos).Mag();
    fDeltaTime+=tof-tofByInitMom((next_pos-fPos).Mag());
    fTime += tof;
    fPos=next_pos;
    fMom.SetMag(mom_out);
  }
  else{
    fFlightLength += (next_pos-fPos).Mag();
    fTime += fFlightLength/(100.*Const);
    fPos=next_pos;
  }
  
  return true;
}

bool ChargeParticle::calcELossBeamTGeo(const double &z, const double &cds_field)
{
  if( (z-fPos.Z())*fMom.Z()<0 ){
    fFlag=false;
    return false;
  }

  if( fPos.Z()<50 ){
    TVector3 pos_out, mom_out;
    double tmp_fl, tmp_tof;
    if( !ELossTools::CalcHelixElossPointToPointTGeo(fPos, fMom, fMass, 50, cds_field, pos_out, mom_out, tmp_fl, tmp_tof) ){
      cout<<"!!!!! CalcHElixElossPointToPointTGeo() return false !!!!!"<<endl;
      fFlag=false;
      return false;
    }
    fPos=pos_out;
    fMom=mom_out;
    fTime += tmp_tof;
    fFlightLength+=tmp_fl;

    if( fMass>0 ) fDeltaTime+=tmp_tof-tofByInitMom(tmp_fl);
  }

  TVector3 next_pos=fPos+((z-fPos.Z())/fMom.Z())*fMom;

  if( fMass!=0 ){
    double mom_out, tof;
    if( !ELossTools::CalcElossBeamTGeo(fPos, next_pos, fMom.Mag(), fMass, mom_out, tof) ){
      fFlag=false;
      return false;
    }
    fFlightLength += (next_pos-fPos).Mag();
    fDeltaTime+=tof-tofByInitMom((next_pos-fPos).Mag());
    fTime += tof;
    fPos=next_pos;
    fMom.SetMag(mom_out);
  }
  else{
    fFlightLength += (next_pos-fPos).Mag();
    fTime += (next_pos-fPos).Mag()/(100.*Const);
    fPos=next_pos;
  }
  
  return true;
}

bool ChargeParticle::nextStepRK4(const double &dt)
{
  TVector3 v0=(100.*Const/sqrt(fMom.Mag2()+fMass*fMass))*fMom;
  TVector3 kmom1=0.01*Const*v0.Cross(ProtonArm::getField(fPos));
  TVector3 kpos1=v0;
  TVector3 mom1=fMom+0.5*kmom1*dt;
  TVector3 pos1=fPos+0.5*kpos1*dt;
  double mom_out, tof;
  //  cout<<"Pos : "<<Form("(%lf, %lf, %lf)   ", fPos.X(), fPos.Y(), fPos.Z())<<ProtonArm::getField(fPos).Y()<<endl;
  //  cout<<MyTools::str(ProtonArm::getField(fPos))<<endl;
  // if( !ELossTools::CalcElossBeamTGeo(fPos, pos1, fMom.Mag(), fMass, mom_out, tof) ) return false;
  // mom1.SetMag(mom_out);

  TVector3 v1=(100.*Const/sqrt(mom1.Mag2()+fMass*fMass))*mom1;
  TVector3 kmom2=0.01*Const*v1.Cross(ProtonArm::getField(pos1));
  TVector3 kpos2=v1;
  TVector3 mom2=fMom+0.5*kmom2*dt;
  TVector3 pos2=fPos+0.5*kpos2*dt;
  //  cout<<MyTools::str(ProtonArm::getField(pos1))<<endl;
  // if( !ELossTools::CalcElossBeamTGeo(fPos, pos2, fMom.Mag(), fMass, mom_out, tof) ) return false;
  // mom2.SetMag(mom_out);

  TVector3 v2=(100.*Const/sqrt(mom2.Mag2()+fMass*fMass))*mom2;
  TVector3 kmom3=0.01*Const*v2.Cross(ProtonArm::getField(pos2));
  TVector3 kpos3=v2;
  TVector3 mom3=fMom+kmom3*dt;
  TVector3 pos3=fPos+kpos3*dt;
  //  cout<<MyTools::str(ProtonArm::getField(pos2))<<endl;
  // if( !ELossTools::CalcElossBeamTGeo(fPos, pos3, fMom.Mag(), fMass, mom_out, tof) ) return false;
  // mom3.SetMag(mom_out);

  TVector3 v3=(100.*Const/sqrt(mom3.Mag2()+fMass*fMass))*mom3;
  TVector3 kmom4=0.01*Const*v3.Cross(ProtonArm::getField(pos3));
  TVector3 kpos4=v3;
  //  cout<<MyTools::str(ProtonArm::getField(pos3))<<endl;

  TVector3 next_pos=fPos+(kpos1+2.*kpos2+2.*kpos3+kpos4)*(dt/6.);
  TVector3 next_mom=fMom+(kmom1+2.*kmom2+2.*kmom3+kmom4)*(dt/6.);

  if( fMass!=0 ){
    if( !ELossTools::CalcElossBeamTGeo(fPos, next_pos, fMom.Mag(), fMass, mom_out, tof) ){
      fFlag=false;
      return false;
    }
    //  cout<<" dmom : "<<fMom.Mag()-mom_out<<endl;
    double tmp_fl=0.5*dt*(v().Mag()+100.*Const*mom_out/sqrt(mom_out*mom_out+fMass*fMass));
    //  double tmp_fl=dt*v().Mag();

    fFlightLength+=tmp_fl;
    //    if( dt<tofByInitMom(tmp_fl) ) cout<<"  !!! dt<calc_tof !!!"<<endl;
    fDeltaTime+=dt-tofByInitMom(tmp_fl);
    fTime += dt;
    fPos=next_pos;
    fMom=next_mom;
    fMom.SetMag(mom_out);
  }
  else{
    fFlightLength+=dt*100.*Const;
    fTime += dt;
    fPos=next_pos;
    fMom=next_mom;
  }

  // TString mat_name;
  // GeomTools::GetIDMat(fPos, mat_name);
  //  if( mat_name!="Air" && mat_name!="BLDCGas" ) cout<<"  Mat : "<<mat_name<<endl;

  return true;
}

void ChargeParticle::dump() const
{
  cout<<"mass   : "<<fMass<<endl;
  cout<<"pos "<<Form("(%6.4lf, %6.4lf, %6.4lf)", fPos.X(), fPos.Y(), fPos.Z())<<" [cm]"<<endl;
  cout<<"mom "<<Form("(%6.4lf, %6.4lf, %6.4lf)", fMom.X(), fMom.Y(), fMom.Z())<<" [GeV/c]"<<endl;
  cout<<"time : "<<fTime<<" [ns]    delta tof : "<<fDeltaTime<<" [ns]"<<endl;
  cout<<"Flight Length : "<<fFlightLength<<" [cm]"<<endl;
}

