#include "MyTools.h"

int MyTools::FindMass2(CDSTrack *track, LocalTrack *bpc, const double tof, const double beammom, const int pid_beam, 
		       double &calc_beta, double &mass2, double &tmptof,
		       TVector3 &vtx_beam, TVector3 &vtx_cds, TVector3 &p, double &out_beam, double &beam_tof)
{
  mass2 = piMass*piMass;
  double param[5];
  TVector3 tmp_vec;
  track-> GetParameters(CID_CDC, param, tmp_vec);
  TVector3 T0pos = bpc-> GetPosatZ(-110.5);
  TVector3 pos1 = MathTools::CalcHelixPosatR(param, 15.4);
  TVector3 cdh_vtx = track-> CDHVertex();
  double cdc_dis = MathTools::CalcHelixArc(param, cdh_vtx, pos1);
  double tmp_dis;

  double mom = track-> Momentum();
  int count = 0;
  while( true ){
    if( mass2<0 ) mass2=0;
    double tmpmass2=mass2;
    double tmpmass;
    if( tmpmass2<0 ) tmpmass = 0.0005;
    else tmpmass = sqrt(tmpmass2);

    double tmpl;
    if( !track->CalcVertexTimeLength(T0pos, bpc->GetMomDir(), tmpmass, vtx_beam, vtx_cds, tmptof, tmpl, true) ){
      //      std::cout<<"  !!! track->CalcVertexTimeLength   return false !!!  dis="<<tmpl<<std::endl;
      //      std::cout<<"      count : "<<count<<std::endl;
      return 1;
    }
    if( !ELossTools::CalcElossBeamTGeo(T0pos, vtx_beam, beammom, parMass[pid_beam], out_beam, beam_tof) ){
      //      std::cout<<"  !!! TrackTools::CalcElossBeamTGeo   return false !!!"<<std::endl;
      return 2;
    }

    tmptof += beam_tof;
    calc_beta = cdc_dis/((tof-tmptof)*Const*100.);
    mass2 = mom*mom*(1/(calc_beta*calc_beta)-1);
    mass2 = (tmpmass2+mass2)/2.;

    if( mass2<=0 ) break;
    else if( abs(tmpmass2-mass2)<0.0001 ) break;
    
    count++;
    if( count>=50 ){
      //      std::cout<<"  !!! MyTools::FindMass2C   count over 50  return false !!!"<<std::endl;
      return 3;
    }
  }
  if( !track-> GetMomentum(vtx_cds, p, true, true) ){
    //    std::cout<<"  !!! MyTools::FoindMass2C  track->GetMomentum() false !!!"<<std::endl;
    return 4;
  }
  //  std::cout<<"tmptof,beta,cdcdis,mass2  "<<tmptof<<"  "<<calc_beta<<"  "<<cdc_dis<<"  "<<mass2<<std::endl;
  return 11;
}

bool MyTools::IsTarget(const TVector3 &pos)
{
  if( -8<pos.Z() && pos.Z()<3 ){
    double r = sqrt(pos.X()*pos.X()+pos.Y()*pos.Y());
    if( r<3 ) return true;
  }
  return false;
}
