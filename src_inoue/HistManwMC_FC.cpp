#include "HistManwMC.h"
#include "MyParam.h"
#include "ProtonArm.h"

using namespace std;
static const double FC_P_MIN=0.5;
static const double FC_P_MAX=1.5;


void HistManwMC::anaFC(EventHeader *header)
{
  if( simReader ) return;

  TH1F *h1;
  TH2F *h2;
  //  cout<<"===== HistManwMC::anaFC ====="<<endl;
  if( bltrackMan->ntrackFDC1()==1 ) fFDC1track=bltrackMan->trackFDC1(0);
  if( fCVC_hit.size()==1 && fPC_hit.size()==0 ) fFChit=fCVC_hit[0];
  if( fCVC_hit.size()==0 && fPC_hit.size()==1 ) fFChit=fPC_hit[0];
  fFC_start=fVtxBeam;

  if( fFDC1track ) fFC_FDC1pos=fFDC1track->GetPosatZ(fFDC1track->hit(0)->gz());
  if( fFChit ) confMan->GetGeomMapManager()->GetGPos(fFChit->cid(), fFChit->seg(), fFC_hitpos);

  // if( header ){
  //   if( header->IsTrig(Trig_Charged) ){
  //     cout<<"> Charge Trigger"<<endl;
  //     cout<<"> nFDC1 : "<<bltrackMan->ntrackFDC1()<<endl;
  //     cout<<"> nCVC  : "<<fCVC_hit.size()<<endl;
  //     cout<<"> nPC   : "<<fPC_hit.size()<<endl;
  //     cout<<"> nBVC  : "<<fBVC_hit.size()<<endl;
  //   }
  //   // else{
  //   //   cout<<"> Not Charge Trigger"<<endl;
  //   //   cout<<"> nFDC1 : "<<bltrackMan->ntrackFDC1()<<endl;
  //   //   cout<<"> nCVC  : "<<fCVC_hit.size()<<endl;
  //   //   cout<<"> nPC   : "<<fPC_hit.size()<<endl;
  //   //   cout<<"> nBVC  : "<<fBVC_hit.size()<<endl;
  //   // }
  // }

  if( fFDC1track && fFChit ){
    TVector3 indir=fFC_FDC1pos-fFC_start;
    TVector3 USWK_pos=fFC_FDC1pos+((250.-fFC_start.Z())/indir.Z())*indir;
    TVector3 outdir=fFC_hitpos-USWK_pos;
    fFC_Angle=indir.Angle(outdir);
    h1 = (TH1F*)rtFile->Get("FC_angle"), h1-> Fill(fFC_Angle);

    TVector3 mom(0, 0, 1.0);
    TVector3 vtx=fFC_start;
    //    cout<<" FPID : "<<fFPID<<endl;
    if( !ProtonArm::fit(particleMass[fFPID], vtx, mom, fFC_FDC1pos, fFChit) ){
      cout<<" !!! ProtonArm::fit return false"<<endl;
    }
    ChargeParticle charged=ProtonArm::shootChargeParticle(particleMass[fFPID], vtx, mom, confMan);
    fFCflag=true;

    // cout<<"> Vtx0    : "<<Form("(%lf, %lf, %lf)", fVtxBeam.x(), fVtxBeam.y(), fVtxBeam.z())<<endl;
    // cout<<"> Vtx     : "<<Form("(%lf, %lf, %lf)", fFC_start.x(), fFC_start.y(), fFC_start.z())<<endl;
    // cout<<"> Fit Vtx : "<<Form("(%lf, %lf, %lf)", vtx.x(), vtx.y(), vtx.z())<<endl;
    // cout<<"> Mom by USWK : "<<mom.Mag()<<"[GeV/c] "<<Form("(%lf, %lf, %lf)", mom.x(), mom.y(), mom.z())<<endl;

    fFC_Mom_USWK=mom.Mag();
    double fl=charged.fl();

    double beam_out, tmp_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, fFC_start, fD5mom, parMass[fBeamPID], beam_out, tmp_tof);
    double tof=fFChit->ctmean()-fT0time-tmp_tof;
    fFC_beta=fl/(100.*Const*tof);
    fFC_mass2=mom.Mag2()*(1./(fFC_beta*fFC_beta)-1.);
    // cout<<">> beta  : "<<fFC_beta<<endl;
    // cout<<">> mass2 : "<<fFC_mass2<<endl;

    TVector3 calc_hitpos=charged.pos();

    h1 = (TH1F*)rtFile->Get("FC_mass2"),     h1-> Fill(fFC_mass2);
    h2 = (TH2F*)rtFile->Get("FC_mass2_mom"), h2-> Fill(fFC_mass2, fFC_Mom_USWK);
    if( FC_P_MIN<fFC_mass2 && fFC_mass2<FC_P_MAX ){
      fFPID=F_Proton;
      vtx=fFC_start;
      if( !ProtonArm::fit(particleMass[fFPID], vtx, mom, fFC_FDC1pos, fFChit) ){
	cout<<" !!! ProtonArm::fit2 return false"<<endl;
      }
      ChargeParticle proton=ProtonArm::shootChargeParticle(particleMass[fFPID], vtx, mom, confMan);

      // cout<<">>> FC proton : "<<fFC_mass2<<endl;
      // cout<<">>> FL  : "<<fl<<"  "<<proton.fl()<<endl;
      // cout<<">>> Time diff : "<<proton.deltaTime()<<endl;
      // cout<<">>> TOF : "<<tof<<"  "<<proton.time()<<endl;
      // cout<<">>> offset : "<<tof-proton.time()<<endl;
      fFC_Mom_USWK=mom.Mag();

      fl=proton.fl();
      tof -= proton.deltaTime();
      double beta=fl/(tof*100.*Const);
      fFC_Mom_TOF=particleMass[fFPID]*beta/sqrt(1-beta*beta);
      //      cout<<">>> mom : "<<fFC_Mom_USWK<<"  "<<fFC_Mom_TOF<<endl;

      //      string str; cin>>str; if(str=="q") exit(0);
    }
  }
}
