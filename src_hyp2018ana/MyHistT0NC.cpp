#include "MyHistT0NC.h"

using namespace std;

void initHistT0NC()
{
  const int nseg=112;
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("T0%d_offset_NC", seg), Form("T0 seg%d offset", seg), 1000, -5.0, 5.0);

    new TH2F(Form("T0%d_eu_offset", seg), Form("T0 seg%d dE_{up} vs offset", seg),   100, 0, 10, 1000, -1, 1);
    new TH2F(Form("T0%d_ed_offset", seg), Form("T0 seg%d dE_{down} vs offset", seg), 100, 0, 10, 1000, -1, 1);
  }

  new TH1F("NC_offset", "NC offset by #gamma-ray", 1000, -5.0, 5.0);
  for( int seg=1; seg<=nseg; seg++ ){
    new TH1F(Form("NC%d_offset", seg), Form("NC seg%d offset", seg), 1000, -5.0, 5.0);
    new TH1F(Form("NC%d_offset2", seg), Form("NC seg%d offset", seg), 1000, -5.0, 5.0);

    new TH2F(Form("NC%d_eu_offset", seg), Form("NC seg%d dE_{up} vs offset", seg),   100, 0, 100, 1000, -1, 1);
    new TH2F(Form("NC%d_ed_offset", seg), Form("NC seg%d dE_{down} vs offset", seg), 100, 0, 100, 1000, -1, 1);

    new TNtuple(Form("NC%d_slewing_info", seg), Form("NC seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");

    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("NC%d_offset_T0%d", seg, seg2), Form("NC seg%d offset by T0%d", seg, seg2), 1000, -5.0, 5.0);
      new TH1F(Form("NC%d_offset2_T0%d", seg, seg2), Form("NC seg%d offset by T0%d", seg, seg2), 1000, -5.0, 5.0);
    }
  }
}

void fillT0NC(BeamLineHitMan *blMan, AnaInfo *info)
{
  if( info->nBeam()!=1 || !info->beam(0)->flag() ) return;
  if( info->nFNeutral()!=1 ) return;
  if( info->forwardNeutral(0)->pid()==F_Gamma ){
    ForwardNeutralInfo *fnInfo=info->forwardNeutral(0);
    BeamInfo *beam = info->beam(0);
    double offset=fnInfo->time()-beam->T0time()-fnInfo->offset();

    // cout<<"T0-NC  NCseg"<<fnInfo->seg()<<"  offset="<<offset<<endl;
    // cout<<Form("T0pos (%lf, %lf, %lf)", beam->T0pos().X(), beam->T0pos().Y(), beam->T0pos().Z())<<endl;
    // cout<<Form("Vtx (%lf, %lf, %lf)", fnInfo->vertex().X(), fnInfo->vertex().Y(), fnInfo->vertex().Z())<<endl;

    HodoscopeLikeHit *hit = fnInfo->NC(blMan);
    HodoscopeLikeHit *T0hit=info->beam(0)->T0(blMan);

    double eu=hit->eu(), ed=hit->ed();
    double tu=hit->tu()-beam->T0time()-fnInfo->offset();
    double td=hit->td()-beam->T0time()-fnInfo->offset();
    double ctu=hit->ctu()-beam->T0time()-fnInfo->offset();
    double ctd=hit->ctd()-beam->T0time()-fnInfo->offset();

    TNtuple *tup=(TNtuple*)gFile->Get(Form("NC%d_slewing_info", hit->seg()));
    tup->Fill(eu, ed, tu, td, ctu, ctd, offset);

    MyHistTools::fillTH(Form("T0%d_eu_offset", T0hit->seg()), T0hit->eu(), -offset);
    MyHistTools::fillTH(Form("T0%d_ed_offset", T0hit->seg()), T0hit->ed(), -offset);

    MyHistTools::fillTH(Form("NC%d_offset", fnInfo->seg()), offset);
    MyHistTools::fillTH(Form("NC%d_offset_T0%d", fnInfo->seg(), info->beam(0)->T0seg()), offset);
    MyHistTools::fillTH(Form("T0%d_offset_NC", info->beam(0)->T0seg()), -offset);

    if( info->minDCA() && GeomTools::GetID(info->minDCA()->vertexBeam())==CID_Fiducial ){
      MyHistTools::fillTH("NC_offset", offset);
    }
  }
}

void fillSlewingNC(BeamLineHitMan *blMan, AnaInfo *info, ForwardNeutralInfo *fnInfo)
{
  if( info->nBeam()!=1 || !info->beam(0)->flag() ) return;
  if( fnInfo->pid()==F_Gamma ){
    BeamInfo *beam = info->beam(0);
    double offset=fnInfo->time()-beam->T0time()-fnInfo->offset();

    HodoscopeLikeHit *hit = fnInfo->NC(blMan);
    double eu=hit->eu(), ed=hit->ed();

    MyHistTools::fillTH(Form("NC%d_eu_offset", fnInfo->seg()), eu, offset);
    MyHistTools::fillTH(Form("NC%d_ed_offset", fnInfo->seg()), ed, offset);

    MyHistTools::fillTH(Form("NC%d_offset2", fnInfo->seg()), offset);
    MyHistTools::fillTH(Form("NC%d_offset2_T0%d", fnInfo->seg(), info->beam(0)->T0seg()), offset);
  }
}
