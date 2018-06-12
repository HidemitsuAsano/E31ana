#include "SimTools.h"

using namespace std;

DetectorHit *SimTools::GetHit(const int &cid, const int &seg, DetectorData* data)
{
  for(  int i=0; i<data->detectorHitSize(); i++ ){
    DetectorHit *hit=data->detectorHit(i);
    if( hit->detectorID()==cid && hit->channelID()+1==seg ) return hit;
  }
  return 0;
}

DetectorHit *SimTools::GetHit(const int &cid, const int &layer, const int &wire, DetectorData* data)
{
  for(  int i=0; i<data->detectorHitSize(); i++ ){
    DetectorHit *hit=data->detectorHit(i);
    if( hit->detectorID()==cid && hit->layerID()+1==layer && hit->channelID()+1==wire ) return hit;
  }
  return 0;
}

Track *SimTools::GetTrack(const int &id, MCData *data)
{
  for( int i=0; i<data->trackSize(); i++ ){
    if( data->track(i)->trackID()==id ) return data->track(i);
  }
  return 0;
}

int SimTools::Generation(const int &id, MCData *data)
{
  Track *track=GetTrack(id, data);
  int generation=0;
  while(track){
    track=GetTrack(track->parentTrackID(), data);
    generation++;
  }
  return generation;
}

vector<Track*> SimTools::GetDaughterTracks(const int &id, MCData *data)
{
  vector<Track*> daughterTracks;
  for( int i=0; i<data->trackSize(); i++ ){
    if( data->track(i)->parentTrackID()==id ) daughterTracks.push_back(data->track(i));
  }
  return daughterTracks;
}

vector<Track*> SimTools::GetParentTracks(const int &id, MCData *data)
{
  vector<Track*> parentTracks;
  Track *track=GetTrack(id, data);

  while(true){
    track=GetTrack(track->parentTrackID(), data);
    if( track ) parentTracks.push_back(track);
    else break;
  }
  return parentTracks;
}

string SimTools::ReactionName(const int &processID)
{
  if( processID==0        ) return "Initial";
  else if( processID==1   ) return "Decay";
  else if( processID==2   ) return "conv";
  else if( processID==3   ) return "Transportation";
  else if( processID==4   ) return "phot";
  else if( processID==5   ) return "annihil";
  else if( processID==6   ) return "compt";
  else if( processID==7   ) return "eBrem";
  else if( processID==8   ) return "hadElastic";
  else if( processID==9   ) return "CoulombScat";
  else if( processID==10  ) return "nKiller";
  else if( processID==11  ) return "photoNuclear";
  else if( processID==100 ) return "pi-Inelastic";
  else if( processID==101 ) return "pi+Inelastic";
  else if( processID==102 ) return "kaon-Inelastic";
  else if( processID==103 ) return "kaon+Inelastic";
  else if( processID==104 ) return "kaon0LInelastic";
  else if( processID==105 ) return "kaon0SInelastic";
  else if( processID==106 ) return "lambdaInelastic";
  else if( processID==107 ) return "sigma+Inelastic";
  else if( processID==108 ) return "sigma-Inelastic";
  else if( processID==109 ) return "sigma0Inelastic";
  else if( processID==110 ) return "protonInelastic";
  else if( processID==111 ) return "neutronInelastic";
  else if( processID==112 ) return "dInelastic";
  else if( processID==113 ) return "tInelastic";
  else if( processID==200 ) return "eIoni";
  else if( processID==201 ) return "hIoni";
  else if( processID==202 ) return "ionIoni";
  else if( processID==203 ) return "muIoni";
  else if( processID==204 ) return "hBertiniCAptureAtRest";
  else if( processID==205 ) return "nCaptrue";
  else if( processID==206 ) return "muMinusCaptrueAtRest";
  else return "unknown";
}

string SimTools::ParticleName(const int &pdg)
{
  if( pdg==11 ) return "e-";
  else if( pdg==22    ) return "gamma";
  else if( pdg==12    ) return "nu_e";
  else if( pdg==13    ) return "mu-";
  else if( pdg==14    ) return "nu_mu";
  else if( pdg==15    ) return "tau-";
  else if( pdg==16    ) return "nu_tau";
  else if( pdg==-11   ) return "e+";
  else if( pdg==-12   ) return "anti_nu_e";
  else if( pdg==-13   ) return "mu+";
  else if( pdg==-14   ) return "anti_nu_mu";
  else if( pdg==-15   ) return "tau+";
  else if( pdg==-16   ) return "anti_nu_tau";
  else if( pdg==323   ) return "k_star+";
  else if( pdg==313   ) return "k_star0";
  else if( pdg==321   ) return "kaon+";
  else if( pdg==311   ) return "kaon0";
  else if( pdg==-323  ) return "k_star-";
  else if( pdg==-313  ) return "anti_k_star0";
  else if( pdg==-321  ) return "kaon-";
  else if( pdg==-311  ) return "anti_kaon0";
  else if( pdg==130   ) return "kaon0L";
  else if( pdg==310   ) return "kaon0S";
  else if( pdg==211   ) return "pi+";
  else if( pdg==-211  ) return "pi-";
  else if( pdg==111   ) return "pi0";
  else if( pdg==113   ) return "rho0";
  else if( pdg==213   ) return "rho+";
  else if( pdg==-213  ) return "rho-";
  else if( pdg==2214  ) return "delta+";
  else if( pdg==1114  ) return "delta-";
  else if( pdg==2114  ) return "delta0";
  else if( pdg==3122  ) return "lambda";
  else if( pdg==13122 ) return "lambda(1405)";
  else if( pdg==3124  ) return "lambda(1520)";
  else if( pdg==23122 ) return "lambda(1600)";
  else if( pdg==33122 ) return "lambda(1670)";
  else if( pdg==13124 ) return "lambda(1690)";
  else if( pdg==43122 ) return "lambda(1800)";
  else if( pdg==53122 ) return "lambda(1810)";
  else if( pdg==3126  ) return "lambda(1820)";
  else if( pdg==13126 ) return "lambda(1830)";
  else if( pdg==23124 ) return "lambda(1890)";
  else if( pdg==2112  ) return "neutron";
  else if( pdg==2212  ) return "proton";
  else if( pdg==3224  ) return "sigma(1385)+";
  else if( pdg==3114  ) return "sigma(1385)-";
  else if( pdg==3214  ) return "sigma(1385)0";
  else if( pdg==3222  ) return "sigma+";
  else if( pdg==3112  ) return "sigma-";
  else if( pdg==3212  ) return "sigma0";
  else if( pdg==1000010020 ) return "deuteron";
  else if( pdg==1000010030 ) return "triton";
  else if( pdg==1000020030 ) return "He3";
  else if( pdg==1000020040 ) return "alpha";
  int code1=pdg%1000000;
  code1/=10;
  int code2=pdg%1000;
  code2/=10;

  if( code1>(int)(sizeof(ElementSymbol)/sizeof(ElementSymbol[0])) ) return Form("%d", pdg);
  return Form("%s%d", ElementSymbol[code1].Data(), code2);

}
