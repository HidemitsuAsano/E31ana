// SimDataMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <new>

#include "SimDataMan.h"

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SimDataMan::SimDataMan()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
  SimDataMan::SimDataMan(TFile *file ,ConfMan *conf, bool tracking_flag)
: outFile(0), tree2(0), tree(0), confMan(0), runHeader(0), eventHeader(0), detectorData(0), mcData(0), reactionData(0)
{
  const std::string funcname = "SimDataMan::Initialize";
  confMan = conf;

  std::cout<<"["<<funcname<<"] Initialization start ...";
  inFile = file;
  if( !inFile->IsOpen() ){
    std::cout<<" !!! TFile("<<inFile->GetName()<<") is not open !!!"<<std::endl;
    exit(-1);
  }

  tree2 = (TTree*)inFile-> Get("tree2");
  if( !tree2 ){
    std::cout<<" !!! not find TTree(tree2) !!!"<<std::endl;
    exit(-1);
  }
  tree2-> SetBranchAddress("RunHeaderMC", &runHeader);
  if( tree2->GetEntries()==1 )  tree2-> GetEntry(0);
  else std::cout<<"tree2 is filled "<<tree2->GetEntries()<<std::endl;

  tree  = (TTree*)inFile-> Get("tree");
  if( !tree ){
    std::cout<<" !!! not find TTree(tree) !!!"<<std::endl;
    exit(-1);
  }
  NumOfEvent = tree-> GetEntries();

  tree-> SetBranchAddress("EventHeaderMC", &eventHeader);
  tree-> SetBranchAddress("DetectorData", &detectorData);
  tree-> SetBranchAddress("MCData", &mcData);
  tree-> SetBranchAddress("ReactionData", &reactionData);

  std::cout<<" finish."<<std::endl;

  Clear();
  //  PrintHeader();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SimDataMan::Get(const int &evn)
{
  if( evn<-1 || NumOfEvent<evn ){
		std::cout<<" !!! SimDataMan::Get("<<evn<<") return false !!!"<<std::endl;
		return false;
	}
	tree-> GetEntry(evn);
	return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SimDataMan::Get(const int &evn, CDSHitMan *cdsMan, BeamLineHitMan *blMan)
{
	if( evn<-1 || NumOfEvent<evn ){
		std::cout<<" !!! SimDataMan::Get("<<evn<<") return false !!!"<<std::endl;
		return false;
	}
	tree-> GetEntry(evn);
	Check();
	bool T0_flag = false;
	for( int i=0; i<detectorData->detectorHitSize(); i++ ){
		DetectorHit *mchit = detectorData-> detectorHit(i); 
		int cid       = mchit->detectorID();
		int layer     = mchit->layerID();
		int channel   = mchit->channelID();
		int pdg_id    = mchit->pdg();
		int parent_id = mchit->parentID();
		double time = mchit->time();
		double dE = mchit->de();
		TVector3 pos = mchit->pos();
		if( pdg_id==321 && parent_id==0 ) time *= -1;

		if( cid==CID_CDH || cid==CID_IH || cid==CID_BHD || cid==CID_DEF || cid==CID_T0 || cid==CID_BPD || 
				cid==CID_LC1 || cid==CID_LC2 || cid==CID_BVC || cid==CID_NC || cid==CID_CVC || cid==CID_PC ){

			if( cid!=CID_NC && (pdg_id==22 || pdg_id==2112) ) continue;
			else if( dE<0.1 ) continue; 

			if( cid==CID_T0 && pdg_id==321 && parent_id==0 ){
				T0Time = time;
			}
			if( cid==CID_CDH ) CDHHitPDG.push_back(pdg_id);
			if( cid==CID_IH  ) nIH++;
			if( cid==CID_DEF ) nDEF++;
			if( cid==CID_T0  ) nT0++;
			if( cid==CID_BPD ) nBPD++;
			if( cid==CID_BVC ) nBVC++;
			if( cid==CID_NC  ) nNC++;
			if( cid==CID_CVC ) nCVC++;
			if( cid==CID_PC  ) nPC++;

			int seg=channel+1;

			HodoscopeLikeHit hit;
			hit.SetCounterID(cid);
			hit.SetSegment(seg);
			if( cid==CID_IH ) hit.SetNSensor(1);
			hit.SetCTMean(time);
			hit.SetEMean(dE);
			double hitpos;
			if( cid==CID_CDH || cid==CID_IH ){ hitpos = pos.z()/10.0;}
			else{ hitpos = pos.Y()/10.0;}
			hit.SetHitPosition(hitpos);

			if( confMan!=0 ){
				int c, n, a;
				for( int at=0; at<2; at++ ){
					for( int ud=0; ud<2; ud++ ){
						confMan-> GetCounterMapManager()-> GetCNA(cid,seg,at,ud,c,n,a);
						hit.SetCrate(at, ud, c);
						hit.SetSlot(at, ud, n);
						hit.SetChannel(at, ud, a);
					}
				}

				double lv;
				confMan-> GetGeomMapManager()-> GetLightVelocity(cid, seg, lv);
				TVector3 gpos;
				confMan-> GetGeomMapManager()-> GetGPos(cid, seg, gpos);

				if( cid==CID_CDH || cid==CID_IH ) gpos.SetZ(hitpos);
				else gpos.SetY(hitpos);

				hit.SetPos(gpos);
				hit.SetCTSub(hitpos/lv);

				hit.SetSimulatedResolution(confMan);
				hit.Reverse(confMan);
				double atime     = hit.ctmean();
				double adE       = hit.emean();
				TVector3 ahitpos = hit.pos();
				//hit.Calc(confMan);
				hit.SetCTMean(atime);
				hit.SetEMean(adE);
				hit.SetPos(ahitpos);
			}
			else std::cout<<" !!! SimDataMan Conffile not find !!!"<<std::endl;
			if( cid==CID_CDH || cid==CID_IH ) cdsMan->AddHit(hit);
			else                              blMan-> AddHit(hit);

		}
		else if( cid==CID_CDC ){
			layer += 1;
			nCDC[layer-1]++;

			int wire = channel+1;
			double dx = fabs(mchit->dx()/10.);
			CDCHit hit;
			hit.SetCounterID(cid);
			hit.SetLayer(layer);
			hit.SetWire(wire);
			hit.SetDriftLength(dx);
			hit.SetTimeOffsetTrue(time);

			if( confMan!=0 ){
				int c,n,a;
				int slayer, asdnum, ttype, asdch;
				double rad, phi, tilt;
				confMan->GetCDCWireMapManager()->GetWire( layer, wire, slayer, asdnum, ttype, asdch, rad, phi, tilt, c, n, a );
				hit.SetCrate(c);
				hit.SetSlot(n);
				hit.SetChannel(a);
				hit.SetSuperLayer(slayer);
				hit.SetASDNum(asdnum);
				hit.SetTransType(ttype);
				hit.SetASDChannel(asdch);
				hit.SetRadius(rad); 
				hit.SetPhi(phi);
				hit.SetTiltAngle(tilt);
				hit.SetWirePosition(rad*cos(phi/180.*3.141592),rad*sin(phi/180.*3.141592),0);
			}
			hit.SetSimulatedResolution(confMan);
			hit.Reverse(confMan);
			double adl = hit.dl();
			hit.Calc(confMan);
			hit.SetDriftLength(adl);

			cdsMan-> AddHit(hit);
		}
		else if( cid==CID_BPC || cid==CID_FDC1 || cid==CID_BLC2a || cid==CID_BLC2b ){
			layer += 1;
			if( cid==CID_BPC )   nBPC[layer-1]++;
			if( cid==CID_BLC2a ) nBLC2a[layer-1]++;
			if( cid==CID_BLC2b ) nBLC2b[layer-1]++;
			if( cid==CID_FDC1  ) nFDC1[layer-1]++;

			int wire = channel+1;
			double dx = fabs(mchit->dx()/10.);
			ChamberLikeHit hit;
			hit.SetCounterID(cid);
			hit.SetLayer(layer);
			hit.SetWire(wire);
			hit.SetDriftLength(dx);
			hit.SetTimeOffsetTrue(time);
			int c,n,a;
			if( confMan!=0 ){
				confMan->GetCounterMapManager()->GetCNA(cid, layer, wire, 1, 0, c, n, a);
				hit.SetCrate(c);
				hit.SetSlot(n);
				hit.SetChannel(a);

				hit.SetSimulatedResolution(confMan);
				hit.Reverse(confMan);
				double adl=hit.dl();
				hit.Calc(confMan);
				hit.SetDriftLength(adl);
			}
			blMan-> AddHit(hit);
		}
		else if( cid==CID_AC ){
		}
		else std::cout<<" !!! SimDataMan::Set() unknown CID="<<cid<<std::endl;
	}

	for( int i=0; i<mcData->trackSize(); i++ ){
		Track *track = mcData-> track(i);
		if( track-> parentTrackID()==0 && track->pdgID()==321 ){
			BeamMom = 0.001*track-> momentum().Mag();
		}
	}

	//  PrintEvent();
	return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::PrintHeader()
{
	tree2-> GetEntry(0);
	std::cout<<"===== Run Header MC Print ================================="<<std::endl;
	std::cout<<"= Seed : "<<runHeader->seed()<<std::endl;
	std::cout<<"= Num Of Event : "<<runHeader-> numEvent()<<std::endl;
	std::cout<<"=     Gnerated : "<<runHeader-> numGenerated()<<std::endl;
	std::cout<<"= Num Of Proces : "<<runHeader->CStable().CSSize()<<std::endl;
	std::cout<<"  print CS table (y/n)?"<<std::endl;
	while( true ){
		char in;
		std::cin>>in;
		if( in=='y' ){
			runHeader->CStable().PrintAllCS();
			break;
		}
		if( in=='n' ){
			break;
		}
	}
	std::cout<<"Num Of Event : "<<tree->GetEntries()<<std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::PrintEvent()
{
	std::cout<<"===== MC Event Print ==================================="<<std::endl;
	if( reactionData->ReactionID()==2 ) std::cout<<" K- 3He -> K- P d_s"<<std::endl;
	else std::cout<<"= ReactionID : "<<reactionData->ReactionID()<<std::endl;

	std::cout<<"= n Particle : "<<reactionData->ParticleSize()<<std::endl;
	for( int i=0; i<reactionData->ParticleSize(); i++ ){
		TVector3 mom = reactionData->GetParticle(i).Vect();
		std::cout<<"  > Particle : "<<PDGToName(reactionData->PDG(i))<<"  mom : ("<<mom.X()<<", "<<mom.Y()<<", "<<mom.Z()<<")"<<std::endl;
	}
	std::cout<<"= n Initi Particle : "<<reactionData->InitParticleSize()<<std::endl;
	for( int i=0; i<reactionData->InitParticleSize(); i++ ){
		TVector3 mom = reactionData->GetInitParticle(i).Vect();
		std::cout<<"  > Particle : "<<PDGToName(reactionData->InitPDG(i))<<"  mom : ("<<mom.X()<<", "<<mom.Y()<<", "<<mom.Z()<<")"<<std::endl;
	}
	std::cout<<"= n Track : "<<mcData->trackSize()<<std::endl;
#if 0
	for( int i=0; i<mcData->trackSize(); i++ ){
		PrintTrack(mcData->track(i));
	}
#endif
	std::cout<<"= n DetectorHit : "<<detectorData->detectorHitSize()<<std::endl;
#if 0
	for( int i=0; i<detectorData-> detectorHitSize(); i++ ){
		PrintDetectorHit(detectorData->detectorHit(i));
	}
#endif
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::PrintTrack(Track *track)
{
	TVector3 pos = track-> vertex();
	TVector3 mom = track-> momentum();
	std::cout<<"  > "<<PDGToName(track->pdgID())<<" trackID : "<<track->trackID()<<" parentID: "<<track-> parentTrackID()<<std::endl;
	std::cout<<"    mom : "<<std::setw(8)<<mom.Mag()<<"[GeV/c] ("<<mom.x()<<", "<<mom.y()<<", "<<mom.z()<<")"<<std::endl;
	std::cout<<"                     Vtx Pos : ("<<pos.x()<<", "<<pos.y()<<", "<<pos.z()<<")"<<std::endl;
}



// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::PrintDetectorHit(DetectorHit *hit)
{
	std::cout<<"  > "<<CIDToName(hit->detectorID())<<"  layer : "<<hit-> layerID()<<"  ch : "<<hit->channelID()<<"  time : "<<hit->time()<<"[ns]  dE : "<<hit->de()<<"[MeV]"<<std::endl;
	TVector3 pos = hit->pos();
	std::cout<<"    Hit Particle "<<PDGToName(hit->pdg())<<"  Hit Pos : ("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")"<<std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::PrintHit()
{
	int nCDClayer = 0;
	for( int i=0; i<NumOfCDCLayers; i++ ){
		if( nCDC[i]>0 ) nCDClayer++;
	}
	int nBLC2alayer = 0;
	int nBLC2blayer = 0;
	int nBPClayer   = 0;
	for( int i=0; i<8; i++ ){
		if( nBLC2a[i]>0 ) nBLC2alayer++;
		if( nBLC2b[i]>0 ) nBLC2blayer++;
		if( nBPC[i]>0 ) nBPClayer++;
	}
	int nFDC1layer = 0;
	for( int i=0; i<6; i++ ){
		if( nFDC1[i]>0 ) nFDC1layer++;
	}

	std::cout<<"===== Print Hodoscope Hit ============================="<<std::endl;
	std::cout<<"= nT0  : "<<nT0<<std::endl;
	std::cout<<"= nBPD : "<<nBPD<<std::endl;
	std::cout<<"= nDEF : "<<nDEF<<std::endl;
	std::cout<<"= nBVC : "<<nBVC<<std::endl;
	std::cout<<"= nCVC : "<<nCVC<<std::endl;
	std::cout<<"= nPC  : "<<nPC<<std::endl;
	std::cout<<"= nNC  : "<<nNC<<std::endl;
	std::cout<<std::endl;
	std::cout<<"= nCDH : "<<CDHHitPDG.size();
	for( int i=0; i<CDHHitPDG.size(); i++ ) std::cout<<" "<<PDGToName(CDHHitPDG[i]);
	std::cout<<std::endl;
	std::cout<<"= nIH  : "<<nIH<<std::endl;
	std::cout<<"===== Print Chamber Hit ==============================="<<std::endl;
	std::cout<<"= n Layer BLC2a : "<<nBLC2alayer<<std::endl;
	std::cout<<"= n Layer BLC2b : "<<nBLC2blayer<<std::endl;
	std::cout<<"= n Layer BPC   : "<<nBPClayer<<std::endl;
	std::cout<<"= n Layer FDC1  : "<<nFDC1layer<<std::endl;
	std::cout<<std::endl;
	std::cout<<"= n Layer CDC   : "<<nCDClayer<<std::endl;
	std::cout<<"======================================================="<<std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
std::string SimDataMan::PDGToName(const int &pdg_id) const
{
	if( pdg_id==11 )              return "e        ";
	else if( pdg_id==-11 )        return "e+       ";
	else if( pdg_id==22  )        return "gamma    ";
	else if( pdg_id==211 )        return "pi+      ";
	else if( pdg_id==-211 )       return "pi-      ";
	else if( pdg_id==111 )       return "pi0      ";
	else if( pdg_id==-321 )       return "K-       ";
	else if( pdg_id==321  )       return "K+       ";
	else if( pdg_id==311  )       return "K0       ";
	else if( pdg_id==130  )       return "K0l      ";
	else if( pdg_id==310  )       return "K0s      ";
	else if( pdg_id==2212 )       return "proton   ";
	else if( pdg_id==2112 )       return "neutron  ";
	else if( pdg_id==3122 )       return "Lambda   ";
	else if( pdg_id==3212 )       return "Sigma0   ";
	else if( pdg_id==3222 )       return "Sigma+   ";
	else if( pdg_id==3112 )       return "Sigma-   ";
	else if( pdg_id==13122 )       return "Lambda(1405)   ";
	else if( pdg_id==1000010020 ) return "deuteron ";
	else if( pdg_id==1000020030 ) return "He3      ";
	else{
		char c_str[512];
		sprintf(c_str,"%d",pdg_id);
		return std::string(c_str);
	}
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
std::string SimDataMan::CIDToName(const int &cid) const
{
	if( cid==CID_T0 )         return "T0   ";
	else if( cid==CID_BPD   ) return "BPD  ";
	else if( cid==CID_CDH   ) return "CDH  ";
	else if( cid==CID_BVC   ) return "BCV  ";
	else if( cid==CID_CVC   ) return "CVC  ";
	else if( cid==CID_PC    ) return "PC   ";
	else if( cid==CID_IH    ) return "IH   ";
	else if( cid==CID_BD    ) return "BD   ";
	else if( cid==CID_NC    ) return "NC   ";
	else if( cid==CID_BLC2b ) return "BLC2b";
	else if( cid==CID_BLC2a ) return "BLC2a";
	else if( cid==CID_BLC2  ) return "BLC2 ";
	else if( cid==CID_BPC   ) return "BPC  ";
	else if( cid==CID_FDC1  ) return "FDC1 ";
	else if( cid==CID_CDC   ) return "CDC  ";
	else if( cid==CID_AC    ) return "AC   ";
	else if( cid==CID_WVC   ) return "WVC  ";
	else{
		char c_str[512];
		sprintf(c_str,"%d",&cid);
		return std::string(c_str);
	}
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::Clear()
{
	CDHHitPDG.clear();
	nIH  = 0;
	nBPD = 0;
	nT0  = 0;
	nDEF = 0;
	nBVC = 0;
	nCVC = 0;
	nPC  = 0;
	nNC  = 0;

	T0Time -999.;
	BeamMom = -999.;

	for( int i=0; i<NumOfCDCLayers; i++ ) nCDC[i]=0;
	for( int i=0; i<8; i++ ){
		nBLC2a[i] = 0;
		nBLC2b[i] = 0;
		nBPC[i] = 0;
	}
	for( int i=0; i<6; i++ ) nFDC1[i] = 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SimDataMan::Check()
{
	TVector3 tot_init;
	TVector3 tot_final;
	for( int i=0; i<reactionData->InitParticleSize(); i++ ){
		tot_init += reactionData->GetInitParticle(i).Vect();
		TLorentzVector lmom = reactionData->GetInitParticle(i);
	}
	for( int i=0; i<reactionData->ParticleSize(); i++ ){
		tot_final += reactionData->GetParticle(i).Vect();
		TLorentzVector lmom = reactionData->GetParticle(i);
	}

	if( (tot_init-tot_final).Mag2()>0.0001 ){
		std::cout<<" !!! SimDataMan::Check not match init &final lorentz momentum !!!"<<std::endl;
		std::cout<<"     diff : "<<(tot_init-tot_final).Mag()<<std::endl;
		std::cout<<"Total Init Momentum : ("<<tot_init.X()<<", "<<tot_init.Y()<<", "<<tot_init.Z()<<")"<<std::endl;
		std::cout<<"Total Final Momentum : ("<<tot_final.X()<<", "<<tot_final.Y()<<", "<<tot_final.Z()<<")"<<std::endl;
		return false;
	}
	return true;
}
