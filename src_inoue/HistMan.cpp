#include "HistMan.h"

HistMan::HistMan(TFile *f) : fFile(f)
{
  clear();
  initTrig();
  initBLDC();
}

void HistMan::initTrig()
{
  fFile->cd();
  std::cout<<"===== HistMan::initTrig START ====="<<std::endl;
  new TH1F("TrigMode", "Trigger Mode", 20, 0.5, 20.5);
  for( int seg=1; seg<=16; seg++ ){
    new TH1F(Form("Trig%d_TDC", seg),  Form("Trigger %d TDC",seg),  4096, -0.5, 4095.5);
    new TH1F(Form("Trig%d_time", seg), Form("Trigger %d time",seg), 200, -100, 100);
  }
}

void HistMan::initBLDC()
{
  new TH1F("Kf_Reduction", "K/f Event Reduction", 20, -0.5, 19.5);
  new TH1F("N_Reduction", "Nutral Event Reduction", 20, -0.5, 19.5);

  new TH2F("ntrackBLC",   "ntrack BLC1 vs BLC2", 10, -0.5, 9.5, 10, -0.5, 9.5);
  new TH1F("ntrackBLC1",  "BLC1 ntrack",  10, -0.5, 9.5);
  new TH1F("ntrackBLC1a", "BLC1b ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1b", "BLC1a ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2",  "BLC2 ntrack",  10, -0.5, 9.5);
  new TH1F("ntrackBLC2a", "BLC2a ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2b", "BLC2b ntrack", 10, -0.5, 9.5);
  new TH1F("ntrackBPC",   "BPC ntrack",   10, -0.5, 9.5);
  new TH1F("ntrackFDC1",  "FDC1 ntrack",  10, -0.5, 9.5);

  new TH1F("statusBLC1",  "BLC1 status",  10, -0.5, 9.5);
  new TH1F("statusBLC1a", "BLC1b status", 10, -0.5, 9.5);
  new TH1F("statusBLC1b", "BLC1a status", 10, -0.5, 9.5);
  new TH1F("statusBLC2",  "BLC2 status",  10, -0.5, 9.5);
  new TH1F("statusBLC2a", "BLC2a status", 10, -0.5, 9.5);
  new TH1F("statusBLC2b", "BLC2b status", 10, -0.5, 9.5);
  new TH1F("statusBPC",   "BPC status",   10, -0.5, 9.5);
  new TH1F("statusFDC1",  "FDC1 stauts",  10, -0.5, 9.5);

  new TH1F("D5mom", "Beam momentum by D5", 1000, 0.5, 1.5);
  new TH1F("D5chi", "D5 Chi-Square", 1000, 0.0, 100);



  std::cout<<"===== HistMan::initTrig FINISH ====="<<std::endl;
}

void HistMan::fill(EventHeader *header, ConfMan *conf)
{
  TH1F *h1;
  h1 = (TH1F*)fFile-> Get("TrigMode"), h1-> Fill(header->trigmode(conf));
  for( int seg=1; seg<=16; seg++ ){
    h1 = (TH1F*)fFile-> Get(Form("Trig%d_TDC", seg)), h1-> Fill(header->pattern(seg));
    h1 = (TH1F*)fFile-> Get(Form("Trig%d_time", seg)), h1-> Fill(header->time(seg));
  }
}

void HistMan::fill(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamSpectrometer *D5)
{
  TH1F *h1;
  TH2F *h2;

  if( fT0hit.size()==1 ){
    fStatus=1;
    if( bltrackMan->ntrackBLC1()==1 ){
      fStatus=2;
      if( bltrackMan->ntrackBLC2()==1 ){
	fStatus=3;
	if( D5->chisquare()>0.0 ){
	  fStatus=4;
	  if( bltrackMan->ntrackBPC()==1 ){
	    fStatus=5;
	  }
	}
      }
    }
  }



  if( header-> IsTrig(Trig_Kf, conf) ){
    h1 = (TH1F*)fFile-> Get("Kf_Reduction"), h1->Fill(0);
    if( fT0hit.size()==1 ){
      h1->Fill(1);
      if( bltrackMan->ntrackBLC1()==1 ){
	h1->Fill(2);
	if( bltrackMan->ntrackBLC2()==1 ){
	  h1->Fill(3);
	  if( D5->chisquare()>0.0 ){
	    h1->Fill(4);
	    if( bltrackMan->ntrackBPC()==1 ){
	      h1->Fill(5);
	    }
	  }
	}
      }
    }
  }

  if( header-> IsTrig(Trig_Kf, conf) ){
    h1 = (TH1F*)fFile-> Get("Kf_Reduction"), h1->Fill(0);
    if( fT0hit.size()==1 ){
      h1->Fill(1);
      if( bltrackMan->ntrackBLC1()==1 ){
	h1->Fill(2);
	if( bltrackMan->ntrackBLC2()==1 ){
	  h1->Fill(3);
	  if( D5->chisquare()>0.0 ){
	    h1->Fill(4);
	    if( bltrackMan->ntrackBPC()==1 ){
	      h1->Fill(5);
	    }
	  }
	}
      }
    }
  }

  if( header-> IsTrig(Trig_Neutral, conf) ){
    h1 = (TH1F*)fFile-> Get("N_Reduction"), h1->Fill(0);
    if( fT0hit.size()==1 ){
      h1->Fill(1);
      if( bltrackMan->ntrackBLC1()==1 ){
	h1->Fill(2);
	if( bltrackMan->ntrackBLC2()==1 ){
	  h1->Fill(3);
	  if( D5->chisquare()>0.0 ){
	    h1->Fill(4);
	    if( bltrackMan->ntrackBPC()==1 ){
	      h1->Fill(5);
	    }
	  }
	}
      }
    }
  }

  h1 = (TH1F*)fFile-> Get("ntrackBLC1"),  h1-> Fill(bltrackMan->ntrackBLC1());
  h1 = (TH1F*)fFile-> Get("ntrackBLC1a"), h1-> Fill(bltrackMan->ntrackBLC1a());
  h1 = (TH1F*)fFile-> Get("ntrackBLC1b"), h1-> Fill(bltrackMan->ntrackBLC1b());
  h1 = (TH1F*)fFile-> Get("ntrackBLC2"),  h1-> Fill(bltrackMan->ntrackBLC2());
  h1 = (TH1F*)fFile-> Get("ntrackBLC2a"), h1-> Fill(bltrackMan->ntrackBLC2a());
  h1 = (TH1F*)fFile-> Get("ntrackBLC2b"), h1-> Fill(bltrackMan->ntrackBLC2b());
  h1 = (TH1F*)fFile-> Get("ntrackBPC"),   h1-> Fill(bltrackMan->ntrackBPC());
  h1 = (TH1F*)fFile-> Get("ntrackFDC1"),  h1-> Fill(bltrackMan->ntrackFDC1());
  h2 = (TH2F*)fFile-> Get("ntrackBLC"),   h2-> Fill(bltrackMan->ntrackBLC1(), bltrackMan->ntrackBLC2());

  h1 = (TH1F*)fFile-> Get("statusBLC1"),  h1-> Fill(bltrackMan->status(CID_BLC1));
  h1 = (TH1F*)fFile-> Get("statusBLC1a"), h1-> Fill(bltrackMan->status(CID_BLC1a));
  h1 = (TH1F*)fFile-> Get("statusBLC1b"), h1-> Fill(bltrackMan->status(CID_BLC1b));
  h1 = (TH1F*)fFile-> Get("statusBLC2"),  h1-> Fill(bltrackMan->status(CID_BLC2));
  h1 = (TH1F*)fFile-> Get("statusBLC2a"), h1-> Fill(bltrackMan->status(CID_BLC2a));
  h1 = (TH1F*)fFile-> Get("statusBLC2b"), h1-> Fill(bltrackMan->status(CID_BLC2b));
  h1 = (TH1F*)fFile-> Get("statusBPC"),   h1-> Fill(bltrackMan->status(CID_BPC));
  h1 = (TH1F*)fFile-> Get("statusFDC1"),  h1-> Fill(bltrackMan->status(CID_FDC1));

  if( D5->chisquare()>0 ){
    h1 = (TH1F*)fFile-> Get("D5mom"), h1->Fill(D5->mom());
    h1 = (TH1F*)fFile-> Get("D5chi"), h1->Fill(D5->chisquare());
  }
}

void HistMan::set(BeamLineHitMan *blMan)
{
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      fBHDhit.push_back(blMan->BHD(i));
    }
  }
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      fT0hit.push_back(blMan->T0(i));
    }
  }
  for( int i=0; i<blMan->nDEF(); i++ ){
    if( blMan->DEF(i)->CheckRange() ){
      fDEFhit.push_back(blMan->DEF(i));
    }
  }
  for( int i=0; i<blMan->nBPD(); i++ ){
    if( blMan->BPD(i)->CheckRange() ){
      fBPDhit.push_back(blMan->BPD(i));
    }
  }
  for( int i=0; i<blMan->nCVC(); i++ ){
    if( blMan->CVC(i)->CheckRange() ){
      fCVChit.push_back(blMan->CVC(i));
    }
  }
  for( int i=0; i<blMan->nPC(); i++ ){
    if( blMan->PC(i)->CheckRange() ){
      fPChit.push_back(blMan->PC(i));
    }
  }
  for( int i=0; i<blMan->nBD(); i++ ){
    if( blMan->BD(i)->CheckRange() ){
      fBDhit.push_back(blMan->BD(i));
    }
  }

}

void HistMan::clear()
{
  fBHDhit.clear();
  fT0hit.clear();
  fBPDhit.clear();
  fDEFhit.clear();
  fCVChit.clear();
  fPChit.clear();
  fBVChit.clear();
  fBDhit.clear();
  for( int lay=0; lay<8; lay++ ) fNChit[lay].clear();

}
