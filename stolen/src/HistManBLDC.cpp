#include "HistManBeamAna.h"

void HistManBeamAna::initBLDC()
{
  fFile-> cd();
  std::cout<<"BLDC Histogram Initialize ... ";

  new TH1F("BLC1a_eff", "BLC1a efficiency", 10, -0.5, 9.5);
  new TH1F("BLC1b_eff", "BLC1a efficiency", 10, -0.5, 9.5);
  new TH1F("BLC2a_eff", "BLC1a efficiency", 10, -0.5, 9.5);
  new TH1F("BLC2b_eff", "BLC1a efficiency", 10, -0.5, 9.5);
  new TH1F("BPC_eff",   "BPC efficiency",   10, -0.5, 9.5);
  new TH1F("FDC1_eff",  "FDC1 efficiency",  10, -0.5, 9.5);
  new TH1F("BLC1a_lay_eff", "BLC1a layer efficiency", 9, -0.5, 8.5);
  new TH1F("BLC1b_lay_eff", "BLC1a layer efficiency", 9, -0.5, 8.5);
  new TH1F("BLC2a_lay_eff", "BLC1a layer efficiency", 9, -0.5, 8.5);
  new TH1F("BLC2b_lay_eff", "BLC1a layer efficiency", 9, -0.5, 8.5);
  new TH1F("BPC_lay_eff",   "BPC layer efficiency", 9 , -0.5, 8.5);
  new TH1F("FDC1_lay_eff",  "FDC1 layer efficiency", 7 , -0.5, 6.5);
  new TH1F("FDC1_status",  "FDC1 status",  11, -0.5, 10.5);

  new TH1F("BLC1a_eff_18", "BLC1a efficiency 1&8 layer hit", 10, -0.5, 9.5);
  new TH1F("BLC1b_eff_18", "BLC1b efficiency 1&8 layer hit", 10, -0.5, 9.5);
  new TH1F("BLC2a_eff_18", "BLC2a efficiency 1&8 layer hit", 10, -0.5, 9.5);
  new TH1F("BLC2b_eff_18", "BLC2b efficiency 1&8 layer hit", 10, -0.5, 9.5);
  new TH1F("BPC_eff_18",   "BPC efficiency 1&8 layer hit", 10, -0.5, 9.5);
  new TH1F("FDC1_eff_18",  "FDC1 efficiency 1&8 layer hit", 8, -0.5, 7.5);

  for( int lay=1; lay<=8; lay++ ){
    new TH1F(Form("mul_BLC1a_%d", lay), Form("BLC1a multiplicity layer%d", lay),  33, -0.5, 32.5);
    new TH1F(Form("hitpat_BLC1a_%d", lay), Form("BLC1a hit pttern layer%d", lay), 32, 0.5, 32.5);
    new TH1F(Form("BLC1a_tdc_%d", lay), Form("BLC1a layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC1a_dt_%d", lay), Form("BLC1a layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC1a_dl_%d", lay), Form("BLC1a layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC1a_res_%d", lay), Form("BLC1a layer%d residual", lay), 1000, -0.5, 0.5);

    new TH1F(Form("BLC1a_tdc_%d_1hit", lay), Form("BLC1a layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC1a_dt_%d_1hit", lay), Form("BLC1a layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC1a_dl_%d_1hit", lay), Form("BLC1a layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC1a_res_%d_1hit", lay), Form("BLC1a layer%d residual", lay), 1000, -0.5, 0.5);
    for( int wire=1; wire<=32; wire++ ){
      new TH1F(Form("BLC1a_tdc_%d_%d", lay, wire), Form("BLC1a layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC1a_dt_%d_%d", lay, wire), Form("BLC1a layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC1a_dl_%d_%d", lay, wire), Form("BLC1a layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC1a_res_%d_%d", lay, wire), Form("BLC1a layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);

      new TH1F(Form("BLC1a_tdc_%d_%d_1hit", lay, wire), Form("BLC1a layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC1a_dt_%d_%d_1hit", lay, wire), Form("BLC1a layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC1a_dl_%d_%d_1hit", lay, wire), Form("BLC1a layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC1a_res_%d_%d_1hit", lay, wire), Form("BLC1a layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);
    }

    new TH1F(Form("mul_BLC1b_%d", lay), Form("BLC1b multiplicity layer%d", lay),  33, -0.5, 32.5);
    new TH1F(Form("hitpat_BLC1b_%d", lay), Form("BLC1b hit pttern layer%d", lay), 32, 0.5, 32.5);
    new TH1F(Form("BLC1b_tdc_%d", lay), Form("BLC1b layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC1b_dt_%d", lay), Form("BLC1b layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC1b_dl_%d", lay), Form("BLC1b layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC1b_res_%d", lay), Form("BLC1b layer%d residual", lay), 1000, -0.5, 0.5);

    new TH1F(Form("BLC1b_tdc_%d_1hit", lay), Form("BLC1b layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC1b_dt_%d_1hit", lay), Form("BLC1b layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC1b_dl_%d_1hit", lay), Form("BLC1b layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC1b_res_%d_1hit", lay), Form("BLC1b layer%d residual", lay), 1000, -0.5, 0.5);
    for( int wire=1; wire<=32; wire++ ){
      new TH1F(Form("BLC1b_tdc_%d_%d", lay, wire), Form("BLC1b layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC1b_dt_%d_%d", lay, wire), Form("BLC1b layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC1b_dl_%d_%d", lay, wire), Form("BLC1b layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC1b_res_%d_%d", lay, wire), Form("BLC1b layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);

      new TH1F(Form("BLC1b_tdc_%d_%d_1hit", lay, wire), Form("BLC1b layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC1b_dt_%d_%d_1hit", lay, wire), Form("BLC1b layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC1b_dl_%d_%d_1hit", lay, wire), Form("BLC1b layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC1b_res_%d_%d_1hit", lay, wire), Form("BLC1b layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);
    }

    new TH1F(Form("mul_BLC2a_%d", lay), Form("BLC2a multiplicity layer%d", lay),  33, -0.5, 32.5);
    new TH1F(Form("hitpat_BLC2a_%d", lay), Form("BLC2a hit pttern layer%d", lay), 32,  0.5, 32.5);
    new TH1F(Form("BLC2a_tdc_%d", lay), Form("BLC2a layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC2a_dt_%d", lay), Form("BLC2a layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC2a_dl_%d", lay), Form("BLC2a layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC2a_res_%d", lay), Form("BLC2a layer%d residual", lay), 1000, -0.5, 0.5);

    new TH1F(Form("BLC2a_tdc_%d_1hit", lay), Form("BLC2a layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC2a_dt_%d_1hit", lay), Form("BLC2a layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC2a_dl_%d_1hit", lay), Form("BLC2a layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC2a_res_%d_1hit", lay), Form("BLC2a layer%d residual", lay), 1000, -0.5, 0.5);
    for( int wire=1; wire<=32; wire++ ){
      new TH1F(Form("BLC2a_tdc_%d_%d", lay, wire), Form("BLC2a layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC2a_dt_%d_%d", lay, wire), Form("BLC2a layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC2a_dl_%d_%d", lay, wire), Form("BLC2a layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC2a_res_%d_%d", lay, wire), Form("BLC2a layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);

      new TH1F(Form("BLC2a_tdc_%d_%d_1hit", lay, wire), Form("BLC2a layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC2a_dt_%d_%d_1hit", lay, wire), Form("BLC2a layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC2a_dl_%d_%d_1hit", lay, wire), Form("BLC2a layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC2a_res_%d_%d_1hit", lay, wire), Form("BLC2a layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);
    }

    new TH1F(Form("mul_BLC2b_%d", lay), Form("BLC2b multiplicity layer%d", lay),  33, -0.5, 32.5);
    new TH1F(Form("hitpat_BLC2b_%d", lay), Form("BLC2b hit pttern layer%d", lay), 32,  0.5, 32.5);
    new TH1F(Form("BLC2b_tdc_%d", lay), Form("BLC2b layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC2b_dt_%d", lay), Form("BLC2b layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC2b_dl_%d", lay), Form("BLC2b layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC2b_res_%d", lay), Form("BLC2b layer%d residual", lay), 1000, -0.5, 0.5);

    new TH1F(Form("BLC2b_tdc_%d_1hit", lay), Form("BLC2b layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BLC2b_dt_%d_1hit", lay), Form("BLC2b layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BLC2b_dl_%d_1hit", lay), Form("BLC2b layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BLC2b_res_%d_1hit", lay), Form("BLC2b layer%d residual", lay), 1000, -0.5, 0.5);
    for( int wire=1; wire<=32; wire++ ){
      new TH1F(Form("BLC2b_tdc_%d_%d", lay, wire), Form("BLC2b layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC2b_dt_%d_%d", lay, wire), Form("BLC2b layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC2b_dl_%d_%d", lay, wire), Form("BLC2b layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC2b_res_%d_%d", lay, wire), Form("BLC2b layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);

      new TH1F(Form("BLC2b_tdc_%d_%d_1hit", lay, wire), Form("BLC2b layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BLC2b_dt_%d_%d_1hit", lay, wire), Form("BLC2b layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BLC2b_dl_%d_%d_1hit", lay, wire), Form("BLC2b layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BLC2b_res_%d_%d_1hit", lay, wire), Form("BLC2b layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);
    }

    new TH1F(Form("mul_BPC_%d", lay), Form("BPC multiplicity layer%d", lay), 16, -0.5, 15.5);
    new TH1F(Form("hitpat_BPC_%d", lay), Form("BPC hit pttern layer%d", lay), 15, 0.5, 15.5);
    new TH1F(Form("BPC_tdc_%d", lay), Form("BPC layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BPC_dt_%d", lay), Form("BPC layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BPC_dl_%d", lay), Form("BPC layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BPC_res_%d", lay), Form("BPC layer%d residual", lay), 1000, -0.5, 0.5);

    new TH1F(Form("BPC_tdc_%d_1hit", lay), Form("BPC layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("BPC_dt_%d_1hit", lay), Form("BPC layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("BPC_dl_%d_1hit", lay), Form("BPC layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("BPC_res_%d_1hit", lay), Form("BPC layer%d residual", lay), 1000, -0.5, 0.5);
    for( int wire=1; wire<=15; wire++ ){
      new TH1F(Form("BPC_tdc_%d_%d", lay, wire), Form("BPC layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BPC_dt_%d_%d", lay, wire), Form("BPC layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BPC_dl_%d_%d", lay, wire), Form("BPC layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BPC_res_%d_%d", lay, wire), Form("BPC layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);

      new TH1F(Form("BPC_tdc_%d_%d_1hit", lay, wire), Form("BPC layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("BPC_dt_%d_%d_1hit", lay, wire), Form("BPC layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("BPC_dl_%d_%d_1hit", lay, wire), Form("BPC layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("BPC_res_%d_%d_1hit", lay, wire), Form("BPC layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);
    }
  }

  for( int lay=1; lay<=6; lay++ ){
    new TH1F(Form("mul_FDC1_%d", lay),     Form("FDC1 multiplicity layer%d", lay), 65, -0.5, 64.5);
    new TH1F(Form("hitpat_FDC1_%d", lay),  Form("FDC1 hit pttern layer%d", lay), 64, 0.5, 64.5);
    new TH1F(Form("FDC1_tdc_%d", lay),     Form("FDC1 layer%d TDC", lay), 4000, 0, 4000);
    new TH1F(Form("FDC1_dt_%d", lay),      Form("FDC1 layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("FDC1_dl_%d", lay),      Form("FDC1 layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("FDC1_res_%d", lay),     Form("FDC1 layer%d residual", lay), 1000, -0.5, 0.5);

    new TH1F(Form("FDC1_tdc_%d_1hit", lay), Form("FDC1 layer%d TDC", lay), 4000, 0, 4000); 
    new TH1F(Form("FDC1_dt_%d_1hit", lay), Form("FDC1 layer%d dt", lay), 500, -100, 400);
    new TH1F(Form("FDC1_dl_%d_1hit", lay), Form("FDC1 layer%d dl", lay), 500, -0.5, 4.5);
    new TH1F(Form("FDC1_res_%d_1hit", lay), Form("FDC1 layer%d residual", lay), 1000, -0.5, 0.5);
    for( int wire=1; wire<=64; wire++ ){
      new TH1F(Form("FDC1_tdc_%d_%d", lay, wire), Form("FDC1 layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("FDC1_dt_%d_%d", lay, wire), Form("FDC1 layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("FDC1_dl_%d_%d", lay, wire), Form("FDC1 layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("FDC1_res_%d_%d", lay, wire), Form("FDC1 layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);

      new TH1F(Form("FDC1_tdc_%d_%d_1hit", lay, wire), Form("FDC1 layer%d wire%d TDC", lay, wire), 4000, 0, 4000); 
      new TH1F(Form("FDC1_dt_%d_%d_1hit", lay, wire), Form("FDC1 layer%d wire%d dt", lay, wire), 500, -100, 400);
      new TH1F(Form("FDC1_dl_%d_%d_1hit", lay, wire), Form("FDC1 layer%d wire%d dl", lay, wire), 500, -0.5, 4.5);
      new TH1F(Form("FDC1_res_%d_%d_1hit", lay, wire), Form("FDC1 layer%d wire%d residual", lay, wire), 1000, -0.5, 0.5);

    }
  }

  new TH2F("prof_atBPC_byBLC2", "Beam Profile at BPC by BLC2", 1000, -25, 25, 1000, -25, 25); 

  new TH2F("BLC2BPC_pos_diff", "BLC2-BPC position diff", 1000, -10, 10, 1000, -10, 10);
  new TH2F("BLC2BPC_dir_diff", "BLC2-BPC direction diff", 1000, -0.1, 0.1, 1000, -0.1, 0.1);

  new TH2F("profFF_byBPC", "Beam Profile FF", 3000, -15, 15, 3000, -15, 15);
  new TH2F("profFF_byBPC_CDH2", "Beam Profile FF", 3000, -15, 15, 3000, -15, 15);
  new TH2F("profFF_byBPC_K", "Beam Profile FF", 3000, -15, 15, 3000, -15, 15);

  std::cout<<"finish"<<std::endl;
}

void HistManBeamAna::fillBLDC(EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  TH1 *h1;
  TH2 *h2;
  bool BLC1a_trig_18 = false;
  bool BLC1b_trig_18 = false;
  bool BLC2a_trig_18 = false;
  bool BLC2b_trig_18 = false;
  bool BPC_trig_18 = false;
  bool FDC1_trig_16 = false;

  h1 = getTH("FDC1_status"), h1-> Fill(bltrackMan->status(CID_FDC1));

  {
    bool BLC1a_all_1hit = true;
    for( int lay=1; lay<=8; lay++ ){
      if( blMan->nBLC1a(lay)!=1 ) BLC1a_all_1hit = false;
    }

    if( blMan-> nBLC1a(1)==1 && blMan-> nBLC1a(8)==1 ){
      int wire1 = blMan->BLC1a(1, 0)-> wire();
      int wire2 = blMan->BLC1a(8, 0)-> wire();
      if( 1<wire1 && wire1<32 && 1<wire2 && wire2<32 ) BLC1a_trig_18 = true;
    }
    TH1 *h1_eff = getTH("BLC1a_eff_18");
    h1_eff-> Fill(0);
    for( int lay=1; lay<=8; lay++ ){
      h1 = getTH(Form("mul_BLC1a_%d", lay)), h1-> Fill(blMan->nBLC1a(lay));
      for( int i=0; i<blMan->nBLC1a(lay); i++ ){
	ChamberLikeHit *hit = blMan-> BLC1a(lay, i);
	int tdc = hit-> tdc(), wire = hit-> wire();
	double dt = hit-> dt();
	double dl = hit-> dl();
	
	h1 = getTH(Form("hitpat_BLC1a_%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("BLC1a_tdc_%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC1a_dt_%d", lay)), h1-> Fill(dt);
	h1 = getTH(Form("BLC1a_dl_%d", lay)), h1-> Fill(dl);

	h1 = getTH(Form("BLC1a_tdc_%d_%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC1a_dt_%d_%d", lay, wire)), h1-> Fill(dt);
	h1 = getTH(Form("BLC1a_dl_%d_%d", lay, wire)), h1-> Fill(dl);
	if( BLC1a_all_1hit ){
	  h1 = getTH(Form("BLC1a_tdc_%d_1hit", lay)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC1a_dt_%d_1hit", lay)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC1a_dl_%d_1hit", lay)), h1-> Fill(dl);

	  h1 = getTH(Form("BLC1a_tdc_%d_%d_1hit", lay, wire)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC1a_dt_%d_%d_1hit", lay, wire)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC1a_dl_%d_%d_1hit", lay, wire)), h1-> Fill(dl);
	}
      }
      if( BLC1a_trig_18 ) if( 0<blMan->nBLC1a(lay) && blMan->nBLC1a(lay)<3 ) h1_eff->Fill(lay);
    }
    if( bltrackMan-> ntrackBLC1a()==1 ){
      LocalTrack *track = bltrackMan->trackBLC1a(0);
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *hit = track->hit(i);
	int layer = hit-> layer(), wire = hit->wire();
	double dltrack;
	if( hit->xy()==0 )      dltrack = hit->x()-hit->wx();
	else if( hit->xy()==1 ) dltrack = hit->y()-hit->wy();

	h1 = getTH(Form("BLC1a_res_%d", layer)), h1-> Fill(hit->resl());
	h1 = getTH(Form("BLC1a_res_%d_%d", layer, wire)), h1-> Fill(hit->resl());
	if( BLC1a_all_1hit ){
	  h1 = getTH(Form("BLC1a_res_%d_1hit", layer)), h1-> Fill(hit->resl());
	  h1 = getTH(Form("BLC1a_res_%d_%d_1hit", layer, wire)), h1-> Fill(hit->resl());
	}
      }
    }

    if( BLC1a_trig_18 ) if( bltrackMan->ntrackBLC1a()>0 ) h1_eff-> Fill(9);
  }
  {
    bool BLC1b_all_1hit = true;
    for( int lay=1; lay<=8; lay++ ){
      if( blMan->nBLC1b(lay)!=1 ) BLC1b_all_1hit = false;
    }

    if( blMan-> nBLC1b(1)==1 && blMan-> nBLC1b(8)==1 ){
      int wire1 = blMan->BLC1b(1, 0)-> wire();
      int wire2 = blMan->BLC1b(8, 0)-> wire();
      if( 1<wire1 && wire1<32 && 1<wire2 && wire2<32 ) BLC1b_trig_18 = true;
    }
    TH1 *h1_eff = getTH("BLC1b_eff_18");
    h1_eff-> Fill(0);
    for( int lay=1; lay<=8; lay++ ){
      h1 = getTH(Form("mul_BLC1b_%d", lay)), h1-> Fill(blMan->nBLC1b(lay));
      for( int i=0; i<blMan->nBLC1b(lay); i++ ){
	ChamberLikeHit *hit = blMan-> BLC1b(lay, i);
	int tdc = hit-> tdc(), wire = hit-> wire();
	double dt = hit-> dt();
	double dl = hit-> dl();
	
	h1 = getTH(Form("hitpat_BLC1b_%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("BLC1b_tdc_%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC1b_dt_%d", lay)), h1-> Fill(dt);
	h1 = getTH(Form("BLC1b_dl_%d", lay)), h1-> Fill(dl);

	h1 = getTH(Form("BLC1b_tdc_%d_%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC1b_dt_%d_%d", lay, wire)), h1-> Fill(dt);
	h1 = getTH(Form("BLC1b_dl_%d_%d", lay, wire)), h1-> Fill(dl);
	if( BLC1b_all_1hit ){
	  h1 = getTH(Form("BLC1b_tdc_%d_1hit", lay)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC1b_dt_%d_1hit", lay)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC1b_dl_%d_1hit", lay)), h1-> Fill(dl);

	  h1 = getTH(Form("BLC1b_tdc_%d_%d_1hit", lay, wire)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC1b_dt_%d_%d_1hit", lay, wire)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC1b_dl_%d_%d_1hit", lay, wire)), h1-> Fill(dl);
	}
      }
      if( BLC1b_trig_18 ) if( 0<blMan->nBLC1b(lay) && blMan->nBLC1b(lay)<3 ) h1_eff->Fill(lay);
    }
    if( bltrackMan-> ntrackBLC1b()==1 ){
      LocalTrack *track = bltrackMan->trackBLC1b(0);
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *hit = track->hit(i);
	int layer = hit-> layer(), wire = hit->wire();
	double dltrack;
	if( hit->xy()==0 )      dltrack = hit->x()-hit->wx();
	else if( hit->xy()==1 ) dltrack = hit->y()-hit->wy();

	h1 = getTH(Form("BLC1b_res_%d", layer)), h1-> Fill(hit->resl());
	h1 = getTH(Form("BLC1b_res_%d_%d", layer, wire)), h1-> Fill(hit->resl());
	if( BLC1b_all_1hit ){
	  h1 = getTH(Form("BLC1b_res_%d_1hit", layer)), h1-> Fill(hit->resl());
	  h1 = getTH(Form("BLC1b_res_%d_%d_1hit", layer, wire)), h1-> Fill(hit->resl());
	}
      }
    }

    if( BLC1b_trig_18 ) if( bltrackMan->ntrackBLC1b()>0 ) h1_eff-> Fill(9);
  }

  {
    bool BLC2a_all_1hit = true;
    for( int lay=1; lay<=8; lay++ ){
      if( blMan->nBLC2a(lay)!=1 ) BLC2a_all_1hit = false;
    }

    if( blMan-> nBLC2a(1)==1 && blMan-> nBLC2a(8)==1 ){
      int wire1 = blMan->BLC2a(1, 0)-> wire();
      int wire2 = blMan->BLC2a(8, 0)-> wire();
      if( 1<wire1 && wire1<32 && 1<wire2 && wire2<32 ) BLC2a_trig_18 = true;
    }
    TH1 *h1_eff = getTH("BLC2a_eff_18");
    h1_eff-> Fill(0);
    for( int lay=1; lay<=8; lay++ ){
      h1 = getTH(Form("mul_BLC2a_%d", lay)), h1-> Fill(blMan->nBLC2a(lay));
      for( int i=0; i<blMan->nBLC2a(lay); i++ ){
	ChamberLikeHit *hit = blMan-> BLC2a(lay, i);
	int tdc = hit-> tdc(), wire = hit-> wire();
	double dt = hit-> dt();
	double dl = hit-> dl();
	
	h1 = getTH(Form("hitpat_BLC2a_%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("BLC2a_tdc_%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC2a_dt_%d", lay)), h1-> Fill(dt);
	h1 = getTH(Form("BLC2a_dl_%d", lay)), h1-> Fill(dl);

	h1 = getTH(Form("BLC2a_tdc_%d_%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC2a_dt_%d_%d", lay, wire)), h1-> Fill(dt);
	h1 = getTH(Form("BLC2a_dl_%d_%d", lay, wire)), h1-> Fill(dl);
	if( BLC2a_all_1hit ){
	  h1 = getTH(Form("BLC2a_tdc_%d_1hit", lay)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC2a_dt_%d_1hit", lay)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC2a_dl_%d_1hit", lay)), h1-> Fill(dl);

	  h1 = getTH(Form("BLC2a_tdc_%d_%d_1hit", lay, wire)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC2a_dt_%d_%d_1hit", lay, wire)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC2a_dl_%d_%d_1hit", lay, wire)), h1-> Fill(dl);
	}
      }
      if( BLC2a_trig_18 ) if( 0<blMan->nBLC2a(lay) && blMan->nBLC2a(lay)<3 ) h1_eff->Fill(lay);
    }
    if( bltrackMan-> ntrackBLC2a()==1 ){
      LocalTrack *track = bltrackMan->trackBLC2a(0);
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *hit = track->hit(i);
	int layer = hit-> layer(), wire = hit->wire();
	double dltrack;
	if( hit->xy()==0 )      dltrack = hit->x()-hit->wx();
	else if( hit->xy()==1 ) dltrack = hit->y()-hit->wy();

	h1 = getTH(Form("BLC2a_res_%d", layer)), h1-> Fill(hit->resl());
	h1 = getTH(Form("BLC2a_res_%d_%d", layer, wire)), h1-> Fill(hit->resl());
	if( BLC2a_all_1hit ){
	  h1 = getTH(Form("BLC2a_res_%d_1hit", layer)), h1-> Fill(hit->resl());
	  h1 = getTH(Form("BLC2a_res_%d_%d_1hit", layer, wire)), h1-> Fill(hit->resl());
	}
      }
    }

    if( BLC2a_trig_18 ) if( bltrackMan->ntrackBLC2a()>0 ) h1_eff-> Fill(9);
  }
  {
    bool BLC2b_all_1hit = true;
    for( int lay=1; lay<=8; lay++ ){
      if( blMan->nBLC2b(lay)!=1 ) BLC2b_all_1hit = false;
    }

    if( blMan-> nBLC2b(1)==1 && blMan-> nBLC2b(8)==1 ){
      int wire1 = blMan->BLC2b(1, 0)-> wire();
      int wire2 = blMan->BLC2b(8, 0)-> wire();
      if( 1<wire1 && wire1<32 && 1<wire2 && wire2<32 ) BLC2b_trig_18 = true;
    }
    TH1 *h1_eff = getTH("BLC2b_eff_18");
    h1_eff-> Fill(0);
    for( int lay=1; lay<=8; lay++ ){
      h1 = getTH(Form("mul_BLC2b_%d", lay)), h1-> Fill(blMan->nBLC2b(lay));
      for( int i=0; i<blMan->nBLC2b(lay); i++ ){
	ChamberLikeHit *hit = blMan-> BLC2b(lay, i);
	int tdc = hit-> tdc(), wire = hit-> wire();
	double dt = hit-> dt();
	double dl = hit-> dl();
	
	h1 = getTH(Form("hitpat_BLC2b_%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("BLC2b_tdc_%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC2b_dt_%d", lay)), h1-> Fill(dt);
	h1 = getTH(Form("BLC2b_dl_%d", lay)), h1-> Fill(dl);

	h1 = getTH(Form("BLC2b_tdc_%d_%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("BLC2b_dt_%d_%d", lay, wire)), h1-> Fill(dt);
	h1 = getTH(Form("BLC2b_dl_%d_%d", lay, wire)), h1-> Fill(dl);
	if( BLC2b_all_1hit ){
	  h1 = getTH(Form("BLC2b_tdc_%d_1hit", lay)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC2b_dt_%d_1hit", lay)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC2b_dl_%d_1hit", lay)), h1-> Fill(dl);

	  h1 = getTH(Form("BLC2b_tdc_%d_%d_1hit", lay, wire)), h1-> Fill(tdc);
	  h1 = getTH(Form("BLC2b_dt_%d_%d_1hit", lay, wire)), h1-> Fill(dt);
	  h1 = getTH(Form("BLC2b_dl_%d_%d_1hit", lay, wire)), h1-> Fill(dl);
	}
      }
      if( BLC2b_trig_18 ) if( 0<blMan->nBLC2b(lay) && blMan->nBLC2b(lay)<3 ) h1_eff->Fill(lay);
    }
    if( bltrackMan-> ntrackBLC2b()==1 ){
      LocalTrack *track = bltrackMan->trackBLC2b(0);
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *hit = track->hit(i);
	int layer = hit-> layer(), wire = hit->wire();
	double dltrack;
	if( hit->xy()==0 )      dltrack = hit->x()-hit->wx();
	else if( hit->xy()==1 ) dltrack = hit->y()-hit->wy();

	h1 = getTH(Form("BLC2b_res_%d", layer)), h1-> Fill(hit->resl());
	h1 = getTH(Form("BLC2b_res_%d_%d", layer, wire)), h1-> Fill(hit->resl());
	if( BLC2b_all_1hit ){
	  h1 = getTH(Form("BLC2b_res_%d_1hit", layer)), h1-> Fill(hit->resl());
	  h1 = getTH(Form("BLC2b_res_%d_%d_1hit", layer, wire)), h1-> Fill(hit->resl());
	}
      }
    }

    if( BLC2b_trig_18 ) if( bltrackMan->ntrackBLC2b()>0 ) h1_eff-> Fill(9);
  }

  {
    bool BPC_all_1hit = true;
    for( int lay=1; lay<=8; lay++ ){
      if( blMan->nBPC(lay)!=1 ) BPC_all_1hit = false;
    }

    if( blMan-> nBPC(1)==1 && blMan-> nBPC(8)==1 ){
      int wire1 = blMan->BPC(1, 0)-> wire();
      int wire2 = blMan->BPC(8, 0)-> wire();
      if( 1<wire1 && wire1<32 && 1<wire2 && wire2<32 ) BPC_trig_18 = true;
    }
    TH1 *h1_eff = getTH("BPC_eff_18");
    h1_eff-> Fill(0);
    for( int lay=1; lay<=8; lay++ ){
      h1 = getTH(Form("mul_BPC_%d", lay)), h1-> Fill(blMan->nBPC(lay));
      for( int i=0; i<blMan->nBPC(lay); i++ ){
	ChamberLikeHit *hit = blMan-> BPC(lay, i);
	int tdc = hit-> tdc(), wire = hit-> wire();
	double dt = hit-> dt();
	double dl = hit-> dl();
	
	h1 = getTH(Form("hitpat_BPC_%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("BPC_tdc_%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("BPC_dt_%d", lay)), h1-> Fill(dt);
	h1 = getTH(Form("BPC_dl_%d", lay)), h1-> Fill(dl);

	h1 = getTH(Form("BPC_tdc_%d_%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("BPC_dt_%d_%d", lay, wire)), h1-> Fill(dt);
	h1 = getTH(Form("BPC_dl_%d_%d", lay, wire)), h1-> Fill(dl);
	if( BPC_all_1hit ){
	  h1 = getTH(Form("BPC_tdc_%d_1hit", lay)), h1-> Fill(tdc);
	  h1 = getTH(Form("BPC_dt_%d_1hit", lay)), h1-> Fill(dt);
	  h1 = getTH(Form("BPC_dl_%d_1hit", lay)), h1-> Fill(dl);

	  h1 = getTH(Form("BPC_tdc_%d_%d_1hit", lay, wire)), h1-> Fill(tdc);
	  h1 = getTH(Form("BPC_dt_%d_%d_1hit", lay, wire)), h1-> Fill(dt);
	  h1 = getTH(Form("BPC_dl_%d_%d_1hit", lay, wire)), h1-> Fill(dl);
	}
      }
      if( BPC_trig_18 ) if( 0<blMan->nBPC(lay) && blMan->nBPC(lay)<3 ) h1_eff->Fill(lay);
    }
    if( bltrackMan-> ntrackBPC()==1 ){
      LocalTrack *track = bltrackMan->trackBPC(0);
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *hit = track->hit(i);
	int layer = hit-> layer(), wire = hit->wire();
	double dltrack;
	if( hit->xy()==0 )      dltrack = hit->x()-hit->wx();
	else if( hit->xy()==1 ) dltrack = hit->y()-hit->wy();

	h1 = getTH(Form("BPC_res_%d", layer)), h1-> Fill(hit->resl());
	h1 = getTH(Form("BPC_res_%d_%d", layer, wire)), h1-> Fill(hit->resl());
	if( BPC_all_1hit ){
	  h1 = getTH(Form("BPC_res_%d_1hit", layer)), h1-> Fill(hit->resl());
	  h1 = getTH(Form("BPC_res_%d_%d_1hit", layer, wire)), h1-> Fill(hit->resl());
	}
      }
    }

    if( BPC_trig_18 ) if( bltrackMan->ntrackBPC()>0 ) h1_eff-> Fill(9);
  }

  {
    bool FDC1_all_1hit = true;
    for( int lay=1; lay<=6; lay++ ){
      if( blMan->nFDC1(lay)!=1 ) FDC1_all_1hit = false;
    }

    if( blMan-> nFDC1(1)==1 && blMan-> nFDC1(6)==1 ){
      int wire1 = blMan->FDC1(1, 0)-> wire();
      int wire2 = blMan->FDC1(6, 0)-> wire();
      if( 1<wire1 && wire1<32 && 1<wire2 && wire2<32 ) FDC1_trig_16 = true;
    }
    TH1 *h1_eff = getTH("FDC1_eff_18");
    h1_eff-> Fill(0);
    for( int lay=1; lay<=6; lay++ ){
      h1 = getTH(Form("mul_FDC1_%d", lay)), h1-> Fill(blMan->nFDC1(lay));
      for( int i=0; i<blMan->nFDC1(lay); i++ ){
	ChamberLikeHit *hit = blMan-> FDC1(lay, i);
	int tdc = hit-> tdc(), wire = hit-> wire();
	double dt = hit-> dt();
	double dl = hit-> dl();
	
	h1 = getTH(Form("hitpat_FDC1_%d", lay)), h1-> Fill(wire);
	h1 = getTH(Form("FDC1_tdc_%d", lay)), h1-> Fill(tdc);
	h1 = getTH(Form("FDC1_dt_%d", lay)), h1-> Fill(dt);
	h1 = getTH(Form("FDC1_dl_%d", lay)), h1-> Fill(dl);

	h1 = getTH(Form("FDC1_tdc_%d_%d", lay, wire)), h1-> Fill(tdc);
	h1 = getTH(Form("FDC1_dt_%d_%d", lay, wire)), h1-> Fill(dt);
	h1 = getTH(Form("FDC1_dl_%d_%d", lay, wire)), h1-> Fill(dl);
	if( FDC1_all_1hit ){
	  h1 = getTH(Form("FDC1_tdc_%d_1hit", lay)), h1-> Fill(tdc);
	  h1 = getTH(Form("FDC1_dt_%d_1hit", lay)), h1-> Fill(dt);
	  h1 = getTH(Form("FDC1_dl_%d_1hit", lay)), h1-> Fill(dl);

	  h1 = getTH(Form("FDC1_tdc_%d_%d_1hit", lay, wire)), h1-> Fill(tdc);
	  h1 = getTH(Form("FDC1_dt_%d_%d_1hit", lay, wire)), h1-> Fill(dt);
	  h1 = getTH(Form("FDC1_dl_%d_%d_1hit", lay, wire)), h1-> Fill(dl);
	}
      }
      if( FDC1_trig_16 ) if( 0<blMan->nFDC1(lay) && blMan->nFDC1(lay)<3 ) h1_eff->Fill(lay);
    }
    if( bltrackMan-> ntrackFDC1()==1 ){
      LocalTrack *track = bltrackMan->trackFDC1(0);
      for( int i=0; i<track->nhit(); i++ ){
	ChamberLikeHit *hit = track->hit(i);
	int layer = hit-> layer(), wire = hit->wire();
	double dltrack;
	if( hit->xy()==0 )      dltrack = hit->x()-hit->wx();
	else if( hit->xy()==1 ) dltrack = hit->y()-hit->wy();

	h1 = getTH(Form("FDC1_res_%d", layer)), h1-> Fill(hit->resl());
	h1 = getTH(Form("FDC1_res_%d_%d", layer, wire)), h1-> Fill(hit->resl());
	if( FDC1_all_1hit ){
	  h1 = getTH(Form("FDC1_res_%d_1hit", layer)), h1-> Fill(hit->resl());
	  h1 = getTH(Form("FDC1_res_%d_%d_1hit", layer, wire)), h1-> Fill(hit->resl());
	}
      }
    }

    if( FDC1_trig_16 ) if( bltrackMan->ntrackFDC1()>0 ) h1_eff-> Fill(7);
  }

  if( bltrackMan->ntrackBLC2()==1 && bltrackMan->ntrackBPC()==1 ){
    LocalTrack *BPCtrack = bltrackMan->trackBPC(0);
    LocalTrack *BLC2track = bltrackMan->trackBLC2(0);
    double center_pos = -55;
    TVector3 pos_byBPC =  BPCtrack->GetPosatZ(center_pos);
    TVector3 pos_byBLC2 = BLC2track-> GetPosatZ(center_pos);
    TVector3 pos_diff = pos_byBLC2-pos_byBPC;
    TVector3 dir_byBPC  = BPCtrack-> GetMomDir();
    dir_byBPC *= 1./dir_byBPC.Z();
    TVector3 dir_byBLC2 = BLC2track-> GetMomDir();
    dir_byBLC2 *= 1./dir_byBLC2.Z();
    TVector3 dir_diff = dir_byBLC2-dir_byBPC;

    h2 = getTH2("BLC2BPC_pos_diff"), h2-> Fill(pos_diff.X(), pos_diff.Y());
    h2 = getTH2("BLC2BPC_dir_diff"), h2-> Fill(dir_diff.X(), dir_diff.Y());
  }

  if( bltrackMan->ntrackBLC1a()==1 ){
    LocalTrack *BLC1a_track = bltrackMan-> trackBLC1a(0);
    double x, y;
    BLC1a_track-> XYLocalPosatZ(15.0, x, y);
    double r = sqrt(x*x+y*y);
    if( r<7.5 ){
      h1 = getTH("BLC1b_eff"), h1-> Fill(0);
      if( bltrackMan->ntrackBLC1b()>0 ) h1-> Fill(1);

      h1 = getTH("BLC1b_lay_eff"), h1-> Fill(0);
      for( int lay=1; lay<=8; lay++ ){
	if( 0<blMan->nBLC1b(lay) && blMan->nBLC1b(lay)<5 ) h1-> Fill(lay);
      }
    }    
  }

  if( bltrackMan->ntrackBLC1b()==1 ){
    LocalTrack *BLC1b_track = bltrackMan-> trackBLC1b(0);
    double x, y;
    BLC1b_track-> XYLocalPosatZ(-15.0, x, y);
    double r = sqrt(x*x+y*y);
    if( r<7.5 ){
      h1 = getTH("BLC1a_eff"), h1-> Fill(0);
      if( bltrackMan->ntrackBLC1a()>0 ) h1-> Fill(1);

      h1 = getTH("BLC1a_lay_eff"), h1-> Fill(0);
      for( int lay=1; lay<=8; lay++ ){
	if( 0<blMan->nBLC1a(lay) && blMan->nBLC1a(lay)<5 ) h1-> Fill(lay);
      }
    }    
  }

  if( bltrackMan->ntrackBLC2a()==1 ){
    LocalTrack *BLC2a_track = bltrackMan-> trackBLC2a(0);
    double x, y;
    BLC2a_track-> XYLocalPosatZ(15.0, x, y);
    double r = sqrt(x*x+y*y);
    if( r<5 ){
      h1 = getTH("BLC2b_eff"), h1-> Fill(0);
      if( bltrackMan->ntrackBLC2b()>0 ) h1-> Fill(1);

      h1 = getTH("BLC2b_lay_eff"), h1-> Fill(0);
      for( int lay=1; lay<=8; lay++ ){
	if( 0<blMan->nBLC2b(lay) && blMan->nBLC2b(lay)<5 ) h1-> Fill(lay);
      }
    }    
  }

  if( bltrackMan->ntrackBLC2b()==1 ){
    LocalTrack *BLC2b_track = bltrackMan-> trackBLC2b(0);
    double x, y;
    BLC2b_track-> XYLocalPosatZ(-15.0, x, y);
    double r = sqrt(x*x+y*y);
    if( r<5 ){
      h1 = getTH("BLC2a_eff"), h1-> Fill(0);
      if( bltrackMan->ntrackBLC2a()>0 ) h1-> Fill(1);

      h1 = getTH("BLC2a_lay_eff"), h1-> Fill(0);
      for( int lay=1; lay<=8; lay++ ){
	if( 0<blMan->nBLC2a(lay) && blMan->nBLC2a(lay)<5 ) h1-> Fill(lay);
      }
    }    
  }

  if( bltrackMan->ntrackBLC2()==1 ){
    LocalTrack *BLC2_track = bltrackMan-> trackBLC2(0);
    double x, y;
    BLC2_track-> XYPosatZ(-19.70, x, y);
    double r = sqrt(x*x+y*y);
    h2 = getTH2("prof_atBPC_byBLC2"), h2-> Fill(x, y);
    if( r<4 ){
      h1 = getTH("BPC_eff"), h1-> Fill(0);
      if( bltrackMan->ntrackBPC()>0 ) h1-> Fill(1);
    }
  }

  if( fBPD_hit.size()==1 && fDEF_hit.size()==1 ){
    int BPDseg = fBPD_hit[0]-> seg();
    int DEFseg = fDEF_hit[0]-> seg();

    h1 = getTH("BPC_eff"), h1-> Fill(2);
    if( bltrackMan->ntrackBPC()>0 ) h1-> Fill(3);

    if( 30<BPDseg && BPDseg<41 ){
      h1 = getTH("BPC_eff"), h1-> Fill(4);
      if( bltrackMan->ntrackBPC()>0 ) h1-> Fill(5);
      if( 1<DEFseg && DEFseg<8 ){
	h1 = getTH("BPC_eff"), h1-> Fill(6);
	if( bltrackMan->ntrackBPC()>0 ) h1-> Fill(7);
      }
    }
  }

  if( fBVC_hit.size()==1 && fCVC_hit.size()==1 ){
    int CVCseg = fCVC_hit[0]-> seg();
    int BVCseg = fBVC_hit[0]-> seg();

    if( 14<CVCseg && CVCseg<25 && 1<BVCseg && BVCseg<8 ){
      h1 = getTH("FDC1_eff"), h1-> Fill(0);
      if( bltrackMan->ntrackFDC1()>0 ){
	h1-> Fill(1);
      }
    }
  }
  if( bltrackMan->ntrackBLC2()==1 && bltrackMan->ntrackBPC()==1 ){
    LocalTrack *trackBLC2 = bltrackMan->trackBLC2(0);
    LocalTrack *trackBPC = bltrackMan->trackBPC(0);
    TVector3 BLC2pos = trackBLC2-> GetPosatZ(-130);
    TVector3 BPCpos = trackBPC-> GetPosatZ(-20);
    TVector3 dir = (BPCpos-BLC2pos).Unit();
    TVector3 FDC1pos = BPCpos+(170-BPCpos.Z())*dir;
    double  r = sqrt(FDC1pos.X()*FDC1pos.X()+FDC1pos.Y()*FDC1pos.Y());
    if( r<4 ){
      h1 = getTH("FDC1_eff"), h1-> Fill(2);
      if( bltrackMan->ntrackFDC1()>0 ){
	h1-> Fill(3);
      }
      if( fBVC_hit.size()==1 ){
	int BVCseg = fBVC_hit[0]->seg();
	if( 1<BVCseg && BVCseg<8 ){
	  h1-> Fill(4);
	  if( bltrackMan->ntrackFDC1()>0 ){
	    h1-> Fill(5);
	  }
	}
	h1 = getTH("FDC1_lay_eff"), h1->Fill(0);
	for( int lay=1; lay<=6; lay++ ){
	  if( 0<blMan->nFDC1(lay) && blMan->nFDC1(lay)<5 ) h1-> Fill(lay);
	}
      }
    }
  }

  if( header-> IsTrig(Trig_Kf) ){
    if( fDEF_hit.size()==1 && bltrackMan->ntrackBLC2()==1 ){
      LocalTrack *BLC2_track = bltrackMan-> trackBLC2(0);
      double x, y;
      BLC2_track-> XYPosatZ(-19.70, x, y);
      double r = sqrt(x*x+y*y);
      
      if( fDEF_hit[0]->emean()>0.2 ){
	if( r<4 ){
	  h1 = getTH("BPC_eff"), h1-> Fill(8);
	  if( bltrackMan->ntrackBPC()>0 ){
	    h1-> Fill(9);
	  }

	  h1 = getTH("BPC_lay_eff"), h1-> Fill(0);
	  for( int lay=1; lay<=8; lay++ ){
	    if( 0<blMan->nBPC(lay) && blMan->nBPC(lay)<5 ) h1-> Fill(lay);
	  }
	}
      }
    }
  }



  if( bltrackMan->ntrackBPC()==1 ){
    LocalTrack *trackBPC = bltrackMan->trackBPC(0);
    TVector3 prof_FF = trackBPC-> GetPosatZ(0);
    h2 = getTH2("profFF_byBPC"), h2-> Fill(prof_FF.X(), prof_FF.Y());
    if( header-> trigmode2(Mode_KCDH2) ) h2 = getTH2("profFF_byBPC_CDH2"), h2-> Fill(prof_FF.X(), prof_FF.Y());
    if( header-> trigmode2(Mode_Kf) ) h2 = getTH2("profFF_byBPC_K"), h2-> Fill(prof_FF.X(), prof_FF.Y());
  }
}
