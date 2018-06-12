#include "HistManBeamAna.h"

// very rough
static const double BEAM_PI_MIN = 24.5; 
static const double BEAM_PI_MAX = 27.25;
static const double BEAM_K_MIN = 27.25; 
static const double BEAM_K_MAX = 30.0;

HistManBeamAna::HistManBeamAna(TFile *f, ConfMan *conf) : fFile(f)
{
  fNumNEv = 0;
  fGoodNEv=0;

  fD5 = new BeamSpectrometer(conf);

  fFile-> cd();
  std::cout<<"Counter Histogram..."<<std::endl;
  new TH1F("trig_mode", "Trigger Mode", 17, 0.5, 17.5);
  new TH1F("Kf_Reduction", "K/f Reduction", 10, -0.5, 9.5);

  new TH1F("trig_1",  "Beam/f Trigger Time",       4096, -0.5, 4095.5);
  new TH1F("trig_2",  "Kaon Trigger Time",         4096, -0.5, 4095.5);
  new TH1F("trig_3",  "KCDH1/f Trigger Time",      4096, -0.5, 4095.5);
  new TH1F("trig_4",  "Pion Time",                 4096, -0.5, 4095.5);
  new TH1F("trig_5",  "Proton Trigger Time",       4096, -0.5, 4095.5);
  new TH1F("trig_6",  "KCDH1 Trigger Time",        4096, -0.5, 4095.5);
  new TH1F("trig_7",  "KCDH2 Trigger Time",        4096, -0.5, 4095.5);
  new TH1F("trig_8",  "#pi#bar{BVC} Trigger Time", 4096, -0.5, 4095.5);
  new TH1F("trig_9",  "#piCDH1 Trigger Time",      4096, -0.5, 4095.5);
  new TH1F("trig_10" , "#piCDH2 Trigger Time",     4096, -0.5, 4095.5);
  new TH1F("trig_11", "Kaon/f Trigger Time",       4096, -0.5, 4095.5);
  new TH1F("trig_12", "1st Mix Trigger Time",      4096, -0.5, 4095.5);
  new TH1F("trig_13", "Charged Trigger Time",      4096, -0.5, 4095.5);
  new TH1F("trig_14", "Neutral Trigger Time",      4096, -0.5, 4095.5);
  new TH1F("trig_15", "Cosmic Trigger Time",       4096, -0.5, 4095.5);
  new TH1F("trig_16", "K#bar{BVC} Trigger Time",   4096, -0.5, 4095.5);

  std::cout<<"  for T0"<<std::endl;
  new TH1F("hitpat_T0", "T0 Hit Pattern",  5,  0.5, 5.5);
  new TH1F("mul_T0",    "T0 Multiplicity", 6, -0.5, 5.5);
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("T0_ADC_up_%d", seg),   Form("T0 seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("T0_ADC_down_%d", seg), Form("T0 seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("T0_AwT_up_%d", seg),   Form("T0 seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("T0_AwT_down_%d", seg), Form("T0 seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("T0_dE_up_%d", seg), Form("T0 seg%d up dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("T0_dE_down_%d", seg), Form("T0 seg%d down dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("T0_dEwT_up_%d", seg), Form("T0 seg%d up dE w/ TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("T0_dEwT_down_%d", seg), Form("T0 seg%d down dE w/TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("T0_dE_%d", seg), Form("T0 dE seg%d", seg), 500, -0.5, 9.5);
    new TH1F(Form("T0_dEwT_%d", seg), Form("T0 dE seg%d", seg), 500, -0.5, 9.5);
    new TH1F(Form("T0_time_%d", seg), Form("T0 ctime seg%d", seg), 1000, -25, 25);

    new TH1F(Form("T0_TDC_up_%d", seg),   Form("T0 seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("T0_TDC_down_%d", seg), Form("T0 seg%d down TDC", seg), 4000, 0, 4000);
  }

  std::cout<<"  for DEF"<<std::endl;
  new TH1F("hitpat_DEF", "DEF Hit Pattern",  8,  0.5, 8.5);
  new TH1F("mul_DEF",    "DEF Multiplicity", 8, -0.5, 8.5);
  for( int seg=1; seg<=8; seg++ ){
    new TH1F(Form("DEF_ADC_up_%d", seg),   Form("DEF seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("DEF_ADC_down_%d", seg), Form("DEF seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("DEF_AwT_up_%d", seg),   Form("DEF seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("DEF_AwT_down_%d", seg), Form("DEF seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("DEF_dE_up_%d", seg), Form("DEF seg%d up dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("DEF_dE_down_%d", seg), Form("DEF seg%d down dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("DEF_dEwT_up_%d", seg), Form("DEF seg%d up dE w/ TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("DEF_dEwT_down_%d", seg), Form("DEF seg%d down dE w/TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("DEF_dE_%d", seg), Form("DEF dE seg%d", seg), 500, -0.5, 9.5);
    new TH1F(Form("DEF_dEwT_%d", seg), Form("DEF dE seg%d", seg), 500, -0.5, 9.5);

    new TH1F(Form("DEF_TDC_up_%d", seg),   Form("DEF seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("DEF_TDC_down_%d", seg), Form("DEF seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("DEF_time_%d", seg), Form("DEF ctime seg%d", seg), 1000, -25, 25);
  }

  std::cout<<"  for BHD"<<std::endl;
  new TH1F("hitpat_BHD", "BHD Hit Pattern",  20,  0.5, 20.5);
  new TH1F("mul_BHD",    "BHD Multiplicity", 21, -0.5, 20.5);
  for( int seg=1; seg<=20; seg++ ){
    new TH1F(Form("BHD_ADC_up_%d", seg),   Form("BHD seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("BHD_ADC_down_%d", seg), Form("BHD seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("BHD_AwT_up_%d", seg),   Form("BHD seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BHD_AwT_down_%d", seg), Form("BHD seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("BHD_dE_up_%d", seg), Form("BHD seg%d up dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("BHD_dE_down_%d", seg), Form("BHD seg%d down dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("BHD_dEwT_up_%d", seg), Form("BHD seg%d up dE w/ TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("BHD_dEwT_down_%d", seg), Form("BHD seg%d down dE w/TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("BHD_dE_%d", seg), Form("BHD dE seg%d", seg), 500, -0.5, 9.5);
    new TH1F(Form("BHD_dEwT_%d", seg), Form("BHD dE seg%d", seg), 500, -0.5, 9.5);

    new TH1F(Form("BHD_TDC_up_%d", seg),   Form("BHD seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BHD_TDC_down_%d", seg), Form("BHD seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("BHD_time_%d", seg), Form("BHD ctime seg%d", seg), 1000, -25, 25);
  }

  std::cout<<"  for BPD"<<std::endl;
  new TH1F("hitpat_BPD", "BPD Hit Pattern", 70, 0.5, 70.5);
  new TH1F("mul_BPD", "BPD Multiplicity", 71, -0.5, 70.5);
  for( int seg=1; seg<=70; seg++ ){
    new TH1F(Form("BPD_ADC_up_%d", seg),   Form("BPD seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_ADC_down_%d", seg), Form("BPD seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("BPD_AwT_up_%d", seg),   Form("BPD seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_AwT_down_%d", seg), Form("BPD seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("BPD_dE_up_%d", seg), Form("BPD seg%d up dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("BPD_dE_down_%d", seg), Form("BPD seg%d down dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("BPD_dEwT_up_%d", seg), Form("BPD seg%d up dE w/ TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("BPD_dEwT_down_%d", seg), Form("BPD seg%d down dE w/TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("BPD_dE_%d", seg), Form("BPD dE seg%d", seg), 500, -0.5, 9.5);
    new TH1F(Form("BPD_dEwT_%d", seg), Form("BPD dE seg%d", seg), 500, -0.5, 9.5);

    new TH1F(Form("BPD_TDC_up_%d", seg),   Form("BPD seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_TDC_down_%d", seg), Form("BPD seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("BPD_time_%d", seg), Form("BPD ctime seg%d", seg), 1000, -25, 25);
  }

  std::cout<<"  for BVC"<<std::endl;
  new TH1F("hitpat_BVC", "BVC Hit Pattern", 70, 0.5, 70.5);
  new TH1F("mul_BVC", "BVC Multiplicity", 71, -0.5, 70.5);
  for( int seg=1; seg<=70; seg++ ){
    new TH1F(Form("BVC_ADC_up_%d", seg),   Form("BVC seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("BVC_ADC_down_%d", seg), Form("BVC seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("BVC_AwT_up_%d", seg),   Form("BVC seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BVC_AwT_down_%d", seg), Form("BVC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("BVC_dE_up_%d", seg), Form("BVC seg%d up dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("BVC_dE_down_%d", seg), Form("BVC seg%d down dE", seg), 500, -0.5, 9.5);
    new TH1F(Form("BVC_dEwT_up_%d", seg), Form("BVC seg%d up dE w/ TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("BVC_dEwT_down_%d", seg), Form("BVC seg%d down dE w/TDC", seg), 500, -0.5, 9.5);
    new TH1F(Form("BVC_dE_%d", seg), Form("BVC dE seg%d", seg), 500, -0.5, 9.5);
    new TH1F(Form("BVC_dEwT_%d", seg), Form("BVC dE seg%d", seg), 500, -0.5, 9.5);

    new TH1F(Form("BVC_TDC_up_%d", seg),   Form("BVC seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BVC_TDC_down_%d", seg), Form("BVC seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("BVC_time_%d", seg), Form("BVC ctime seg%d", seg), 1000, -25, 25);
  }

  std::cout<<"  for CVC"<<std::endl;
  new TH1F("hitpat_CVC", "CVC Hit Pattern",  34,  0.5, 34.5);
  new TH1F("mul_CVC",    "CVC Multiplicity", 35, -0.5, 34.5);
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("CVC_ADC_up_%d", seg),   Form("CVC seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("CVC_ADC_down_%d", seg), Form("CVC seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("CVC_AwT_up_%d", seg),   Form("CVC seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("CVC_AwT_down_%d", seg), Form("CVC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("CVC_dE_up_%d", seg), Form("CVC seg%d up dE", seg), 500, -0.5, 14.5);
    new TH1F(Form("CVC_dE_down_%d", seg), Form("CVC seg%d down dE", seg), 500, -0.5, 14.5);
    new TH1F(Form("CVC_dEwT_up_%d", seg), Form("CVC seg%d up dE w/ TDC", seg), 500, -0.5, 14.5);
    new TH1F(Form("CVC_dEwT_down_%d", seg), Form("CVC seg%d down dE w/TDC", seg), 500, -0.5, 14.5);
    new TH1F(Form("CVC_dE_%d", seg), Form("CVC dE seg%d", seg), 500, -0.5, 14.5);
    new TH1F(Form("CVC_dEwT_%d", seg), Form("CVC dE seg%d", seg), 500, -0.5, 14.5);

    new TH1F(Form("CVC_TDC_up_%d", seg),   Form("CVC seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("CVC_TDC_down_%d", seg), Form("CVC seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("CVC_time_%d", seg), Form("CVC ctime seg%d", seg), 1000, -25, 25);
  }

  std::cout<<"  for PC"<<std::endl;
  new TH1F("hitpat_PC", "PC Hit Pattern",  27,  0.5, 27.5);
  new TH1F("mul_PC",    "PC Multiplicity", 28, -0.5, 27.5);
  for( int seg=1; seg<=27; seg++ ){
    new TH1F(Form("PC_ADC_up_%d", seg),   Form("PC seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("PC_ADC_down_%d", seg), Form("PC seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("PC_AwT_up_%d", seg),   Form("PC seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("PC_AwT_down_%d", seg), Form("PC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("PC_dE_up_%d", seg), Form("PC seg%d up dE", seg), 500, -0.5, 14.5);
    new TH1F(Form("PC_dE_down_%d", seg), Form("PC seg%d down dE", seg), 500, -0.5, 14.5);
    new TH1F(Form("PC_dEwT_up_%d", seg), Form("PC seg%d up dE w/ TDC", seg), 500, -0.5, 14.5);
    new TH1F(Form("PC_dEwT_down_%d", seg), Form("PC seg%d down dE w/TDC", seg), 500, -0.5, 14.5);
    new TH1F(Form("PC_dE_%d", seg), Form("PC dE seg%d", seg), 500, -0.5, 14.5);
    new TH1F(Form("PC_dEwT_%d", seg), Form("PC dE seg%d", seg), 500, -0.5, 14.5);

    new TH1F(Form("PC_TDC_up_%d", seg),   Form("PC seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("PC_TDC_down_%d", seg), Form("PC seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("PC_time_%d", seg), Form("PC ctime seg%d", seg), 1000, -25, 25);
  }

  std::cout<<"  for NC"<<std::endl;
  new TH1F("hitpat_NC", "NC Hit Pattern",  112,  0.5, 112.5);
  new TH1F("mul_NC",    "NC Multiplicity", 113, -0.5, 112.5);
  for( int seg=1; seg<=112; seg++ ){
    new TH1F(Form("NC_ADC_up_%d", seg),   Form("NC seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("NC_ADC_down_%d", seg), Form("NC seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("NC_AwT_up_%d", seg),   Form("NC seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("NC_AwT_down_%d", seg), Form("NC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("NC_dE_up_%d", seg), Form("NC seg%d up dE", seg), 500, -0.5, 49.5);
    new TH1F(Form("NC_dE_down_%d", seg), Form("NC seg%d down dE", seg), 500, -0.5, 49.5);
    new TH1F(Form("NC_dEwT_up_%d", seg), Form("NC seg%d up dE w/ TDC", seg), 500, -0.5, 49.5);
    new TH1F(Form("NC_dEwT_down_%d", seg), Form("NC seg%d down dE w/TDC", seg), 500, -0.5, 49.5);
    new TH1F(Form("NC_dE_%d", seg), Form("NC dE seg%d", seg), 500, -0.5, 49.5);
    new TH1F(Form("NC_dEwT_%d", seg), Form("NC dE seg%d", seg), 500, -0.5, 49.5);

    new TH1F(Form("NC_TDC_up_%d", seg),   Form("NC seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("NC_TDC_down_%d", seg), Form("NC seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("NC_time_%d", seg), Form("NC ctime seg%d", seg), 1000, -25, 25);

    new TH1F(Form("NC_tsub_%d", seg),      Form("NC ut-dt seg%d", seg), 1000, -50, 50);
    new TH1F(Form("NC_tsub_%d_wbar", seg), Form("NC ut-dt seg%d", seg), 1000, -50, 50);

    new TH1F(Form("NC_hitpos_%d", seg),      Form("NC hitpos seg%d", seg), 2000, -1000, 1000);
    new TH1F(Form("NC_hitpos_%d_wbar", seg), Form("NC hitpos seg%d", seg), 2000, -1000, 1000);
  }

  //*** for Calibration ***//
  std::cout<<"  for BHD-T0"<<std::endl;
  new TH1F("BHDT0_TOF",      "BHD-T0 TOF", 10000, -50, 50);
  new TH1F("BHDT0_TOF_1hit", "BHD-T0 TOF BHD 1hit", 10000, -50, 50);
  new TH1F("BHDT0_TOF_k",      "BHD-T0 TOF", 10000, -50, 50);
  new TH1F("BHDT0_TOF_k_1hit", "BHD-T0 TOF BHD 1hit", 10000, -50, 50);
  new TH1F("BHDT0_TOF_pi",      "BHD-T0 TOF", 10000, -50, 50);
  new TH1F("BHDT0_TOF_pi_1hit", "BHD-T0 TOF BHD 1hit", 10000, -50, 50);
  for( int seg=1; seg<=5; seg++ ){

  }

  for( int seg1=1; seg1<=5; seg1++ ){
    new TH1F(Form("BHDT0%d_TOF",seg1), Form("BHD-T0 seg%d TOF",seg1), 10000, -50, 50);
    new TH1F(Form("BHDT0%d_TOF_1hit",seg1), Form("BHD-T0 seg%d TOF BHD 1hit",seg1), 10000, -50, 50);

    new TH1F(Form("BHDT0%d_TOF_k",seg1), Form("BHD-T0 seg%d TOF",seg1), 10000, -50, 50);
    new TH1F(Form("BHDT0%d_TOF_k_1hit",seg1), Form("BHD-T0 seg%d TOF BHD 1hit",seg1), 10000, -50, 50);
    new TH1F(Form("BHDT0%d_TOF_pi",seg1), Form("BHD-T0 seg%d TOF",seg1), 10000, -50, 50);
    new TH1F(Form("BHDT0%d_TOF_pi_1hit",seg1), Form("BHD-T0 seg%d TOF BHD 1hit",seg1), 10000, -50, 50);

    for( int seg2=1; seg2<=20; seg2++ ){
      new TH1F(Form("BHD%dT0%d_TOF", seg2, seg1), Form("BHD seg%d-T0 seg%d TOF", seg2, seg1), 10000, -50, 50);
      new TH1F(Form("BHD%dT0%d_TOF_1hit", seg2, seg1), Form("BHD seg%d-T0 seg%d TOF BHD 1hit", seg2, seg1), 10000, -50, 50);
      new TH1F(Form("BHD%dT0%d_TOF_k", seg2, seg1), Form("BHD seg%d-T0 seg%d TOF", seg2, seg1), 10000, -50, 50);
      new TH1F(Form("BHD%dT0%d_TOF_k_1hit", seg2, seg1), Form("BHD seg%d-T0 seg%d TOF BHD 1hit", seg2, seg1), 10000, -50, 50);
      new TH1F(Form("BHD%dT0%d_TOF_pi", seg2, seg1), Form("BHD seg%d-T0 seg%d TOF", seg2, seg1), 10000, -50, 50);
      new TH1F(Form("BHD%dT0%d_TOF_pi_1hit", seg2, seg1), Form("BHD seg%d-T0 seg%d TOF BHD 1hit", seg2, seg1), 10000, -50, 50);
    }
  }
  for( int seg2=1; seg2<=20; seg2++ ){
    new TH1F(Form("BHD%dT0_TOF",seg2), Form("BHD seg%d-T0 TOF",seg2), 10000, -50, 50);
    new TH1F(Form("BHD%dT0_TOF_1hit",seg2), Form("BHD seg%d-T0 TOF BHD 1hit",seg2), 10000, -50, 50);

    new TH1F(Form("BHD%dT0_TOF_k",seg2), Form("BHD seg%d-T0 TOF",seg2), 10000, -50, 50);
    new TH1F(Form("BHD%dT0_TOF_k_1hit",seg2), Form("BHD seg%d-T0 TOF BHD 1hit",seg2), 10000, -50, 50);

    new TH1F(Form("BHD%dT0_TOF_pi",seg2), Form("BHD seg%d-T0 TOF",seg2), 10000, -50, 50);
    new TH1F(Form("BHD%dT0_TOF_pi_1hit",seg2), Form("BHD seg%d-T0 TOF BHD 1hit",seg2), 10000, -50, 50);
  }

  std::cout<<"  for T0-BPD"<<std::endl;
  new TH1F("T0BPD_TOF",         "T0-BPD TOF", 10000, -50, 50);
  new TH1F("T0BPD_TOF_1hit",    "T0-BPD TOF BHD 1hit", 10000, -50, 50);
  new TH1F("T0BPD_TOF_k",       "T0-BPD TOF", 10000, -50, 50);
  new TH1F("T0BPD_TOF_k_1hit",  "T0-BPD TOF BHD 1hit", 10000, -50, 50);
  new TH1F("T0BPD_TOF_pi",      "T0-BPD TOF", 10000, -50, 50);
  new TH1F("T0BPD_TOF_pi_1hit", "T0-BPD TOF BHD 1hit", 10000, -50, 50);

  for( int seg1=1; seg1<=5; seg1++ ){
    new TH1F(Form("T0%dBPD_TOF",seg1), Form("T0 seg%d-BPD TOF",seg1), 10000, -50, 50);
    new TH1F(Form("T0%dBPD_TOF_1hit",seg1), Form("T0 seg%d-BPD TOF BPD 1hit",seg1), 10000, -50, 50);

    new TH1F(Form("T0%dBPD_TOF_k",seg1), Form("T0-BPD seg%d TOF",seg1), 10000, -50, 50);
    new TH1F(Form("T0%dBPD_TOF_k_1hit",seg1), Form("T0-BPD seg%d TOF BPD 1hit",seg1), 10000, -50, 50);
    new TH1F(Form("T0%dBPD_TOF_pi",seg1), Form("T0-BPD seg%d TOF",seg1), 10000, -50, 50);
    new TH1F(Form("T0%dBPD_TOF_pi_1hit",seg1), Form("T0-BPD seg%d TOF BPD 1hit",seg1), 10000, -50, 50);
    for( int seg2=1; seg2<=70; seg2++ ){
      new TH1F(Form("T0%dBPD%d_TOF", seg1, seg2),         Form("T0 seg%d-BPD seg%d TOF", seg1, seg2), 10000, -50, 50);
      new TH1F(Form("T0%dBPD%d_TOF_1hit", seg1, seg2),    Form("T0 seg%d-BPD seg%d TOF BPD 1hit", seg1, seg2), 10000, -50, 50);
      new TH1F(Form("T0%dBPD%d_TOF_k", seg1, seg2),       Form("T0 seg%d-BPD seg%d TOF", seg1, seg2), 10000, -50, 50);
      new TH1F(Form("T0%dBPD%d_TOF_k_1hit", seg1, seg2),  Form("T0 seg%d-BPD seg%d TOF BPD 1hit", seg1, seg2), 10000, -50, 50);
      new TH1F(Form("T0%dBPD%d_TOF_pi", seg1, seg2),      Form("T0 seg%d-BPD seg%d TOF", seg1, seg2), 10000, -50, 50);
      new TH1F(Form("T0%dBPD%d_TOF_pi_1hit", seg1, seg2), Form("T0 seg%d-BPD seg%d TOF BPD 1hit", seg1, seg2), 10000, -50, 50);
    }
  }
  for( int seg2=1; seg2<=70; seg2++ ){
    new TH1F(Form("T0BPD%d_TOF",seg2),      Form("T0-BPD seg%d TOF",seg2), 10000, -50, 50);
    new TH1F(Form("T0BPD%d_TOF_1hit",seg2), Form("T0-BPD seg%d TOF BPD 1hit",seg2), 10000, -50, 50);

    new TH1F(Form("T0BPD%d_TOF_k",seg2),      Form("T0-BPD seg%d TOF",seg2), 10000, -50, 50);
    new TH1F(Form("T0BPD%d_TOF_k_1hit",seg2), Form("T0-BPD seg%d TOF BPD 1hit",seg2), 10000, -50, 50);

    new TH1F(Form("T0BPD%d_TOF_pi",seg2),      Form("T0-BPD seg%d TOF",seg2), 10000, -50, 50);
    new TH1F(Form("T0BPD%d_TOF_pi_1hit",seg2), Form("T0-BPD seg%d TOF BPD 1hit",seg2), 10000, -50, 50);
  }

  std::cout<<"  for T0-Forward"<<std::endl;
  for( int lay=1; lay<=8; lay++ ){
    new TH1F(Form("hitpat_NC_%d", lay), Form("NC layer%d hit pattern", lay), 14, 1, 15);
  }

  new TH1F("T0NC_TOF",    "T0-NC seg%d TOF",          40000, -100, 300);
  new TH1F("T0NC_TOF_pi", "T0-NC seg%d TOF trig #pi", 40000, -100, 300);
  new TH1F("T0NC_TOF_k",  "T0-NC seg%d TOF trig K",   40000, -100, 300);
  new TH2F("hitpat_NC2D", "NC hit pattern", 14, 0.5, 14.5, 8, 0.5, 8.5);
  for( int seg=1; seg<=114; seg++ ){
    new TH1F(Form("T0NC%d_TOF", seg),    Form("T0-NC seg%d TOF", seg),          40000, -100, 300);
    new TH1F(Form("T0NC%d_TOF_pi", seg), Form("T0-NC seg%d TOF trig #pi", seg), 40000, -100, 300);
    new TH1F(Form("T0NC%d_TOF_k", seg),  Form("T0-NC seg%d TOF trig K", seg),   40000, -100, 300);

    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("T0%dNC%d_TOF", seg2, seg),    Form("T0 seg%d-NC seg%d TOF", seg2, seg),          40000, -100, 300);
      new TH1F(Form("T0%dNC%d_TOF_pi", seg2, seg), Form("T0 seg%d-NC seg%d TOF trig #pi", seg2, seg), 40000, -100, 300);
      new TH1F(Form("T0%dNC%d_TOF_k", seg2, seg),  Form("T0 seg%d-NC seg%d TOF trig K", seg2, seg),   40000, -100, 300);
    }
  }

  new TH1F("T0CVC_TOF", "T0-CVC TOF",    40000, -100, 300);
  new TH1F("T0CVC_TOF_pi", "T0-CVC TOF", 40000, -100, 300);
  new TH1F("T0CVC_TOF_k", "T0-CVC TOF",  40000, -100, 300);
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("T0CVC%d_TOF", seg),    Form("T0-CVC seg%d_TOF", seg), 4000, -100, 300);
    new TH1F(Form("T0CVC%d_TOF_pi", seg), Form("T0-CVC seg%d_TOF", seg), 4000, -100, 300);
    new TH1F(Form("T0CVC%d_TOF_k", seg),  Form("T0-CVC seg%d_TOF", seg), 4000, -100, 300);
    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("T0%dCVC%d_TOF", seg2, seg),    Form("T0 seg%d-CVC seg%d_TOF", seg2, seg), 4000, -100, 300);
      new TH1F(Form("T0%dCVC%d_TOF_pi", seg2, seg), Form("T0 seg%d-CVC seg%d_TOF", seg2, seg), 4000, -100, 300);
      new TH1F(Form("T0%dCVC%d_TOF_k", seg2, seg),  Form("T0 seg%d-CVC seg%d_TOF", seg2, seg), 4000, -100, 300);
    }
  }

  new TH1F("T0PC_TOF", "T0-PC TOF",    40000, -100, 300);
  new TH1F("T0PC_TOF_pi", "T0-PC TOF", 40000, -100, 300);
  new TH1F("T0PC_TOF_k", "T0-PC TOF",  40000, -100, 300);
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("T0PC%d_TOF", seg),    Form("T0-PC seg%d_TOF", seg), 4000, -100, 300);
    new TH1F(Form("T0PC%d_TOF_pi", seg), Form("T0-PC seg%d_TOF", seg), 4000, -100, 300);
    new TH1F(Form("T0PC%d_TOF_k", seg),  Form("T0-PC seg%d_TOF", seg), 4000, -100, 300);
    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("T0%dPC%d_TOF", seg2, seg),    Form("T0 seg%d-PC seg%d_TOF", seg2, seg), 4000, -100, 300);
      new TH1F(Form("T0%dPC%d_TOF_pi", seg2, seg), Form("T0 seg%d-PC seg%d_TOF", seg2, seg), 4000, -100, 300);
      new TH1F(Form("T0%dPC%d_TOF_k", seg2, seg),  Form("T0 seg%d-PC seg%d_TOF", seg2, seg), 4000, -100, 300);
    }
  }

  std::cout<<"  for T0 stable"<<std::endl;
  for( int seg=1; seg<=5; seg++ ){
    new TH2F(Form("evnum_BHDT0%d_TOF",seg), Form("event number vs BHD-T0 seg%d TOF",seg), 2500, 0, 500000, 100, 25, 35);
  }

  for( int seg=1; seg<=20; seg++ ){
    new TH2F(Form("evnum_BHD%dT0_TOF",seg), Form("event number vs BHD seg%d-T0 TOF",seg), 2500, 0, 500000, 100, 25, 35);
  }

  std::cout<<"  for BeamLineSpectrometer"<<std::endl;
  new TH1F("BeamMomByD5", "Beam Mom by D5", 1300, 0.0, 1.3);
  new TH1F("D5_Chi2",     "D5 Chi-Square",  500,    0, 1000);

  std::cout<<"Counter Histogram... finish"<<std::endl;

  new TTree("TreeBHDT0", "TreeBHDT0");
  new TTree("TreeT0NC",  "TreeT0NC");
  new TTree("TreeT0CVC", "TreeT0CVC");
  new TTree("TreeT0PC",  "TreeT0PC");
  new TTree("TreeT0NC_p",  "TreeT0NC_production");
#if SLWEING_DATA
  fFile-> cd();
  fSlewingBHDT0 = new SlewingData(CID_BHD, CID_T0);
  TTree *tree = (TTree*)fFile-> Get("TreeBHDT0");
  tree-> Branch("SlewingBHDT0", &fSlewingBHDT0);

  fSlewingT0NC = new SlewingData(CID_T0, CID_NC);
  tree = (TTree*)fFile-> Get("TreeT0NC");
  tree-> Branch("SlewingT0NC", &fSlewingT0NC);

  fSlewingT0CVC = new SlewingData(CID_T0, CID_CVC);
  tree = (TTree*)fFile-> Get("TreeT0CVC");
  tree-> Branch("SlewingT0CVC", &fSlewingT0CVC);

  fSlewingT0PC = new SlewingData(CID_T0, CID_PC);
  tree = (TTree*)fFile-> Get("TreeT0PC");
  tree-> Branch("SlewingT0PC", &fSlewingT0PC);
#endif
  initBLDC();
  initFC();

  clear();
}

HistManBeamAna::~HistManBeamAna()
{
  if( fD5 ) delete fD5;
}

void HistManBeamAna::fill(ConfMan *conf, EventHeader *header)
{
  TH1 *h1;
  h1 = getTH("trig_mode");
  int trigmode = header-> trigmode(conf);
  if( trigmode==Mode_Beam        )  h1->Fill(1);
  else if( trigmode==Mode_Kf     )  h1->Fill(2);
  else if( trigmode==Mode_KCDH1f )  h1->Fill(3);
  else if( trigmode==Mode_PiN    )  h1->Fill(4);
  else if( trigmode==Mode_PiC    )  h1->Fill(5);
  else if( trigmode==Mode_KCDH1N )  h1->Fill(6);
  else if( trigmode==Mode_KCDH1C )  h1->Fill(7);
  else if( trigmode==Mode_KCDH2  )  h1->Fill(8);
  else if( trigmode==Mode_KvBVCN )  h1->Fill(9);
  else if( trigmode==Mode_KvBVCC )  h1->Fill(10);
  else if( trigmode==Mode_PiCDH1N ) h1->Fill(11);
  else if( trigmode==Mode_PiCDH1C ) h1->Fill(12);
  else if( trigmode==Mode_Cosmic  ) h1->Fill(13);
  else if( trigmode==Mode_Reject  ) h1->Fill(14);
  else if( trigmode==Mode_Unknown ) h1->Fill(15);

  for( int i=1; i<=16; i++ ){
    h1 = getTH(Form("trig_%d", i)), h1-> Fill(header->time(i));
  }
}

void HistManBeamAna::fill(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  TH1 *h1;
  TH2 *h2;

  for( int i=0; i<blMan->nT0(); i++ ){
    HodoscopeLikeHit *hit = blMan-> T0(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("T0_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("T0_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("T0_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("T0_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("T0_dE_%d", seg)), h1->Fill(dE);

    if( hit-> CheckRange() ){
      fT0_hit.push_back(hit);
      h1 = getTH("hitpat_T0"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctmean = hit-> ctmean();

      h1 = getTH(Form("T0_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("T0_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("T0_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("T0_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("T0_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("T0_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("T0_time_%d", seg)), h1->Fill(ctmean);
      h1 = getTH(Form("T0_dEwT_%d", seg)), h1->Fill(dE);
    }
  }
  h1 = getTH("mul_T0"), h1-> Fill(fT0_hit.size());

  for( int i=0; i<blMan->nDEF(); i++ ){
    HodoscopeLikeHit *hit = blMan-> DEF(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("DEF_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("DEF_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("DEF_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("DEF_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("DEF_dE_%d", seg)), h1-> Fill(dE);

    if( hit-> CheckRange() ){
      fDEF_hit.push_back(hit);
      h1 = getTH("hitpat_DEF"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctmean = hit-> ctmean();

      h1 = getTH(Form("DEF_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("DEF_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("DEF_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("DEF_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("DEF_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("DEF_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("DEF_time_%d", seg)), h1-> Fill(ctmean);
      h1 = getTH(Form("DEF_dEwT_%d", seg)), h1-> Fill(dE);
    }
  }
  h1 = getTH("mul_DEF"), h1-> Fill(fDEF_hit.size());

  for( int i=0; i<blMan->nBHD(); i++ ){
    HodoscopeLikeHit *hit = blMan-> BHD(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("BHD_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("BHD_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("BHD_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("BHD_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("BHD_dE_%d", seg)), h1-> Fill(dE);

    if( hit-> CheckRange() ){
      fBHD_hit.push_back(hit);
      h1 = getTH("hitpat_BHD"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctu = hit-> ctu(), ctd=hit-> ctd();
      double ctmean = hit-> ctmean();

      h1 = getTH(Form("BHD_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("BHD_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("BHD_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("BHD_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("BHD_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("BHD_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("BHD_time_%d", seg)), h1->Fill(ctmean);
      h1 = getTH(Form("BHD_dEwT_%d", seg)), h1->Fill(dE);
    }
  }
  h1 = getTH("mul_BHD"), h1-> Fill(fBHD_hit.size());

  for( int i=0; i<blMan->nBPD(); i++ ){
    HodoscopeLikeHit *hit = blMan-> BPD(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("BPD_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("BPD_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("BPD_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("BPD_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("BPD_dE_%d", seg)),h1-> Fill(dE);

    if( hit-> CheckRange() ){
      fBPD_hit.push_back(hit);
      h1 = getTH("hitpat_BPD"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctu = hit-> ctu(), ctd=hit-> ctd();
      double ctmean = hit-> ctmean();

      h1 = getTH(Form("BPD_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("BPD_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("BPD_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("BPD_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("BPD_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("BPD_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("BPD_time_%d", seg)), h1-> Fill(ctmean);
      h1 = getTH(Form("BPD_dEwT_%d", seg)), h1-> Fill(dE);
    }
  }
  h1 = getTH("mul_BPD"), h1-> Fill(fBPD_hit.size());

  for( int i=0; i<blMan->nBVC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> BVC(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("BVC_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("BVC_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("BVC_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("BVC_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("BVC_dE_%d", seg)),h1-> Fill(dE);

    if( hit-> CheckRange() ){
      fBVC_hit.push_back(hit);
      h1 = getTH("hitpat_BVC"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctu = hit-> ctu(), ctd=hit-> ctd();
      double ctmean = hit-> ctmean();

      h1 = getTH(Form("BVC_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("BVC_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("BVC_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("BVC_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("BVC_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("BVC_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("BVC_time_%d", seg)), h1-> Fill(ctmean);
      h1 = getTH(Form("BVC_dEwT_%d", seg)), h1-> Fill(dE);
    }
  }
  h1 = getTH("mul_BVC"), h1-> Fill(fBVC_hit.size());

  for( int i=0; i<blMan->nCVC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> CVC(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("CVC_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("CVC_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("CVC_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("CVC_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("CVC_dE_%d", seg)),h1-> Fill(dE);

    if( hit-> CheckRange() ){
      fCVC_hit.push_back(hit);
      h1 = getTH("hitpat_CVC"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctu = hit-> ctu(), ctd=hit-> ctd();
      double ctmean = hit-> ctmean();

      h1 = getTH(Form("CVC_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("CVC_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("CVC_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("CVC_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("CVC_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("CVC_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("CVC_time_%d", seg)), h1-> Fill(ctmean);
      h1 = getTH(Form("CVC_dEwT_%d", seg)), h1-> Fill(dE);
    }
  }
  h1 = getTH("mul_CVC"), h1-> Fill(fCVC_hit.size());

  for( int i=0; i<blMan->nPC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> PC(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("PC_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("PC_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("PC_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("PC_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("PC_dE_%d", seg)),h1-> Fill(dE);

    if( hit-> CheckRange() ){
      fPC_hit.push_back(hit);
      h1 = getTH("hitpat_PC"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctu = hit-> ctu(), ctd=hit-> ctd();
      double ctmean = hit-> ctmean();

      h1 = getTH(Form("PC_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("PC_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("PC_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("PC_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("PC_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("PC_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("PC_time_%d", seg)), h1-> Fill(ctmean);
      h1 = getTH(Form("PC_dEwT_%d", seg)), h1-> Fill(dE);
    }
  }
  h1 = getTH("mul_PC"), h1-> Fill(fPC_hit.size());

  for( int i=0; i<blMan->nLB(); i++ ){
    HodoscopeLikeHit *hit = blMan->LB(i);
    if( hit-> CheckRange() ) fLB_hit.push_back(hit);
  }

  for( int i=0; i<blMan->nNC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> NC(i);
    int seg = hit-> seg();
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();
    double dE = hit-> emean();

    h1 = getTH(Form("NC_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("NC_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("NC_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("NC_dE_down_%d", seg)), h1-> Fill(ed);
    h1 = getTH(Form("NC_dE_%d", seg)),h1-> Fill(dE);

    if( hit-> CheckRange() ){
      h1 = getTH("hitpat_NC"), h1-> Fill(seg);
      int tdcu = hit-> tdcu(),tdcd = hit-> tdcd();
      double tu = hit-> tu(), td=hit-> td();
      double ctu = hit-> ctu(), ctd=hit-> ctd();
      double ctmean = hit-> ctmean();
      double tsub = hit-> ctsub();
      double hitpos = hit-> hitpos();

      h1 = getTH(Form("NC_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("NC_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("NC_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("NC_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("NC_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("NC_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("NC_time_%d", seg)), h1-> Fill(ctmean);
      h1 = getTH(Form("NC_dEwT_%d", seg)), h1-> Fill(dE);

      h1 = getTH(Form("NC_tsub_%d", seg)), h1-> Fill(tsub);
      h1 = getTH(Form("NC_hitpos_%d", seg)), h1-> Fill(hitpos);
      if( fLB_hit.size()>0 ){
	h1 = getTH(Form("NC_tsub_%d_wbar", seg)), h1-> Fill(tsub);
	h1 = getTH(Form("NC_hitpos_%d_wbar", seg)), h1-> Fill(hitpos);
      }
    }
  }

  for( int i=0; i<blMan->nNC(); i++ ){
    HodoscopeLikeHit *hit = blMan-> NC(i);
    int seg = hit-> seg();
    int lay = 1+(seg-1)/14;
    int seg2 = 1+(seg-1)%14;
    int adcu = hit-> adcu(), adcd = hit-> adcd();
    double eu = hit-> eu(), ed=hit-> ed();

    if( hit-> CheckRange() ){
      fNC_hit[lay-1].push_back(hit);
      h1 = getTH(Form("hitpat_NC_%d", lay)), h1-> Fill(seg2);
      h2 = getTH2("hitpat_NC2D"), h2-> Fill(seg2, lay);
    }
  }

  if( fT0_hit.size()==1 ){
    if( fBHD_hit.size()==1 ){
#if SLWEING_DATA
      fSlewingBHDT0-> Set(fBHD_hit[0], fT0_hit[0]);
      TTree *tree = (TTree*)fFile-> Get("TreeBHDT0");
      tree-> Fill();
#endif
      double tof = fT0_hit[0]-> ctmean()- fBHD_hit[0]-> ctmean();
      int T0seg = fT0_hit[0]->seg(), BHDseg = fBHD_hit[0]->seg();
      h1 = getTH("BHDT0_TOF_1hit"), h1-> Fill(tof);
      h1 = getTH(Form("BHD%dT0_TOF_1hit", BHDseg)), h1-> Fill(tof);
      h1 = getTH(Form("BHDT0%d_TOF_1hit", T0seg)), h1-> Fill(tof);
      h1 = getTH(Form("BHD%dT0%d_TOF_1hit", BHDseg, T0seg)), h1-> Fill(tof);

      h2 = getTH2(Form("evnum_BHD%dT0_TOF", BHDseg)), h2-> Fill(header->ev(), tof);
      h2 = getTH2(Form("evnum_BHDT0%d_TOF", T0seg)), h2-> Fill(header->ev(), tof);
      if( header->kaon(conf) ){
	h1 = getTH("BHDT0_TOF_k_1hit"), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0_TOF_k_1hit", BHDseg)), h1-> Fill(tof);
	h1 = getTH(Form("BHDT0%d_TOF_k_1hit", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0%d_TOF_k_1hit", BHDseg, T0seg)), h1-> Fill(tof);
      }
      if( header->pion(conf) ){
	h1 = getTH("BHDT0_TOF_pi_1hit"), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0_TOF_pi_1hit", BHDseg)), h1-> Fill(tof);
	h1 = getTH(Form("BHDT0%d_TOF_pi_1hit", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0%d_TOF_pi_1hit", BHDseg, T0seg)), h1-> Fill(tof);
      }
    }

    for( int i=0; i<fBHD_hit.size(); i++ ){
      int T0seg = fT0_hit[0]->seg(), BHDseg = fBHD_hit[i]->seg();
      double tof = fT0_hit[0]-> ctmean()- fBHD_hit[i]-> ctmean();
      h1 = getTH("BHDT0_TOF"), h1-> Fill(tof);
      h1 = getTH(Form("BHD%dT0_TOF", BHDseg)), h1-> Fill(tof);
      h1 = getTH(Form("BHDT0%d_TOF", T0seg)), h1-> Fill(tof);
      h1 = getTH(Form("BHD%dT0%d_TOF", BHDseg, T0seg)), h1-> Fill(tof);
      if( header->kaon(conf) ){
	h1 = getTH("BHDT0_TOF_k"), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0_TOF_k", BHDseg)), h1-> Fill(tof);
	h1 = getTH(Form("BHDT0%d_TOF_k", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0%d_TOF_k", BHDseg, T0seg)), h1-> Fill(tof);
      }
      if( header->pion(conf) ){
	h1 = getTH("BHDT0_TOF_pi"), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0_TOF_pi", BHDseg)), h1-> Fill(tof);
	h1 = getTH(Form("BHDT0%d_TOF_pi", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("BHD%dT0%d_TOF_pi", BHDseg, T0seg)), h1-> Fill(tof);
      }

      if( fBeamPID==Beam_Other ){
	if( BEAM_PI_MIN<tof && tof<BEAM_PI_MAX ) fBeamPID=Beam_Pion;
      }
      if( fBeamPID!=Beam_Kaon ){
	if( BEAM_K_MIN<tof && tof<BEAM_K_MAX ) fBeamPID=Beam_Kaon;
      }
    }

    if( fBPD_hit.size()==1 ){
      double tof = fT0_hit[0]-> ctmean()- fBPD_hit[0]-> ctmean();
      int T0seg = fT0_hit[0]->seg(), BPDseg = fBPD_hit[0]->seg();
      h1 = getTH("T0BPD_TOF_1hit"), h1-> Fill(tof);
      h1 = getTH(Form("T0BPD%d_TOF_1hit", BPDseg)), h1-> Fill(tof);
      h1 = getTH(Form("T0%dBPD_TOF_1hit", T0seg)), h1-> Fill(tof);
      h1 = getTH(Form("T0%dBPD%d_TOF_1hit", T0seg, BPDseg)), h1-> Fill(tof);
      if( header->kaon(conf) ){
	h1 = getTH("T0BPD_TOF_k_1hit"), h1-> Fill(tof);
	h1 = getTH(Form("T0BPD%d_TOF_k_1hit", BPDseg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD_TOF_k_1hit", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD%d_TOF_k_1hit", T0seg, BPDseg)), h1-> Fill(tof);
      }
      if( header->pion(conf) ){
	h1 = getTH("T0BPD_TOF_pi_1hit"), h1-> Fill(tof);
	h1 = getTH(Form("T0BPD%d_TOF_pi_1hit", BPDseg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD_TOF_pi_1hit", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD%d_TOF_pi_1hit", T0seg, BPDseg)), h1-> Fill(tof);
      }
    }

    for( int i=0; i<fBPD_hit.size(); i++ ){
      int T0seg = fT0_hit[0]->seg(), BPDseg = fBPD_hit[i]->seg();
      double tof = fT0_hit[0]-> ctmean()- fBPD_hit[i]-> ctmean();
      h1 = getTH("T0BPD_TOF"), h1-> Fill(tof);
      h1 = getTH(Form("T0BPD%d_TOF", BPDseg)), h1-> Fill(tof);
      h1 = getTH(Form("T0%dBPD_TOF", T0seg)), h1-> Fill(tof);
      h1 = getTH(Form("T0%dBPD%d_TOF", T0seg, BPDseg)), h1-> Fill(tof);
      if( header->kaon(conf) ){
	h1 = getTH("T0BPD_TOF_k"), h1-> Fill(tof);
	h1 = getTH(Form("T0BPD%d_TOF_k", BPDseg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD_TOF_k", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD%d_TOF_k", T0seg, BPDseg)), h1-> Fill(tof);
      }
      if( header->pion(conf) ){
	h1 = getTH("T0BPD_TOF_pi"), h1-> Fill(tof);
	h1 = getTH(Form("T0BPD%d_TOF_pi", BPDseg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD_TOF_pi", T0seg)), h1-> Fill(tof);
	h1 = getTH(Form("T0%dBPD%d_TOF_pi", T0seg, BPDseg)), h1-> Fill(tof);
      }
    }

    for( int lay=1; lay<=8; lay++ ){
      for( int i=0; i<fNC_hit[lay-1].size(); i++ ){
	int T0seg = fT0_hit[0]-> seg();
	int NCseg = fNC_hit[lay-1][i]-> seg();
	double TOF = fNC_hit[lay-1][i]-> ctmean()-fT0_hit[0]->ctmean();
	h1= getTH("T0NC_TOF"), h1-> Fill(TOF);
	h1= getTH(Form("T0NC%d_TOF", NCseg)), h1-> Fill(TOF);
	h1= getTH(Form("T0%dNC%d_TOF", T0seg, NCseg)), h1-> Fill(TOF);
	if( header->kaon(conf) ){
	  h1= getTH("T0NC_TOF_k"), h1-> Fill(TOF);
	  h1= getTH(Form("T0NC%d_TOF_k", NCseg)), h1-> Fill(TOF);
	  h1= getTH(Form("T0%dNC%d_TOF_k", T0seg, NCseg)), h1-> Fill(TOF);
	}
	if( header->pion(conf) ){
	  h1= getTH("T0NC_TOF_pi"), h1-> Fill(TOF);
	  h1= getTH(Form("T0NC%d_TOF_pi", NCseg)), h1-> Fill(TOF);
	  h1= getTH(Form("T0%dNC%d_TOF_pi", T0seg, NCseg)), h1-> Fill(TOF);
	}
      }
    }

    for( int i=0; i<fCVC_hit.size(); i++ ){
      double TOF = fCVC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
      int CVCseg= fCVC_hit[0]->seg();
      int T0seg = fT0_hit[0]-> seg();
      h1 = getTH("T0CVC_TOF"), h1-> Fill(TOF);
      h1 = getTH(Form("T0CVC%d_TOF", CVCseg)), h1-> Fill(TOF);
      h1 = getTH(Form("T0%dCVC%d_TOF", T0seg, CVCseg)), h1-> Fill(TOF);
    }

    for( int i=0; i<fPC_hit.size(); i++ ){
      double TOF = fPC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
      int PCseg= fPC_hit[0]->seg();
      int T0seg = fT0_hit[0]-> seg();
      h1 = getTH("T0PC_TOF"), h1-> Fill(TOF);
      h1 = getTH(Form("T0PC%d_TOF", PCseg)), h1-> Fill(TOF);
      h1 = getTH(Form("T0%dPC%d_TOF", T0seg, PCseg)), h1-> Fill(TOF);
    }
  }

  //  fillFC(conf, blMan, bltrackMan);
  //  if( bltrackMan->ntrackBLC1()==1 && bltrackMan->ntrackBLC2()==1 ){
  int ntrackBLC1=0, ntrackBLC2=0;
  LocalTrack *trackBLC1=0;
  LocalTrack *trackBLC2=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    double track_time = bltrackMan->trackBLC1(i)->GetTrackTime();
    if( -10<track_time && track_time<10 ){
      ntrackBLC1++;
      trackBLC1 = bltrackMan-> trackBLC1(i);
    }
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    double track_time = bltrackMan->trackBLC2(i)->GetTrackTime();
    if( -10<track_time && track_time<10 ){
      ntrackBLC2++;
      trackBLC2 = bltrackMan-> trackBLC2(i);
    }
  }

  if( ntrackBLC1==1 && ntrackBLC2==1 ){
    //    fD5-> TMinuitFit(bltrackMan->trackBLC1(0), bltrackMan->trackBLC2(0), conf);
    fD5-> TMinuitFit(trackBLC1, trackBLC2, conf);
    fBeamMom = fD5-> mom();
    h1 = getTH("BeamMomByD5"), h1-> Fill(fD5->mom());
    h1 = getTH("D5_Chi2"), h1-> Fill(fD5->chisquare());
  }
  fillBLDC(header, blMan, bltrackMan);

  if( header->trigmode2(Mode_KCDH1N) ){
    h1 = getTH("Neutral_Reduction"), h1-> Fill(0);
    if( fT0_hit.size()==1 ){
      h1-> Fill(1);
      if( fBeamPID==Beam_Kaon ){
        h1-> Fill(2);
        if( fBeamMom>0.1 ){
          h1-> Fill(3);
          if( bltrackMan->ntrackBPC()==1 ){
            h1->Fill(4);
	  }
	}
      }
    }
  }

  if( header->IsTrig(Trig_Neutral) ){
    h1 = getTH("Neutral_Reduction"), h1-> Fill(0);
    if( fT0_hit.size()==1 ){
      h1-> Fill(1);
      if( fBeamPID==Beam_Kaon ){
        h1-> Fill(2);
        if( fBeamMom>0.1 ){
          h1-> Fill(3);
          if( bltrackMan->ntrackBPC()==1 ){
            h1->Fill(4);
	  }
	}
      }
    }
  }

  if( header->trigmode2(Mode_Kf) ){
    h1 = getTH("Kf_Reduction"), h1-> Fill(0);
    if( fT0_hit.size()==1 ){
      h1-> Fill(1);
      if( fBeamPID==Beam_Kaon ){
        h1-> Fill(2);
        if( fBeamMom>0.1 ){
          h1-> Fill(3);
          if( bltrackMan->ntrackBPC()==1 ){
            h1->Fill(4);
	  }
	}
      }
    }
  }

#if SLWEING_DATA
  if( fT0_hit.size()==1 ){
    if( fBHD_hit.size()==1 ){
      fSlewingBHDT0-> SetBeamPID(fBeamPID);
      fSlewingBHDT0-> Set(fBHD_hit[0], fT0_hit[0]);
      TTree *tree = (TTree*)fFile-> Get("TreeBHDT0");
      tree-> Fill();
    }
  }
#endif
}

TH1* HistManBeamAna::getTH(const char *name)
{
  TH1* h1 = (TH1*)fFile->Get(name);
  if( h1==0 ){
    std::cout<<" !!! getTH("<<name<<") !!!"<<std::endl;
    exit(-1);
  }
  std::string classname = h1-> ClassName();
  if( classname.find("TH")==std::string::npos ){
    std::cout<<"  !!! HistManBPC::getTH1 error name="<<name<<" !!!"<<std::endl;
    exit(-1);
  }
  else return h1;
}

TH2* HistManBeamAna::getTH2(const char *name)
{
  TH2* h2 = (TH2*)fFile->Get(name);
  if( h2==0 ){
    std::cout<<" !!! getTH2("<<name<<") !!!"<<std::endl;
    exit(-1);
  }
  std::string classname = h2-> ClassName();
  if( classname.find("TH2")==std::string::npos ){
    std::cout<<"  !!! HistManBPC::getTH1 error name="<<name<<" !!!"<<std::endl;
    exit(-1);
  }
  else return h2;
}

void HistManBeamAna::clear()
{
  fLflag = false;
  fK0flag = false;
  fVtxDis = DBL_MAX;
  fVtxBeam.SetXYZ(-999, -999, -999);
  fVtxCDS.SetXYZ(-999, -999, -999);
  fVtxBeam_w.SetXYZ(-999, -999, -999);
  fVtxCDS_w.SetXYZ(-999, -999, -999);

  fBeamPID = Beam_Other;
  fBeamMom = 0.0;
  fBHD_hit.clear();
  fT0_hit.clear();
  fBPD_hit.clear();
  fDEF_hit.clear();
  fCVC_hit.clear();
  fBVC_hit.clear();
  fPC_hit.clear();
  for( int i=0; i<8; i++ ) fNC_hit[i].clear();
  fLB_hit.clear();

  fCDH_hit.clear();
  fCDSpim.clear();
  fCDSkm.clear();
  fCDSpip.clear();
  fCDSp.clear();
  fCDSd.clear();

  fD5-> Clear();
#if SLWEING_DATA
  fSlewingBHDT0-> clear();
  fSlewingT0CVC-> clear();
  fSlewingT0PC-> clear();
  fSlewingT0NC-> clear();
#endif
}
