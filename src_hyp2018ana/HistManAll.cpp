#include "HistManAll.h"

// very rough
static const double BEAM_PI_MIN = 24.5; 
static const double BEAM_PI_MAX = 27.25;
static const double BEAM_K_MIN = 27.25; 
static const double BEAM_K_MAX = 30.0;

HistManAll::HistManAll(TFile *f, ConfMan *conf) : fFile(f)
{
  fD5 = new BeamSpectrometer(conf);
  clear();

  fFile-> cd();
  std::cout<<"Counter Histogram..."<<std::endl;
  new TH1F("trig_1",  "Beam/f Trigger Time",       2000, -100, 100);
  new TH1F("trig_2",  "Kaon Trigger Time",         2000, -100, 100);
  new TH1F("trig_3",  "KCDH1/f Trigger Time",      2000, -100, 100);
  new TH1F("trig_4",  "Pion Time",                 2000, -100, 100);
  new TH1F("trig_5",  "Proton Trigger Time",       2000, -100, 100);
  new TH1F("trig_6",  "KCDH1 Trigger Time",        2000, -100, 100);
  new TH1F("trig_7",  "KCDH2 Trigger Time",        2000, -100, 100);
  new TH1F("trig_8",  "#pi#bar{BVC} Trigger Time", 2000, -100, 100);
  new TH1F("trig_9",  "#piCDH1 Trigger Time",      2000, -100, 100);
  new TH1F("trig_10" , "#piCDH2 Trigger Time",     2000, -100, 100);
  new TH1F("trig_11", "Kaon/f Trigger Time",       2000, -100, 100);
  new TH1F("trig_12", "1st Mix Trigger Time",      2000, -100, 100);
  new TH1F("trig_13", "Charged Trigger Time",      2000, -100, 100);
  new TH1F("trig_14", "Neutral Trigger Time",      2000, -100, 100);
  new TH1F("trig_15", "Cosmic Trigger Time",       2000, -100, 100);
  new TH1F("trig_16", "K#bar{BVC} Trigger Time",   2000, -100, 100);

  std::cout<<"  for T0"<<std::endl;
  new TH1F("hitpat_T0", "T0 Hit Pattern",  5,  0.5, 5.5);
  new TH1F("mul_T0",    "T0 Multiplicity", 6, -0.5, 5.5);
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("T0_ADC_up_%d", seg),   Form("T0 seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("T0_ADC_down_%d", seg), Form("T0 seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("T0_AwT_up_%d", seg),   Form("T0 seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("T0_AwT_down_%d", seg), Form("T0 seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("T0_dE_up_%d", seg), Form("T0 seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("T0_dE_down_%d", seg), Form("T0 seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("T0_dEwT_up_%d", seg), Form("T0 seg%d up dE w/ TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("T0_dEwT_down_%d", seg), Form("T0 seg%d down dE w/TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("T0_dE_%d", seg), Form("T0 dE seg%d", seg), 500, -0.5, 4.5);
    new TH1F(Form("T0_dEwT_%d", seg), Form("T0 dE seg%d", seg), 500, -0.5, 4.5);
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

    new TH1F(Form("DEF_dE_up_%d", seg), Form("DEF seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("DEF_dE_down_%d", seg), Form("DEF seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("DEF_dEwT_up_%d", seg), Form("DEF seg%d up dE w/ TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("DEF_dEwT_down_%d", seg), Form("DEF seg%d down dE w/TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("DEF_dE_%d", seg), Form("DEF dE seg%d", seg), 500, -0.5, 4.5);
    new TH1F(Form("DEF_dEwT_%d", seg), Form("DEF dE seg%d", seg), 500, -0.5, 4.5);

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

    new TH1F(Form("BHD_dE_up_%d", seg), Form("BHD seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("BHD_dE_down_%d", seg), Form("BHD seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("BHD_dEwT_up_%d", seg), Form("BHD seg%d up dE w/ TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("BHD_dEwT_down_%d", seg), Form("BHD seg%d down dE w/TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("BHD_dE_%d", seg), Form("BHD dE seg%d", seg), 500, -0.5, 4.5);
    new TH1F(Form("BHD_dEwT_%d", seg), Form("BHD dE seg%d", seg), 500, -0.5, 4.5);

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

    new TH1F(Form("BPD_dE_up_%d", seg), Form("BPD seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("BPD_dE_down_%d", seg), Form("BPD seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("BPD_dEwT_up_%d", seg), Form("BPD seg%d up dE w/ TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("BPD_dEwT_down_%d", seg), Form("BPD seg%d down dE w/TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("BPD_dE_%d", seg), Form("BPD dE seg%d", seg), 500, -0.5, 4.5);
    new TH1F(Form("BPD_dEwT_%d", seg), Form("BPD dE seg%d", seg), 500, -0.5, 4.5);

    new TH1F(Form("BPD_TDC_up_%d", seg),   Form("BPD seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("BPD_TDC_down_%d", seg), Form("BPD seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("BPD_time_%d", seg), Form("BPD ctime seg%d", seg), 1000, -25, 25);
  }

  std::cout<<"  for CVC"<<std::endl;
  new TH1F("hitpat_CVC", "CVC Hit Pattern",  34,  0.5, 34.5);
  new TH1F("mul_CVC",    "CVC Multiplicity", 35, -0.5, 34.5);
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("CVC_ADC_up_%d", seg),   Form("CVC seg%d up ADC", seg),   4000, 0, 4000);
    new TH1F(Form("CVC_ADC_down_%d", seg), Form("CVC seg%d down ADC", seg), 4000, 0, 4000);
    new TH1F(Form("CVC_AwT_up_%d", seg),   Form("CVC seg%d up ADC w/ TDC", seg),   4000, 0, 4000);
    new TH1F(Form("CVC_AwT_down_%d", seg), Form("CVC seg%d down ADC w/ TDC", seg), 4000, 0, 4000);

    new TH1F(Form("CVC_dE_up_%d", seg), Form("CVC seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("CVC_dE_down_%d", seg), Form("CVC seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("CVC_dEwT_up_%d", seg), Form("CVC seg%d up dE w/ TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("CVC_dEwT_down_%d", seg), Form("CVC seg%d down dE w/TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("CVC_dE_%d", seg), Form("CVC dE seg%d", seg), 500, -0.5, 4.5);
    new TH1F(Form("CVC_dEwT_%d", seg), Form("CVC dE seg%d", seg), 500, -0.5, 4.5);

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

    new TH1F(Form("PC_dE_up_%d", seg), Form("PC seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("PC_dE_down_%d", seg), Form("PC seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("PC_dEwT_up_%d", seg), Form("PC seg%d up dE w/ TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("PC_dEwT_down_%d", seg), Form("PC seg%d down dE w/TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("PC_dE_%d", seg), Form("PC dE seg%d", seg), 500, -0.5, 4.5);
    new TH1F(Form("PC_dEwT_%d", seg), Form("PC dE seg%d", seg), 500, -0.5, 4.5);

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

    new TH1F(Form("NC_dE_up_%d", seg), Form("NC seg%d up dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("NC_dE_down_%d", seg), Form("NC seg%d down dE", seg), 500, -0.5, 4.5);
    new TH1F(Form("NC_dEwT_up_%d", seg), Form("NC seg%d up dE w/ TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("NC_dEwT_down_%d", seg), Form("NC seg%d down dE w/TDC", seg), 500, -0.5, 4.5);
    new TH1F(Form("NC_dE_%d", seg), Form("NC dE seg%d", seg), 500, -0.5, 4.5);
    new TH1F(Form("NC_dEwT_%d", seg), Form("NC dE seg%d", seg), 500, -0.5, 4.5);

    new TH1F(Form("NC_TDC_up_%d", seg),   Form("NC seg%d up TDC", seg),   4000, 0, 4000);
    new TH1F(Form("NC_TDC_down_%d", seg), Form("NC seg%d down TDC", seg), 4000, 0, 4000);
    new TH1F(Form("NC_time_%d", seg), Form("NC ctime seg%d", seg), 1000, -25, 25);
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

      new TH2F(Form("BPDudE_T0%dBPD%d_time",      seg1, seg2), Form("BPD up dE vs T0 seg%d-BPC seg%d TOF",   seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
      new TH2F(Form("BPDddE_T0%dBPD%d_time",      seg1, seg2), Form("BPD down dE vs T0 seg%d-BPC seg%d TOF", seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
      new TH2F(Form("BPDudE_T0%dBPD%d_corrtime",  seg1, seg2), Form("BPD up dE vs T0 seg%d-BPC seg%d TOF",   seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
      new TH2F(Form("BPDddE_T0%dBPD%d_corrtime",  seg1, seg2), Form("BPD down dE vs T0 seg%d-BPC seg%d TOF", seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
      new TH2F(Form("BPDudE_T0%dBPD%d_TOFwoCorr", seg1, seg2), Form("BPD up dE vs T0 seg%d-BPC seg%d TOF",   seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
      new TH2F(Form("BPDddE_T0%dBPD%d_TOFwoCorr", seg1, seg2), Form("BPD down dE vs T0 seg%d-BPC seg%d TOF", seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
      new TH2F(Form("BPDudE_T0%dBPD%d_TOF",       seg1, seg2), Form("BPD up dE vs T0 seg%d-BPC seg%d TOF",   seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
      new TH2F(Form("BPDddE_T0%dBPD%d_TOF",       seg1, seg2), Form("BPD down dE vs T0 seg%d-BPC seg%d TOF", seg1, seg2), 100, -0.5, 4.5, 100, -50, 50);
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

  std::cout<<"Counter Histogram... finish"<<std::endl;

  initBLDC();
  initFC();
}

HistManAll::~HistManAll()
{
  if( fD5 ) delete fD5;
}

void HistManAll::fill(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  TH1 *h1;
  TH2 *h2;
  for( int i=1; i<=16; i++ ){
    h1 = getTH(Form("trig_%d", i)), h1-> Fill(header->time(i));
  }

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

      h1 = getTH(Form("NC_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("NC_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("NC_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("NC_dEwT_down_%d", seg)), h1-> Fill(ed);

      h1 = getTH(Form("NC_TDC_up_%d", seg)), h1-> Fill(tdcu);
      h1 = getTH(Form("NC_TDC_down_%d", seg)), h1-> Fill(tdcd);

      h1 = getTH(Form("NC_time_%d", seg)), h1-> Fill(ctmean);
      h1 = getTH(Form("NC_dEwT_%d", seg)), h1-> Fill(dE);
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

  fillBLDC(header, blMan, bltrackMan);
  fillFC(conf, blMan, bltrackMan);
}

TH1* HistManAll::getTH(const char *name)
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

TH2* HistManAll::getTH2(const char *name)
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

void HistManAll::clear()
{
  fBeamPID = Beam_Other;
  fBHD_hit.clear();
  fT0_hit.clear();
  fBPD_hit.clear();
  fDEF_hit.clear();
  fCVC_hit.clear();
  fPC_hit.clear();
  for( int i=0; i<8; i++ ) fNC_hit[i].clear();
}
