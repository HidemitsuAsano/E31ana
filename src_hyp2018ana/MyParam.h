#ifndef MyParam_h
#define MyParam_h 1

#include <iostream>
#include "TLorentzVector.h"

//*** for CDS PID ***//
static const double PID_Param[3] = { 0.00381414, 0.000119896, 0.0113647 };
static const double Kpi_mid_mass2 = 0.1031;
static const double Ppi_mid_mass2 = 0.76247;

static const TLorentzVector D_LMOM = TLorentzVector(0., 0., 0., dMass);
static const TLorentzVector P_LMOM = TLorentzVector(0., 0., 0., pMass);
static const TLorentzVector N_LMOM = TLorentzVector(0., 0., 0., nMass);
static const double NC_THRE = 8.0;//MeV

//*** for IM parameter ***//
static const double K0_peak = 0.497742;
static const double K0_sigma = 0.00665363;
static const double L_peak = 1.11569;
static const double L_sigma = 0.00196312;
static const double Sp_peak = 1.1887;
static const double Sp_sigma = 0.0041595;
static const double Sm_peak = 1.19683;
static const double Sm_sigma = 0.00460906;

static const double K0_MIN = K0_peak-2.0*K0_sigma;
static const double K0_MAX = K0_peak+2.0*K0_sigma;
static const double K0_MIN3 = K0_peak-3.0*K0_sigma;
static const double K0_MAX3 = K0_peak+3.0*K0_sigma;
static const double K0_MIN4 = K0_peak-4.0*K0_sigma;
static const double K0_MAX4 = K0_peak+4.0*K0_sigma;
static const double K0_MIN5 = K0_peak-5.0*K0_sigma;
static const double K0_MAX5 = K0_peak+5.0*K0_sigma;
static const double L_MIN  = L_peak-2.0*L_sigma;
static const double L_MAX  = L_peak+2.0*L_sigma;
static const double Sp_MIN = Sp_peak-2.0*Sp_sigma;
static const double Sp_MAX = Sp_peak+2.0*Sp_sigma;
static const double Sm_MIN = Sm_peak-2.0*Sm_sigma;
static const double Sm_MAX = Sm_peak+2.0*Sm_sigma;

//*** for p(K-, pi+ pi-) Event ***//
static const double Kpipi_N_MIN = 0.88;
static const double Kpipi_N_MAX = 1.00;
static const double Kpipi_L_peak = 1.11285;
static const double Kpipi_L_sigma = 0.0125546;
static const double Kpipi_L_MIN = Kpipi_L_peak-2.*Kpipi_L_sigma;
static const double Kpipi_L_MAX = Kpipi_L_peak+2.*Kpipi_L_sigma;

//*** for d(K-, n pi+ pi-) Event ***//
static const double Npipi_K0_peak = 0.497689;
static const double Npipi_K0_sigma = 0.00536856;
static const double Npipi_Sp_peak = 1.18901;
static const double Npipi_Sp_sigma = 0.00380568;
static const double Npipi_Sm_peak = 1.19742;
static const double Npipi_Sm_sigma = 0.0048656;

static const double Npipi_K0_MIN = Npipi_K0_peak-2.*Npipi_K0_sigma;
static const double Npipi_K0_MAX = Npipi_K0_peak+2.*Npipi_K0_sigma;
static const double Npipi_Sp_MIN = Npipi_Sp_peak-2.*Npipi_Sp_sigma;
static const double Npipi_Sp_MAX = Npipi_Sp_peak+2.*Npipi_Sp_sigma;
static const double Npipi_Sm_MIN = Npipi_Sm_peak-2.*Npipi_Sm_sigma;
static const double Npipi_Sm_MAX = Npipi_Sm_peak+2.*Npipi_Sm_sigma;

static const double Npipi_K0_MIN25 = Npipi_K0_peak-2.5*Npipi_K0_sigma;
static const double Npipi_K0_MAX25 = Npipi_K0_peak+2.5*Npipi_K0_sigma;
static const double Npipi_Sp_MIN25 = Npipi_Sp_peak-2.5*Npipi_Sp_sigma;
static const double Npipi_Sp_MAX25 = Npipi_Sp_peak+2.5*Npipi_Sp_sigma;
static const double Npipi_Sm_MIN25 = Npipi_Sm_peak-2.5*Npipi_Sm_sigma;
static const double Npipi_Sm_MAX25 = Npipi_Sm_peak+2.5*Npipi_Sm_sigma;

static const double Npipi_K0_MIN3 = Npipi_K0_peak-3.*Npipi_K0_sigma;
static const double Npipi_K0_MAX3 = Npipi_K0_peak+3.*Npipi_K0_sigma;
static const double Npipi_Sp_MIN3 = Npipi_Sp_peak-3.*Npipi_Sp_sigma;
static const double Npipi_Sp_MAX3 = Npipi_Sp_peak+3.*Npipi_Sp_sigma;
static const double Npipi_Sm_MIN3 = Npipi_Sm_peak-3.*Npipi_Sm_sigma;
static const double Npipi_Sm_MAX3 = Npipi_Sm_peak+3.*Npipi_Sm_sigma;

static const double Npipi_K0_MIN4 = Npipi_K0_peak-4.*Npipi_K0_sigma;
static const double Npipi_K0_MAX4 = Npipi_K0_peak+4.*Npipi_K0_sigma;
static const double Npipi_Sp_MIN4 = Npipi_Sp_peak-4.*Npipi_Sp_sigma;
static const double Npipi_Sp_MAX4 = Npipi_Sp_peak+4.*Npipi_Sp_sigma;
static const double Npipi_Sm_MIN4 = Npipi_Sm_peak-4.*Npipi_Sm_sigma;
static const double Npipi_Sm_MAX4 = Npipi_Sm_peak+4.*Npipi_Sm_sigma;

static const double Npipi_K0_SB0_MIN = Npipi_K0_peak-7*Npipi_K0_sigma;
static const double Npipi_K0_SB0_MAX = Npipi_K0_peak-4*Npipi_K0_sigma;
static const double Npipi_K0_SB1_MIN = Npipi_K0_peak+4*Npipi_K0_sigma;
static const double Npipi_K0_SB1_MAX = Npipi_K0_peak+7*Npipi_K0_sigma;
static const double Npipi_K0_SB2_MIN = Npipi_K0_peak-8*Npipi_K0_sigma;
static const double Npipi_K0_SB2_MAX = Npipi_K0_peak-5*Npipi_K0_sigma;
static const double Npipi_K0_SB3_MIN = Npipi_K0_peak+5*Npipi_K0_sigma;
static const double Npipi_K0_SB3_MAX = Npipi_K0_peak+8*Npipi_K0_sigma;
static const double Npipi_K0_SB4_MIN = Npipi_K0_peak-9*Npipi_K0_sigma;
static const double Npipi_K0_SB4_MAX = Npipi_K0_peak-6*Npipi_K0_sigma;
static const double Npipi_K0_SB5_MIN = Npipi_K0_peak+6*Npipi_K0_sigma;
static const double Npipi_K0_SB5_MAX = Npipi_K0_peak+9*Npipi_K0_sigma;

//*** d(K-, n pi+ pi-) MM ***/
static const double KNpim_MM_Sp_peak  = 1.20033;
static const double KNpim_MM_Sp_sigma = 0.0176479;

static const double KNpip_MM_Sm_peak  = 1.18896;
static const double KNpip_MM_Sm_sigma = 0.0159155;

static const double KNpim_MM_Sp_MIN = KNpim_MM_Sp_peak - 2.*KNpim_MM_Sp_sigma;
static const double KNpim_MM_Sp_MAX = KNpim_MM_Sp_peak + 2.*KNpim_MM_Sp_sigma;

static const double KNpip_MM_Sm_MIN = KNpip_MM_Sm_peak - 2.*KNpip_MM_Sm_sigma;
static const double KNpip_MM_Sm_MAX = KNpip_MM_Sm_peak + 2.*KNpip_MM_Sm_sigma;

static const double Npipi_N_MIN = 0.9;
static const double Npipi_N_MAX = 0.98;

static const double Npipi_N_MIN1 = 0.9;
static const double Npipi_N_MAX1 = 0.96;

static const double Npipi_N_MIN2 = 0.9;
static const double Npipi_N_MAX2 = 0.94;

static const double Npipi_N_MIN0 = 0.87;
static const double Npipi_N_MAX0 = 1.05;

//*** for "n" pi -> S ***//
static const double MM_Npim_MIN = 1.165;
static const double MM_Npim_MAX = 1.225;

static const double MM_Npim_MIN2 = 1.155;
static const double MM_Npim_MAX2 = 1.235;

static const double MM_Npim_MIN3 = 1.145;
static const double MM_Npim_MAX3 = 1.245;

static const double MM_Npim_MIN4 = 1.135;
static const double MM_Npim_MAX4 = 1.255;

static const double MM_Npim_MIN5 = 1.125;
static const double MM_Npim_MAX5 = 1.265;


static const double MM_Npip_MIN = 1.155;
static const double MM_Npip_MAX = 1.215;

static const double MM_Npip_MIN2 = 1.145;
static const double MM_Npip_MAX2 = 1.225;

static const double MM_Npip_MIN3 = 1.135;
static const double MM_Npip_MAX3 = 1.235;

static const double MM_Npip_MIN4 = 1.125;
static const double MM_Npip_MAX4 = 1.245;

static const double MM_Npip_MIN5 = 1.115;
static const double MM_Npip_MAX5 = 1.255;

static const double MM_Sp_peak = 1.19625;
static const double MM_Sp_sigma = 0.0256449;

static const double MM_Sp_MIN2 = MM_Sp_peak-2.*MM_Sp_sigma;
static const double MM_Sp_MIN3 = MM_Sp_peak-3.*MM_Sp_sigma;
static const double MM_Sp_MIN4 = MM_Sp_peak-4.*MM_Sp_sigma;
static const double MM_Sp_MIN5 = MM_Sp_peak-5.*MM_Sp_sigma;
static const double MM_Sp_MAX2 = MM_Sp_peak+2.*MM_Sp_sigma;
static const double MM_Sp_MAX3 = MM_Sp_peak+3.*MM_Sp_sigma;
static const double MM_Sp_MAX4 = MM_Sp_peak+4.*MM_Sp_sigma;
static const double MM_Sp_MAX5 = MM_Sp_peak+5.*MM_Sp_sigma;

static const double MM_Sm_peak = 1.18739;
static const double MM_Sm_sigma = 0.024746;

static const double MM_Sm_MIN2 = MM_Sm_peak-2.*MM_Sm_sigma;
static const double MM_Sm_MIN3 = MM_Sm_peak-3.*MM_Sm_sigma;
static const double MM_Sm_MIN4 = MM_Sm_peak-4.*MM_Sm_sigma;
static const double MM_Sm_MIN5 = MM_Sm_peak-5.*MM_Sm_sigma;
static const double MM_Sm_MAX2 = MM_Sm_peak+2.*MM_Sm_sigma;
static const double MM_Sm_MAX3 = MM_Sm_peak+3.*MM_Sm_sigma;
static const double MM_Sm_MAX4 = MM_Sm_peak+4.*MM_Sm_sigma;
static const double MM_Sm_MAX5 = MM_Sm_peak+5.*MM_Sm_sigma;

static const double Lpim_P_MIN=0.89;
static const double Lpim_P_MAX=0.99;

//### for p(K-, n pi) ###
static const double KP_Npim_MM2_mean=0.019744;
static const double KP_Npim_MM2_sigma=0.0103519;
static const double KP_Npim_MM2_MIN=KP_Npim_MM2_mean-2.*KP_Npim_MM2_sigma;
static const double KP_Npim_MM2_MAX=KP_Npim_MM2_mean+2.*KP_Npim_MM2_sigma;

static const double KP_Npip_MM2_mean=0.0193675;
static const double KP_Npip_MM2_sigma=0.00981288;
static const double KP_Npip_MM2_MIN=KP_Npip_MM2_mean-2.*KP_Npip_MM2_sigma;
static const double KP_Npip_MM2_MAX=KP_Npip_MM2_mean+2.*KP_Npip_MM2_sigma;

static const double KP_K0_mean=0.499489;
static const double KP_K0_sigma=0.0152203;
static const double KP_K0_MAX=KP_K0_mean+2.*KP_K0_sigma;
static const double KP_K0_MIN=KP_K0_mean-2.*KP_K0_sigma;

//### for p(K-, L) ###
static const double KP_L_MM2_mean=0.0140127;
static const double KP_L_MM2_sigma=0.0173946;
static const double KP_L_MM2_MIN=KP_L_MM2_mean-2.*KP_L_MM2_sigma;
static const double KP_L_MM2_MAX=KP_L_MM2_mean+2.*KP_L_MM2_sigma;

#endif
