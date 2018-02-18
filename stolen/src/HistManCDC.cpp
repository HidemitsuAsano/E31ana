#include "HistManBeamAna.h"

void HistManBeamAna::initCDC()
{
  for(int lay=1; lay<=NumOfCDCLayers; lay++ ){
    new TH1F(Form("CDC_resi_%d", lay),    Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH1F(Form("CDC_resi_pi_%d", lay), Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH1F(Form("CDC_resi_p_%d", lay),  Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH1F(Form("CDC_resi_k_%d", lay),  Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH2F(Form("CDC_dt_resi_%d", lay), Form("CDC dt vs residual layer%d", lay), 300, -20, 280, 200, -0.1, 0.1);
    for( int wire=1; wire<=NumOfCDCWiresInLayer[lay-1]; wire++ ){
      new TH1F(Form("CDC_resi_%d_%d", lay, wire),    Form("CDC residual l%d w%d", lay, wire), 200, -0.1, 0.1);
      new TH1F(Form("CDC_resi_pi_%d_%d", lay, wire), Form("CDC residual l%d %d", lay, wire), 200, -0.1, 0.1);
      new TH1F(Form("CDC_resi_p_%d_%d", lay, wire),  Form("CDC residual l%d w%d", lay, wire), 200, -0.1, 0.1);
      new TH1F(Form("CDC_resi_k_%d_%d", lay, wire),  Form("CDC residual l%d w%d", lay, wire), 200, -0.1, 0.1);
      new TH2F(Form("CDC_dt_resi_%d_%d", lay, wire), Form("CDC dt vs residual l%d w%d", lay, wire), 300, -20, 280,  200, -0.1, 0.1);
      new TH2F(Form("CDC_dt_dl_%d_%d", lay, wire)  , Form("CDC dt vs dl l%d w%d", lay, wire),       300, -20, 280, 1000, -0.1, 0.9);
    }
  }
}
