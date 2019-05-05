void LightV(const char* filename)
{
  TFile *file = new TFile(filename,"READ");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(6,6);

  TH2F* CDH_diffpos_ctsub_pi_z_seg[36];
 
  gSystem->Load("./lib/libAll.so");
  
  ConfMan *conf = new ConfMan("/gpfs/home/had/hiasano/ana/k18ana/conf/Run78/analyzer_kwsk.conf");
  conf->SetRunNumber(800);
  conf->Initialize();

  //TF1* func[36];
  for(int iseg=0;iseg<36;iseg++){
    c1->cd(iseg+1);
    CDH_diffpos_ctsub_pi_z_seg[iseg] = (TH2F*)file->Get(Form("CDH_diffpos_ctsub_pi_z_seg%d",iseg+1));
    CDH_diffpos_ctsub_pi_z_seg[iseg]->Draw("colz");
    //func[iseg] = TF1("pol1",-4,4);
    CDH_diffpos_ctsub_pi_z_seg[iseg]->Fit("pol1","","",-2.5,2.5);
    double deltaV = pol1->GetParameter(1);
    std::cout << "seg " << iseg+1 << "  " << deltaV << std::endl;
    double lvcurrent=0;
    conf->GetGeomMapManager()->GetLightVelocity( CID_CDH, iseg+1, lvcurrent);
    double lvupdate = lvcurrent + deltaV;
    conf->GetGeomMapManager()->SetLightVelocity( CID_CDH, iseg+1, lvupdate);
  }
  
  ofstream ofs("geom.map");
  conf->GetGeomMapManager()->PrintMap(CID_CDH,ofs);
  ofs.close();
}
