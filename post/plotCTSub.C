#include <iostream>
#include <TSystem.h>
#include <TH1.h>
#include <TFile.h>

class ConfMan;

void plotCTSub()
{
  gSystem->Load("./lib/libAll.so");
  
  ConfMan *conf = new ConfMan("/gpfs/home/had/hiasano/ana/k18ana/conf/Run78/analyzer_kwsk.conf");
  std::ofstream fout("ctsublog");
  for(int irun=100;irun<812;irun++){
    //{int irun=601;
    char filename[1024];
    sprintf(filename,"/gpfs/group/had/knucl/e15/asano/Run78/IMpisigmav140/evanaIMpisigma_0%03d.root",irun);
    std::cout << filename << std::endl;
    TFile *fin = new TFile(filename,"READ");
    if(!fin->IsOpen()){ 
      delete fin;
      continue;
    }
    std::cout << "run # " << irun << std::endl;
    fout << "run # " << irun << std::endl;
    conf->SetRunNumber(irun);
    conf->Initialize();
    std::cout << "couter Map " << conf->GetCounterMapManager()->GetFileName().c_str() << std::endl;
    std::string filenamecds = conf->GetGainMapManager()->GetFileNameCDS();
    std::cout << "cdsfile name " << filenamecds.c_str() << std::endl;
    filenamecds.replace(0,28,"");
    filenamecds.replace(filenamecds.size()-8,2,"36");
    std::cout << "outfile name " << filenamecds.c_str() << std::endl;
    //TCanvas *cctsub = new TCanvas("cctsub","cctsub",1800,1000);
    //cctsub->Divide(6,6);
    for(int icdh=0;icdh<36;icdh++){
      //cctsub->cd(icdh+1);
      TH1F *hctsub = (TH1F*)fin->Get(Form("CTSub%d",icdh+1));
      
      //hctsub->Draw();
      double ctsuboff = 0;
      if(hctsub !=NULL){
        ctsuboff = hctsub->GetMean();
        //std::cout << "icdh " << icdh+1 << " " << hctsub->GetMean() << std::endl;
        fout << "icdh " << icdh+1 << " " << hctsub->GetMean() << std::endl;
      }else{
        //std::cout << "icdh " << icdh+1 << " " << 0 << std::endl;
        fout << "icdh " << icdh+1 << " " << 0 << std::endl;
      }
      delete hctsub;
      int cr[2],sl[2],ch[2];//up,down
      const int UP=0;
      const int DOWN=1;
      conf->GetCounterMapManager()->GetCNA(CID_CDH,icdh+1,1,UP,cr[UP],sl[UP],ch[UP]);
      conf->GetCounterMapManager()->GetCNA(CID_CDH,icdh+1,1,DOWN,cr[DOWN],sl[DOWN],ch[DOWN]);
      double temp[2]={0.,0.};
      double offs[2]={0.,0.};
      //up 
      conf->GetGainMapManager()->GetParam( cr[UP],sl[UP],ch[UP], temp[UP], offs[UP] );
      conf->GetGainMapManager()->SetParam( cr[UP],sl[UP],ch[UP], temp[UP], offs[UP]+ctsuboff);
      //down
      conf->GetGainMapManager()->GetParam( cr[DOWN],sl[DOWN],ch[DOWN], temp[DOWN], offs[DOWN] );
      conf->GetGainMapManager()->SetParam( cr[DOWN],sl[DOWN],ch[DOWN], temp[DOWN], offs[DOWN]-ctsuboff);
      
    }
    //delete conf;
    fin->Close();
    std::ofstream ofs(filenamecds.c_str());
    conf->GetGainMapManager()->PrintMapCDS(ofs);
    delete fin;
    ofs.close();
    //conf->Clear();
  }
}

