MakeDtDx2(){
  gSystem->Load("~/k18ana/macro/HistTools_C.so");
  Main(164,173,"62","bldcraw");
}
Main(int run,int run2,char* aaa,TString prefix){
  std::cout<<"start run"<<run<<std::endl;
  double thre[]={0,0.3,2,5,10,17,25,32,
		 40,50,60,68,75,83,90,95,
		 98,99.7,100};
  int nthre=sizeof(thre)/sizeof(thre[0]);
  DetectorList *dlist=DetectorList::GetInstance();
  dlist->Initialize("param/Run62/Detectors.list");

  TFile *f[100];
  int nfile=run2-run+1;
  for(int irun=run;irun<=run2;irun++){
    f[irun-run]=new TFile(Form("~/k18ana/root/run%s/%s%d.root",aaa,prefix.Data(),irun));
  }
  if(!f[0]->IsOpen()) return;

  TH1F* h2;
  TTimeStamp *tm=new TTimeStamp();
  TString outname=Form("xt_%d_run%d_%d.param",tm->GetDate(),run,run2);
  if(run==run2)
    outname=Form("xt_%d_run%d.param",tm->GetDate(),run);
  ofstream ofs(outname.Data());
  Int_t cid[6]={CID_BPC,CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_FDC1};
  double dlmax[6]={3.6,4.,4.,2.5,2.5,3.0};
  //  double dlmaxfdc[6]={2.898,2.898,3.,3.,2.898,2.898};
  double tdctype[6]={-1,1,1,-1,1,-1};
  TH1F* temp;
  for(int ic=0;ic<6;ic++){
    char *name=dlist->GetName(cid[ic]).data();
    int nlayers=dlist->GetNlayers(cid[ic]);
    int nwires=dlist->GetNwires(cid[ic]);
    double tmpdlmax=dlmax[ic];
    std::cout<<name<<"  "<<nlayers<<"  "<<nwires<<std::endl;
    ofs<<"XTParam: "<<cid[ic]<<"  "<<tmpdlmax<<"  "<<nthre<<std::endl;
    ofs<<"Threshold: ";
    for(int i=0;i<nthre;i++)
      ofs<<thre[i]<<"  ";
    ofs<<std::endl;
    for(int lay=1;lay<=nlayers;lay++){
      for(int wire=1;wire<=nwires;wire++){
	TString hname=Form("h_dt_%s_%d_%d",name,lay,wire);
	h1=(TH1F*)SumHist(f,hname,nfile);
	//	std::cout<<hname<<std::endl;
	if(!h1){
	  std::cout<<" !!! "<<hname<<std::endl;
	  continue;
	}
	double sum=0;
	double sumbefore=0;
	double nev=h1->GetEntries();
	double tmp;
	if(nev<1) continue;

	int tmpn=0;
	int tmpn2=0;
	double x[30],y[30];
	double x2[1000],y2[1000];
	x[tmpn]=h1->GetBinCenter(0);
	y[tmpn]=0;
	tmpn++;
	double tmpth=thre[tmpn]/100.;
	for(int ib=0;ib<h1->GetNbinsX();ib++){
	  sumbefore=sum;
	  sum += h1->GetBinContent(ib);
	  double tmpdx=(double)sum/(double)nev;
	  if(tmpn<nthre&&tmpdx>tmpth){
	    double tmpx1=h1->GetBinCenter(ib-1);
	    double tmpx2=h1->GetBinCenter(ib);
	    double tmpy1=sumbefore;
	    double tmpy2=sum;
	    if(tmpy1==tmpy2) x[tmpn]=tmpx1;
	    else
	      x[tmpn]=tmpx1+h1->GetBinWidth(ib)*(tmpth*nev-tmpy1)/(tmpy2-tmpy1);
	    y[tmpn]=tmpth*tmpdlmax;
	    //	    std::cout<<x[tmpn]<<"  "<<y[tmpn]<<std::endl;
	    tmpn++;
	    if(tmpn<nthre)   tmpth=thre[tmpn]/100.;
	  }
	}
	x[tmpn]=h1->GetBinCenter(h1->GetNbinsX()+2);
	y[tmpn]=tmpdlmax;
	ofs<<std::setw(5)<<cid[ic]
	   <<std::setw(5)<<lay
	   <<std::setw(5)<<wire;
	ofs.setf(std::ios_base::fixed,std::ios_base::floatfield);
	for(int j=0;j<nthre;j++)
	  ofs<<std::setw(8)<<std::setprecision(2)<<x[j];
	ofs<<std::endl;
      }
    }
  }
  ofs.close();
  for(int irun=run;irun<=run2;irun++){
    f[irun-run]->Close();
  }
}
