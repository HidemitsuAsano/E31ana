//asano memo
//macro to create X-T curve ? at first.
//
//Input : output of evanaraw
//Output : x-t curve 

MakeDtDx2(int run_start,int run_end){
  gSystem->Load("../lib/libAll.so");
  std::cout<<"start run"<<run_start<<std::endl;
  //asano memo
  //
  double thre[]={0,0.3,2,5,10,17,25,32,
		 40,50,60,68,75,83,90,95,
		 98,99.7,100};
  int nthre=sizeof(thre)/sizeof(thre[0]);
  DetectorList *dlist=DetectorList::GetInstance();
  dlist->Initialize("../param/Run74/Detectors.list");

  TFile *f[100];
  int nfile=run_end-run_start+1;
  for(int irun=run_start;irun<=run_end;irun++){
    f[irun-run_start]=new TFile(Form("../rawroot/run74raw_00%02d.root",irun));
  }
  if(!f[0]->IsOpen()) return;

  TH1F* h2;
  TTimeStamp *tm=new TTimeStamp();
  TString outname=Form("xt_%d_run%d_%d.param",tm->GetDate(),run_start,run_end);
  if(run_start==run_end)
    outname=Form("xt_%d_run%d.param",tm->GetDate(),run_start);
  ofstream ofs(outname.Data());

  //CID is defined in GlobalVariables.h
  Int_t cid[6]={CID_BPC,CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_FDC1};
  //max drift length (mm)
  double dlmax[6]={3.0,4.,4.,2.5,2.5,3.0};
  //  double dlmaxfdc[6]={2.898,2.898,3.,3.,2.898,2.898};
  // ?
  double tdctype[6]={-1,1,1,-1,1,-1};
  TH1F* htemp=NULL;
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
	TString hname=Form("T%s_%d_%d",name,lay,wire);
	//htemp=(TH1F*)SumHist(f,hname,nfile);
	htemp=(TH1F*)f->Get(hname);
	//	std::cout<<hname<<std::endl;
	if(!htemp){
	  std::cout<<" !!! "<<hname<<std::endl;
	  continue;
	}
	double sum=0;
	double sumbefore=0;
	double nev=htemp->GetEntries();
	double tmp;
	if(nev<1) continue;

	int tmpn=0;
	int tmpn2=0;
	double x[30],y[30];
	double x2[1000],y2[1000];
	x[tmpn]=htemp->GetBinCenter(0);
	y[tmpn]=0;
	tmpn++;
	double tmpth=thre[tmpn]/100.;
	for(int ib=0;ib<htemp->GetNbinsX();ib++){
	  sumbefore=sum;
	  sum += htemp->GetBinContent(ib);
	  double tmpdx=(double)sum/(double)nev;
	  if(tmpn<nthre&&tmpdx>tmpth){
	    double tmpx1=htemp->GetBinCenter(ib-1);
	    double tmpx2=htemp->GetBinCenter(ib);
	    double tmpy1=sumbefore;
	    double tmpy2=sum;
	    if(tmpy1==tmpy2) x[tmpn]=tmpx1;
	    else
	      x[tmpn]=tmpx1+htemp->GetBinWidth(ib)*(tmpth*nev-tmpy1)/(tmpy2-tmpy1);
	    y[tmpn]=tmpth*tmpdlmax;
	    //	    std::cout<<x[tmpn]<<"  "<<y[tmpn]<<std::endl;
	    tmpn++;
	    if(tmpn<nthre)   tmpth=thre[tmpn]/100.;
	  }
	}
	x[tmpn]=htemp->GetBinCenter(htemp->GetNbinsX()+2);
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
  for(int irun=run_start;irun<=run_end;irun++){
    f[irun-run_start]->Close();
  }
}
