const int CID = 2;

void bhdall( TString fname = "root/out730.root" )
{
  gSystem->Load("lib/libAll.so");
  ConfMan *conf = new ConfMan("conf/730.conf");
  conf->Initialize();
  int cid, seg, at, ud;
  int c,n,a;
  double gain, offset;

  TFile *f = new TFile( fname );
  cid=2; seg= 7; at=1; ud=0;
  bhd( f, seg, ud, -5,  5 ); conf->GetCounterMapManager()->GetCNA(cid,seg,at,ud,c,n,a); conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
  cid=2; seg= 8; at=1; ud=0;
  bhd( f, seg, ud, -5,  5 ); conf->GetCounterMapManager()->GetCNA(cid,seg,at,ud,c,n,a); conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
//   bhd( f, seg, ud, -5, -2 );
//   bhd( f, seg, ud, -5, -2 );

}

void bhd( TFile *f, int seg=1, int ud=0, double xmin=-5, double xmax=5 )
{
  TString name;
  if( ud==0 ) name = Form( "BHD%dtimeu", seg );
  else       name = Form( "BHD%dtimed", seg );
  TH1F *h1=0;
  h1 = (TH1F*)f->Get(name);
  if(h1==0) return; 
  h1->Fit("gaus","","",xmin,xmax);
}
