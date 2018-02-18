void tof( TString hname = "tmp.root" )
{
  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetLogy();

  TOF_T0_TOF->SetAxisRange(-5,40);
  TOF_T0_TOF->SetStats(0);
  TOF_T0_TOF->Draw();
  double maxy = TOF_T0_TOF->GetMaximum();
  TOF_T0_TOF_E->SetLineColor(7);
  TOF_T0_TOF_E->Draw("same");
  TOF_T0_TOF_Pi->SetLineColor(4);
  TOF_T0_TOF_Pi->Draw("same");
  TOF_T0_TOF_K->SetLineColor(3);
  TOF_T0_TOF_K->Draw("same");
  TOF_T0_TOF_P->SetLineColor(2);
  TOF_T0_TOF_P->Draw("same");

  TLatex latex;
  latex.SetTextColor(1); latex.DrawLatex( 10, 0.8*maxy, "no selection" );
  latex.SetTextColor(7); latex.DrawLatex( 10, 0.4*maxy, "e trigger" );
  latex.SetTextColor(4); latex.DrawLatex( 10, 0.2*maxy, "pi trigger" );
  latex.SetTextColor(3); latex.DrawLatex( 10, 0.1*maxy, "K trigger" );
  latex.SetTextColor(2); latex.DrawLatex( 10, 0.05*maxy, "p trigger" );

}

void toff2_t0( TString hname = "tmp.root" )
{
  double tgatemin = -1;
  double tgatemax = 1;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );

  TH1F *h1,*h2,*h3;
  TF1 *fun1, *fun2;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(2,3);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->Divide(2,3);

  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    double xTOF,yTOF,zTOF, gxTOF,gyTOF,gzTOF;
    double xT0,yT0,zT0, gxT0,gyT0,gzT0;
    conf->GetGeomMapManager()->GetParam( CID_TOFstop, 19, xTOF, yTOF, zTOF, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );
    conf->GetGeomMapManager()->GetGParam( CID_TOFstop, gxTOF, gyTOF, gzTOF, tmp1, tmp2, tmp3 );
    conf->GetGeomMapManager()->GetParam( CID_T0, seg, xT0, yT0, zT0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );
    conf->GetGeomMapManager()->GetGParam( CID_T0, gxT0, gyT0, gzT0, tmp1, tmp2, tmp3 );
    xTOF += gxTOF; yTOF += gyTOF; zTOF += gzTOF;
    xT0 += gxT0; yT0 += gyT0; zT0 += gzT0;
    double Length = sqrt( (xTOF-xT0)*(xTOF-xT0) + (yTOF-yT0)*(yTOF-yT0) + (zTOF-zT0)*(zTOF-zT0) );

    h1 = (TH1F*)f->Get( Form("TOF_T0%d_TOF19",seg) );
    c1->cd(seg);
    //h1->Draw();
    if( h1->Integral( h1->GetXaxis()->FindBin(tgatemin), h1->GetXaxis()->FindBin(tgatemax) )<5 ){
      //std::cout << " seg" << seg << " skip " << std::endl;
      continue;
    }
    h1->Fit("gaus","q","",tgatemin,tgatemax);
    fun1 = (TF1*)h1->GetFunction("gaus");
    
    h2 = (TH1F*)f->Get( Form("ctsubT0%d",seg) );
    c2->cd(seg);
    //h2->Draw();
    if( h2->Integral( h2->GetXaxis()->FindBin(tgatemin), h2->GetXaxis()->FindBin(tgatemax) )<5 ){
      //std::cout << " seg" << seg << " skip " << std::endl;
      continue;
    }
    h2->Fit("gaus","q","",tgatemin,tgatemax);
    fun2 = (TF1*)h2->GetFunction("gaus");

    double p1 = fun1->GetParameter(1);
    double p2 = fun2->GetParameter(1);

    double of1, of2;
    of1 = p1+p2;
    of2 = p1-p2;

    int c,n,a;
    double gain, offset;
    int at = 1, ud = 0;
    conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset -= of1;
    std::cout << "u seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " offset:" << offset << " gain:" << gain << " len:" << Length/30. << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
    ud = 1;
    conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset -= of2;
    std::cout << "d seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " offset:" << offset << " gain:" << gain << " len:" << Length/30. << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
  }

  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void toff2_tof( TString hname = "tmp.root" )
{
  double tgatemin1 = -5;
  double tgatemax1 =  5;
  double tgatemin2 = -5;
  double tgatemax2 =  5;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );

  TH1F *h1,*h2,*h3;
  TF1 *fun1, *fun2;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);
  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  c2->Divide(8,4);

  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){

    if( seg==8 || seg>=28 ) continue;
    
    h1 = (TH1F*)f->Get( Form("TOF_T03_TOF%d",seg) );
    c1->cd(seg);
    h1->SetAxisRange(tgatemin1,tgatemax1,"X");
    //h1->Draw();
    if( h1->Integral( h1->GetXaxis()->FindBin(tgatemin1), h1->GetXaxis()->FindBin(tgatemax1) )<10 ){
      //std::cout << " seg" << seg << " skip " << std::endl;
      continue;
    }
    h1->Fit("gaus","q","",tgatemin1,tgatemax1);
    fun1 = (TF1*)h1->GetFunction("gaus");
    
    h2 = (TH1F*)f->Get( Form("ctsubTOF%d",seg) );
    c2->cd(seg);
    h2->SetAxisRange(50,4000,"X");
    //h2->Draw();
    if( h2->Integral( h2->GetXaxis()->FindBin(tgatemin2), h2->GetXaxis()->FindBin(tgatemax2) )<10 ){
      //std::cout << " seg" << seg << " skip " << std::endl;
      continue;
    }
    h2->Fit("gaus","q","",tgatemin2,tgatemax2);
    fun2 = (TF1*)h2->GetFunction("gaus");

    double p1 = fun1->GetParameter(1);
    double p2 = fun2->GetParameter(1);

    double of1, of2;
    of1 = p1+p2;
    of2 = p1-p2;

    int c,n,a;
    double gain, offset;
    int at = 1, ud = 0;
    conf->GetCounterMapManager()->GetCNA( CID_TOFstop, seg, at, ud, c, n, a );
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset += of1;
    std::cout << "u seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " of1:" << of1 << " offset:" << offset << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
    ud = 1;
    conf->GetCounterMapManager()->GetCNA( CID_TOFstop, seg, at, ud, c, n, a );
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset += of2;
    std::cout << "d seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " of2:" << of2 << " offset:" << offset << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
  }
  
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();

}

void check_e( TString hname = "tmp.root" )
{
  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("eTOFU%d",seg) );
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("eTOFD%d",seg) );
    h1->SetLineColor(2); h1->Draw("same");
  }
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->Divide(2,3);
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){ 
    c2->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("eT0U%d",seg) );
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("eT0D%d",seg) );
    h1->SetLineColor(2); h1->Draw("same");
  }
}

void check_t( TString hname = "tmp.root" )
{
  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    c1->cd(seg);
    h1 = (TH1F*)f->Get( Form("tTOFU%d",seg) );
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("tTOFD%d",seg) );
    h1->SetLineColor(2); h1->Draw("same");
  }
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->Divide(2,3);
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){ 
    c2->cd(seg);
    h1 = (TH1F*)f->Get( Form("tT0U%d",seg) );
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("tT0D%d",seg) );
    h1->SetLineColor(2); h1->Draw("same");
  }
}

void slew_t0u( TString hname = "tmp.root" )
{
  double egatemin = 0.5;
  double egatemax = 2.;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(4,4);

  TF1 *slewfunc = new TF1( "slewfunc", "[0]/sqrt(x)+[1]", egatemin, egatemax );
  slewfunc->SetNpx(999);
  //slewfunc->SetParameters(0,0);
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    h1 = (TH1F*)f->Get( Form("ectT0U%d",seg) );
    if( h1->Integral( h1->GetXaxis()->FindBin(egatemin), h1->GetXaxis()->FindBin(egatemax) )<10 ) continue;
    c1->cd(seg);
    std::cout << "AAA" << slewfunc << std::endl;
    h1->Fit("slewfunc");
    int ud=0, ith=0, type=1, npar=2;
    double p0 = -h1->GetFunction("slewfunc")->GetParameter(0);
    double p1 = -h1->GetFunction("slewfunc")->GetParameter(1);
    //double par[2] = { p0,p1 };
    double par[2] = { p0,0 };
    conf->GetSlewingMapManager()->SetParam( CID_T0, seg, ud, ith, type, npar, par );
  }
  ofstream ofs("tmp.dat");
  conf->GetSlewingMapManager()->PrintMap(CID_T0,ofs);
  ofs.close();
}

void slew_t0d( TString hname = "tmp.root" )
{
  double egatemin = 0.5;
  double egatemax = 2.;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(4,4);

  TF1 *slewfunc = new TF1( "slewfunc", "[0]/sqrt(x)+[1]", egatemin, egatemax );
  slewfunc->SetNpx(999);
  //slewfunc->SetParameters(0,0);
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    h1 = (TH1F*)f->Get( Form("ectT0D%d",seg) );
    if( h1->Integral( h1->GetXaxis()->FindBin(egatemin), h1->GetXaxis()->FindBin(egatemax) )<10 ) continue;
    c1->cd(seg);
    std::cout << "AAA" << slewfunc << std::endl;
    h1->Fit("slewfunc");
    int ud=1, ith=0, type=1, npar=2;
    double p0 = -h1->GetFunction("slewfunc")->GetParameter(0);
    double p1 = -h1->GetFunction("slewfunc")->GetParameter(1);
    //double par[2] = { p0,p1 };
    double par[2] = { p0,0 };
    conf->GetSlewingMapManager()->SetParam( CID_T0, seg, ud, ith, type, npar, par );
  }
  ofstream ofs("tmp.dat");
  conf->GetSlewingMapManager()->PrintMap(CID_T0,ofs);
  ofs.close();
}

void slew_tofu( TString hname = "tmp.root" )
{
  double egatemin = 0.5;
  double egatemax = 49.5;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);

  TF1 *slewfunc = new TF1( "slewfunc", "[0]/sqrt(x)+[1]", egatemin, egatemax );
  slewfunc->SetNpx(999);
  //slewfunc->SetParameters(0,0);
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    //h1 = (TH1F*)f->Get( Form("etTOFU%d",seg) );
    h1 = (TH1F*)f->Get( Form("dETOF_T03_TOFU%d",seg) );
    if( h1->Integral( h1->GetXaxis()->FindBin(egatemin), h1->GetXaxis()->FindBin(egatemax) )<10 ) continue;
    c1->cd(seg);
    std::cout << "AAA" << slewfunc << std::endl;
    h1->Fit("slewfunc");
    int ud=0, ith=0, type=1, npar=2;
    double p0 = -h1->GetFunction("slewfunc")->GetParameter(0);
    double p1 = -h1->GetFunction("slewfunc")->GetParameter(1);
    //double par[2] = { p0,p1 };
    double par[2] = { p0,0 };
    conf->GetSlewingMapManager()->SetParam( CID_TOFstop, seg, ud, ith, type, npar, par );
  }
  ofstream ofs("tmp.dat");
  conf->GetSlewingMapManager()->PrintMap(CID_TOFstop,ofs);
  ofs.close();
}

void slew_tofd( TString hname = "tmp.root" )
{
  double egatemin = 0.5;
  double egatemax = 49.5;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);

  TF1 *slewfunc = new TF1( "slewfunc", "[0]/sqrt(x)+[1]", egatemin, egatemax );
  slewfunc->SetNpx(999);
  //slewfunc->SetParameters(0,0);
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    //h1 = (TH1F*)f->Get( Form("etTOFU%d",seg) );
    h1 = (TH1F*)f->Get( Form("dETOF_T03_TOFD%d",seg) );
    if( h1->Integral( h1->GetXaxis()->FindBin(egatemin), h1->GetXaxis()->FindBin(egatemax) )<10 ) continue;
    //if( seg==9 ) continue;
    c1->cd(seg);
    std::cout << "AAA" << slewfunc << std::endl;
    h1->Fit("slewfunc");
    int ud=1, ith=0, type=1, npar=2;
    double p0 = -h1->GetFunction("slewfunc")->GetParameter(0);
    double p1 = -h1->GetFunction("slewfunc")->GetParameter(1);
    //double par[2] = { p0,p1 };
    double par[2] = { p0,0 };
    conf->GetSlewingMapManager()->SetParam( CID_TOFstop, seg, ud, ith, type, npar, par );
  }
  ofstream ofs("tmp.dat");
  conf->GetSlewingMapManager()->PrintMap(CID_TOFstop,ofs);
  ofs.close();
}

void toff_tofu( TString hname = "tmp.root" )
{
  double tgatemin = 1000;
  double tgatemax = 2000;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);

  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    h1 = (TH1F*)f->Get( Form("TTOFU%d",seg) );
    c1->cd(seg);
    h1->SetAxisRange(50,4000,"X");
    //h1->Draw();
    if( h1->Integral( h1->GetXaxis()->FindBin(tgatemin), h1->GetXaxis()->FindBin(tgatemax) )==0 ) continue;
    h1->Fit("gaus","q","",tgatemin, tgatemax );
    int c,n,a;
    int at = 1, ud = 0;
    conf->GetCounterMapManager()->GetCNA( CID_TOFstop, seg, at, ud, c, n, a );
    double gain, offset;
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset = h1->GetFunction("gaus")->GetParameter(1)*gain;
    std::cout << " seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " offset:" << offset << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
  }
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void toff_tofd( TString hname = "tmp.root" )
{
  double tgatemin = 1000;
  double tgatemax = 2000;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);

  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    h1 = (TH1F*)f->Get( Form("TTOFD%d",seg) );
    c1->cd(seg);
    h1->SetAxisRange(50,4000,"X");
    //h1->Draw();
    if( h1->Integral( h1->GetXaxis()->FindBin(tgatemin), h1->GetXaxis()->FindBin(tgatemax) )==0 ) continue;
    h1->Fit("gaus","q","",tgatemin, tgatemax );
    int c,n,a;
    int at = 1, ud = 1;
    conf->GetCounterMapManager()->GetCNA( CID_TOFstop, seg, at, ud, c, n, a );
    double gain, offset;
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset = h1->GetFunction("gaus")->GetParameter(1)*gain;
    std::cout << " seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " offset:" << offset << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
  }
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void toff_t0u( TString hname = "tmp.root" )
{
  double tgatemin = 500;
  double tgatemax = 1000;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(2,3);

  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    h1 = (TH1F*)f->Get( Form("TT0U%d",seg) );
    c1->cd(seg);
    h1->SetAxisRange(50,4000,"X");
    //h1->Draw();
    if( h1->Integral( h1->GetXaxis()->FindBin(tgatemin), h1->GetXaxis()->FindBin(tgatemax) )==0 ) continue;
    h1->Fit("gaus","q","",tgatemin, tgatemax );
    int c,n,a;
    int at = 1, ud = 0;
    conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
    double gain, offset;
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset = h1->GetFunction("gaus")->GetParameter(1)*gain;
    std::cout << " offset:" << offset << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
  }
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void toff_t0d( TString hname = "tmp.root" )
{
  double tgatemin = 600;
  double tgatemax = 700;

  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(2,3);

  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    h1 = (TH1F*)f->Get( Form("TT0D%d",seg) );
    c1->cd(seg);
    h1->SetAxisRange(50,4000,"X");
    //h1->Draw();
    if( h1->Integral( h1->GetXaxis()->FindBin(tgatemin), h1->GetXaxis()->FindBin(tgatemax) )==0 ) continue;
    h1->Fit("gaus","q","",tgatemin, tgatemax );
    int c,n,a;
    int at = 1, ud = 1;
    conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
    double gain, offset;
    conf->GetGainMapManager()->GetParam(c,n,a,gain,offset);
    offset = h1->GetFunction("gaus")->GetParameter(1)*gain;
    std::cout << " offset:" << offset << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, offset );
  }
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void ped_tofu( TString hname = "tmp.root" )
{
  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    h1 = (TH1F*)f->Get( Form("ATOFU%d",seg) );
    h2 = (TH1F*)f->Get( Form("AwTTOFU%d",seg) );
    h3 = new TH1F(*h1);
    h3->Add(h1,h2,1,-1);
    c1->cd(seg);
    h2->Fit("landau","q","",0,1500);
    h3->Fit("gaus","q","",0,200);
    gPad->SetLogy();
    h1->SetAxisRange(0,1500,"X");
    h1->Draw();
    h2->SetLineColor(2); h2->Draw("same");
    h3->SetLineColor(3); h3->Draw("same");
    //std::cout << "MPV:" << h2->GetFunction("landau")->GetParameter(1) << " Mean:" << h2->GetMean() << std::endl;
    double ped = h3->GetFunction("gaus")->GetParameter(1);
    double gain = 6./(h2->GetFunction("landau")->GetParameter(1)-h3->GetFunction("gaus")->GetParameter(1));
    int c,n,a;
    int at = 0, ud = 0;
    conf->GetCounterMapManager()->GetCNA( CID_TOFstop, seg, at, ud, c, n, a );
    std::cout << " seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " ped:" << ped << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, ped );
  }
  
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void ped_tofd( TString hname = "tmp.root" )
{
  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(8,4);
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    h1 = (TH1F*)f->Get( Form("ATOFD%d",seg) );
    h2 = (TH1F*)f->Get( Form("AwTTOFD%d",seg) );
    h3 = new TH1F(*h1);
    h3->Add(h1,h2,1,-1);
    c1->cd(seg);
    h2->Fit("landau","q","",0,1500);
    h3->Fit("gaus","q","",0,200);
    gPad->SetLogy();
    h1->SetAxisRange(0,1500,"X");
    h1->Draw();
    h2->SetLineColor(2); h2->Draw("same");
    h3->SetLineColor(3); h3->Draw("same");
    //std::cout << "MPV:" << h2->GetFunction("landau")->GetParameter(1) << " Mean:" << h2->GetMean() << std::endl;
    double ped = h3->GetFunction("gaus")->GetParameter(1);
    double gain = 6./(h2->GetFunction("landau")->GetParameter(1)-h3->GetFunction("gaus")->GetParameter(1));
    int c,n,a;
    int at = 0, ud = 1;
    conf->GetCounterMapManager()->GetCNA( CID_TOFstop, seg, at, ud, c, n, a );
    std::cout << " seg:" << seg << " c,n,a:" << c << "," << n << "," << a << " ped:" << ped << " gain:" << gain << std::endl;
    conf->GetGainMapManager()->SetParam( c, n, a, gain, ped );
  }
  
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void ped_t0u( TString hname = "tmp.root" )
{
  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(2,2);
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    h1 = (TH1F*)f->Get( Form("AT0U%d",seg) );
    h2 = (TH1F*)f->Get( Form("AwTT0U%d",seg) );
    h3 = new TH1F(*h1);
    h3->Add(h1,h2,1,-1);
    c1->cd(seg);
    h2->Fit("landau","q");
    h3->Fit("gaus","q");
    gPad->SetLogy();
    h1->SetAxisRange(0,500,"X");
    h1->Draw();
    h2->SetLineColor(2); h2->Draw("same");
    h3->SetLineColor(3); h3->Draw("same");
    //std::cout << "MPV:" << h2->GetFunction("landau")->GetParameter(1) << " Mean:" << h2->GetMean() << std::endl;
    double ped = h3->GetFunction("gaus")->GetParameter(1);
    double gain = 1./(h2->GetFunction("landau")->GetParameter(1)-h3->GetFunction("gaus")->GetParameter(1));
    std::cout << " ped:" << ped << " gain:" << gain << std::endl;
    
    int c,n,a;
    int at = 0, ud = 0;
    conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
    conf->GetGainMapManager()->SetParam( c, n, a, gain, ped );
  }
  
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

void ped_t0d( TString hname = "tmp.root" )
{
  gSystem->Load("./lib/libAll.so");
  TFile *f = new TFile( hname );
  TH1F *h1,*h2,*h3;
  
  ConfMan *conf = new ConfMan("./conf/Oct2010/analyzer.conf");
  conf->Initialize();

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->Divide(2,2);
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    h1 = (TH1F*)f->Get( Form("AT0D%d",seg) );
    h2 = (TH1F*)f->Get( Form("AwTT0D%d",seg) );
    h3 = new TH1F(*h1);
    h3->Add(h1,h2,1,-1);
    c1->cd(seg);
    h2->Fit("landau","q");
    h3->Fit("gaus","q");
    gPad->SetLogy();
    h1->SetAxisRange(0,500,"X");
    h1->Draw();
    h2->SetLineColor(2); h2->Draw("same");
    h3->SetLineColor(3); h3->Draw("same");
    //std::cout << "MPV:" << h2->GetFunction("landau")->GetParameter(1) << " Mean:" << h2->GetMean() << std::endl;
    double ped = h3->GetFunction("gaus")->GetParameter(1);
    double gain = 1./(h2->GetFunction("landau")->GetParameter(1)-h3->GetFunction("gaus")->GetParameter(1));
    std::cout << " ped:" << ped << " gain:" << gain << std::endl;
    
    int c,n,a;
    int at = 0, ud = 1;
    conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
    conf->GetGainMapManager()->SetParam( c, n, a, gain, ped );
  }
  
  ofstream ofs("tmp.dat");
  conf->GetGainMapManager()->PrintMapHeader(ofs);
  conf->GetGainMapManager()->PrintMap(2,ofs);
  conf->GetGainMapManager()->PrintMap(1,ofs);
  conf->GetGainMapManager()->PrintMap(0,ofs);
  ofs.close();
}

