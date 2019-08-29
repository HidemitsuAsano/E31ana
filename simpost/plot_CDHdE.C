void plot_CDHdE()
{

    TFile *_file0 = TFile::Open("simIMpisigma_nSmpip_v55.root");
    CDHdE->ProjectionY("pyiron");
    pyiron->SetLineColor(4);
    pyiron->SetXTitle("CDH dE [MeVee]");
    pyiron->Draw("HE");

    TFile *_file1 = TFile::Open("simIMpisigma_nSmpip_v54.root");
    CDHdE->ProjectionY("pynew");
    pynew->SetLineColor(2);
    pynew->Scale(0.95);
    pynew->Draw("HEsame");


    TFile *_file2 = TFile::Open("simIMpisigma_nSmpip_DoraAir_v70.root");
    CDHdE->ProjectionY("pyair");
    pyair->SetLineColor(1);
    pyair->Draw("HEsame");

    TFile *_file3 = TFile::Open("../post/evanaIMpisigma_v159.root");
    //CDHdE_woK0_wSid_n->ProjectionY("pydata");
    CDHdE_pippim->ProjectionY("pydata");
    pydata->SetLineColor(5);
    double nair = pyair->Integral();
    double ndata = pydata->Integral();
    pydata->Scale(nair/ndata);
    pydata->Scale(1.3);
    //pydata->SetMarkerStyle(22);
    //pydata->SetMarkerColor(5);
    pydata->GetXaxis()->SetLimits(0,8.8);
    pydata->Scale(1.13);
    pydata->Draw("HEsame");








}
