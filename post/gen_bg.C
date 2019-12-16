void gen_bg()
{

   if (!gROOT->GetClass("TGenPhaseSpace")) gSystem.Load("libPhysics");

   const double dMass = 1.87561;
   const double kpMass = 0.4936;
   const double piMass = 0.13957;
   const double nMass = 0.939565;
   
   TLorentzVector target(0.0, 0.0, 0.0, dMass);
   TLorentzVector beam(0.0, 0.0, 1.00,kpMass);
   TLorentzVector W = beam + target;
   
   double XMass=0.1;
   Double_t masses[3] = {piMass,piMass,XMass};

   TGenPhaseSpace event;
   event.SetDecay(W,3,masses);
   
   TH2F *h1 = new TH2F("h1","h1", 50,1.1,1.8, 50,1.1,1.8);

   for(Int_t i=0;i<10000;i++){
     Double_t weight = event.Generate();

     TLorentzVector *pim = event.GetDecay(0);
     TLorentzVector *pip = event.GetDecay(1);
     TLorentzVector *X   = event.GetDecay(2);

     TLorentzVector pippim = *pim+*pip;
     TLorentzVector pipX = *pip+*X;
     TLorentzVector pimX = *pim+*X;
     TLorentzVector Miss = W-*pip-*pim-*X;


   }








}

