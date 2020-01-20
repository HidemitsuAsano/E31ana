//Missing mass
Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(x[0]<1.116){
    return par[0]*exp(-0.5*pow((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2.0)); 
  }else if(1.116<=x[0] && x[0]<1.5){
    return par[4]*exp(-0.5*pow(((x[0]-par[5])/par[6]),2.0)); 
  }else{
    return 1.;
  }
}


Double_t param_MMnmiss[7]={
  1.64736,
  0.906359,
  0.235915,
  -0.929384,
  2.58101,
  1.48951,
  0.265225};

Double_t param_MMnmiss_corr[7]={
  1.07364,
  2.01487,
  45.9014,
  22.4603,
  1.06126,
  1.01389,
  0.528322};




//IM(npip+)
Double_t func_IMnpip(Double_t *x,Double_t *par)
{
  if(x[0]<1.08){
    return 1.;
  } else if(1.08 <= x[0] && x[0]<1.10){
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
  } else if(1.10 <= x[0] && x[0]<1.25){
    return exp(par[3]+par[4]*x[0]); 
  } else if(1.25 <= x[0] && x[0]<2.0){
    return par[5]*exp(-0.5*pow(((x[0]-par[6])/par[7]),2.0)); 
  } else {
    return 1.;
  }
}

Double_t param_IMnpip[8]={
     1.10885,
     1.09604,
     0.0172468,
     0.815124,
    -0.715797,
     1.23375,
     1.56849,
     0.43139
};

Double_t func_IMnpim(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 1.;
  } else if(1.07 <= x[0] && x[0]<1.10){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  } else if(1.10 <= x[0] && x[0]<1.70){
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0);
  } else {
    return 1.;
  }
}

Double_t param_IMnpim[8]={
   170945.,
    -471604.,
    433669.,
    -132922.,
    22.0792,
    -45.7144,
    33.0676,
    -8.02205
};


Double_t func_MMom(Double_t *x,Double_t *par)
{
  if(x[0]<0.4){
    return 1.0;
  }else if(0.4<= x[0] && x[0]<1.5){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return 1.0;
  }
}

/*
Double_t param_MMom[6]={
-7.57649,  
 56.2289,    
 -141.95,    
 169.847,   
-96.1952,    
 20.7417 
};*/

Double_t param_MMom[6]={
  -9.68238,
    71.0961,
    -179.177,
    213.259,
    -119.984,
    25.6942
};


Double_t func_nmom(Double_t *x,Double_t *par)
{
  if(x[0]<0.14){
    return 1.0;
  }else if(0.14<= x[0] && x[0]<1.0){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_nmom[9]={
 6.2837,
 -96.7466,
 689.093,
 -2498.12,
  5095.65,
 -6017.46,
 3971.57,
 -1278.91,
 129.626
};

Double_t func_pipmom(Double_t *x,Double_t *par)
{
  if(x[0]<0.08){
    return 1.0;
  }else if(0.08<= x[0] && x[0]<0.70){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_pipmom[9]={
 8.84121,   
 -261.15,   
 3382.08,   
-22296.8,   
 81893.6,   
 -171437,   
  197687,   
 -110674,   
 20336.8  
};

Double_t func_pimmom(Double_t *x,Double_t *par)
{
  if(x[0]<0.06){
    return 1.0;
  }else if(0.06<= x[0] && x[0]<0.73){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_pimmom[9]={
  8.02094,    
 -229.097,    
  2832.75,    
 -17372.6,    
  58512.6,    
  -110909,    
   112181,    
 -49442.4,    
  3313.79   
};


Double_t func_Mompippim(Double_t *x,Double_t *par)
{
  if(x[0]<0.98){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  }else{
    return 1.0;
  }
}

Double_t param_Mompippim[4]={
 0.98013,
 1.84325,
  -6.81408,
 5.05349
};

Double_t func_IMpippim(Double_t *x,Double_t *par)
{
  if(x[0]<0.28){
    return 1.0;
  }else if(0.28 <= x[0] && x[0] < 0.97){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0);
  }else{
    return 1.0;
  }
}

Double_t param_IMpippim[7]={
 -3.25082,
 76.8531,
 -527.225,
 1791.19,
 -3231.95,
  2959.35,
  -1073.98
};


Double_t func_IMnpipi(Double_t *x,Double_t *par)
{
  if(x[0]<1.22){
    return 1;
  }else if(1.22<=x[0] && x[0]<2.00) {
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return 1;
  }
}


Double_t param_IMnpipi[6]={
  971.716,
  -3118.7,
   4009.6,
  -2576.06,
  826.215,
  -105.722
};


