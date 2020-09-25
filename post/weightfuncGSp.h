//Missing mass



Double_t func_MMnmiss_mod(Double_t *x,Double_t *par)
{
   //connection using woods-saxon
   if(0.0<x[0] && x[0]<1.5){
     double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0))*(1./(1.0+exp((x[0]-1.075)/par[10])))
       +(1.0-1./(1.0+exp((x[0]-1.075)/par[10])))*(par[11]*exp(-0.5*pow(((x[0]-par[12])/par[13]),2.0)))*(1./(1.0+exp((x[0]-1.14)/par[14])))
       +(1.0-1./(1.0+exp((x[0]-1.14)/par[14])))*
        (par[15]+par[16]*x[0]+par[17]*pow(x[0],2.0)+par[18]*pow(x[0],3.0)+par[19]*pow(x[0],4.0));
     if(ret<0) return 0;
     else      return ret;
   }else{
     return 1.0;
   }
}



Double_t param_MMnmiss_mod[20]={
0.00416195,
2.25759,
-2.27426,
-4.60196,
13.5346,
139.641,
-634.334,
1029.75,
-736.461,
194.433,
0.018,
36.9335,
0.340912,
0.302993,
0.01,
299.414,
-864.313,
947.832,
-465.967,
86.3392
};



Double_t func_MMnmiss_wK0_mod(Double_t *x,Double_t *par)
{
  double ret= (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0))
       *(1./(1.0+exp((x[0]-1.08)/par[6])))
       +(1. - 1./(1.0+exp((x[0]-1.08)/par[6])))
       *(par[7]+par[8]*x[0]+par[9]*pow(x[0],3.0)+par[10]*pow(x[0],4.0))
       ;
   if(ret<0) return 0;
   else return ret;
}

Double_t param_MMnmiss_wK0_mod[11]={
-0.00714053,
0.991907,
4.4545,
-20.0309,
35.2696,
-19.0939,
0.01,
-7.59543,
8.25753,
4.54261,
-4.65572
};


Double_t func_IMnpip_mod(Double_t *x,Double_t *par)
{
  if(1.06<=x[0] && x[0]<2.00){
    //gaus (1.06-1.23) 
    double ret= // (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))
                (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
                *(1./(1.0+exp((x[0]-1.23)/par[3])))
                +(1.0-1./(1.0+exp((x[0]-1.23)/par[3])))
                *(exp(par[4]+par[5]*x[0]));
    if(ret<0) return 0;
    else return ret;
  }else{
    return 1.;
  }
}

Double_t param_IMnpip_mod[6]={
4.71493,
1.16827,
0.0451448,
0.01,
13.5979,
-10.2407
};


Double_t func_Mompippim(Double_t *x,Double_t *par)
{
  double ret =  
    par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)
    +par[5]*pow(x[0],5.0)    
    +par[6]*pow(x[0],6.0);    
  if(ret<0 ) return 0;
  else return ret;
}

Double_t param_Mompippim[7]={
0.807614,
-3.16567,
30.9999,
-58.2724,
-3.15916,
68.6104,
-35.6868
};

Double_t func_Momnpip(Double_t *x,Double_t *par)
{
  double ret =  
    par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)
    +par[5]*pow(x[0],5.0);
  if(ret<0 ) return 0;
  else return ret;
}

Double_t param_Momnpip[6]={
1.53735,
3.38526,
-26.817,
51.5213,
-40.917,
11.6677
};


Double_t func_Momnpim(Double_t *x,Double_t *par)
{
  double ret =  
    par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)
    +par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0);
  if(ret<0 ) return 0;
  else return ret;
}

Double_t param_Momnpim[8]={
1.08198,
-6.6943,
72.6287,
-338.881,
780.79,
-922.573,
538.291,
-122.778
};

Double_t func_IMnpip_wK0_mod(Double_t *x,Double_t *par)
{
  if(1.08<x[0] && x[0]<1.63){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))
    *(1./(1.0+exp((x[0]-1.19)/par[4])))
    +(1.0-1./(1.0+exp((x[0]-1.19)/par[4])))
    //*(par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)+par[9]*pow(x[0],4.0))
    *(par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0))
    //*(1./(1.0+exp((x[0]-1.25)/par[10])))
    *(1./(1.0+exp((x[0]-1.25)/par[8])))
    //+(1.0-1./(1.0+exp((x[0]-1.25)/par[10])))
    +(1.0-1./(1.0+exp((x[0]-1.25)/par[8])))
    //*(par[11]+par[12]*x[0]+par[13]*pow(x[0],2.0)+par[14]*pow(x[0],3.0)+par[15]*pow(x[0],4.0));
    *(par[9]+par[10]*x[0]+par[11]*pow(x[0],2.0)+par[12]*pow(x[0],3.0)+par[13]*pow(x[0],4.0));
  
    if(ret<0) return 0;
    else return ret;
  }else{
    return 0.0;
  }
};

Double_t param_IMnpip_wK0_mod[14]={
-232.137,
179.349,
244.028,
-194.096,
0.01,
-14.2499,
24.9294,
-9.69115,
0.01,
37.9034,
-47.0971,
14.3774,
0.154323,
-0.00105421,
};




Double_t func_IMnpim_mod(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 0.0;
  }else if(1.07 <= x[0]){
    double ret = (par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0)))
          *(1./(1.0+exp((x[0]-1.12)/par[3])))
          +(1.0 - 1./(1.0+exp((x[0]-1.12)/par[3])))
          *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)
           +par[9]*pow(x[0],5.0)); 
    if(ret>0) return ret;
    else      return 0;
 }else{
   return 0;
  }
}


Double_t param_IMnpim_mod[10]={
1.15104,
1.11053,
0.0228043,
0.010,
730.135,
-2459.45,
3271.54,
-2144.66,
693.116,
-88.4319
};


Double_t func_IMnpim_wK0_mod(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<2.0) {
   //1.07-1.11  gaus    1.11-1.24 pol5            1.24-2.0 exp
    return (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
         *(1./(1.0+exp((x[0]-1.092)/par[3])))
         +(1.0 - (1./(1.0+exp((x[0]-1.092)/par[3]))))
         //*(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0))
         //*(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0))
         *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0))
         //*(1./(1.0+exp((x[0]-1.22)/par[8])))
         *(1./(1.0+exp((x[0]-1.22)/par[7])))
         //+(1.0 - (1./(1.0+exp((x[0]-1.22)/par[8]))))
         +(1.0 - (1./(1.0+exp((x[0]-1.22)/par[7]))))
         //*exp(par[9]+par[10]*x[0]);
         *exp(par[8]+par[9]*x[0]);
  } else {
    return 0.;
  }
}


Double_t param_IMnpim_wK0_mod[10]={
0.0858578,
1.10463,
0.00904288,
0.005,
424.378,
-705.599,
293.465,
0.005,
8.48006,
-6.98101
};




Double_t func_nmom_mod(Double_t *x,Double_t *par)
{
  if(x[0]<0.140){
    return 1.0;
  }else if(0.14<= x[0] && x[0]<1.0){
 //   return exp(par[0]+par[1]*x[0])*((1./(1.0+exp((x[0]-0.2)/par[2]))))
 //   + (par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0)+par[6]*pow(x[0],3.0)+par[7]*pow(x[0],4.0)+par[8]*pow(x[0],5.0)
 //   +par[9]*pow(x[0],6.0))*(1.-1./(1.0+exp((x[0]-0.2)/par[2])));
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))*((1./(1.0+exp((x[0]-0.2)/par[4]))))
    + (par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)+par[9]*pow(x[0],4.0)+par[10]*pow(x[0],5.0)
    +par[11]*pow(x[0],6.0))*(1.-1./(1.0+exp((x[0]-0.2)/par[4])));
    if(ret<0) return 0;
    else return ret;
  }else{
    return 1.0;
  }
}

Double_t param_nmom_mod[12]={
67.6781,
-804.696,
3497.47,
-5382.05,
0.01,
18.656,
-148.512,
516.165,
-973.885,
1042.07,
-596.986,
142.917
};


Double_t param_nmom_wK0_mod[12]={
2946.78,
-50380.9,
290158,
-557585,
0.01,
106.424,
-612.966,
1283.35,
-651.591,
-1466.61,
2230.77,
-891.24
};




Double_t func_IMpippim_mod(Double_t *x,Double_t *par)
{
  if(0.27<x[0]){
  double ret = (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
        *(1./(1.0+exp((x[0]-0.30)/par[3])))            
      +(1.0 - 1./(1.0+exp((x[0]-0.30)/par[3])))
       *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))
       *(1./(1.0+exp((x[0]-0.73)/par[11])))
      +(1.0 - 1./(1.0+exp((x[0]-0.73)/par[11])))
       *(par[12]*exp(-0.5*pow(((x[0]-par[13])/par[14]),2.0)))
       ;
  if(ret<0) return 0;
  else      return ret;
  }else{
    return 0.0;
  }
};


Double_t param_IMpippim_mod[15]={
0.347053,
0.296448,
0.00910623,
0.015,
5.41219,
-28.1255,
33.0591,
76.4173,
-64.7074,
-248.039,
263.232,
0.03,
0.644541,
0.785028,
0.0522259
};




Double_t func_q_mod(Double_t *x,Double_t *par)
{
  if(0<x[0]) {  //&& x[0]<1.17){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0));
    if(ret>0) return ret;
    else      return 0;
  }else{
   return 0;
  }
}

Double_t param_q_mod[8]={
3.0229,
-2.85293,
3.70457,
-2.08965,
-2.97093,
1.50674,
2.97921,
-1.88548
};

Double_t func_q_wK0_mod(Double_t *x,Double_t *par)
{
  if(0<x[0] && x[0]<1.3) {  //&& x[0]<1.17){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0));
    
    
    if(ret <0 ) return 0;
    else return ret;
  }else{
    return 0;
  }
}



Double_t param_q_wK0_mod[8]={
0.287913,
6.54975,
-34.3679,
64.8546,
-38.1852,
-29.9088,
47.7638,
-16.515
};




