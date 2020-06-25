//Missing mass



Double_t func_MMnmiss_mod(Double_t *x,Double_t *par)
{
   //connection using woods-saxon
   if(0.0<x[0] && x[0]<1.5){
  //   return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
  //     +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
  //     +par[7]*pow(x[0],7.0))*(1./(1.0+exp((x[0]-1.10)/par[8])))
     return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0))*(1./(1.0+exp((x[0]-1.08)/par[10])))
       +(1.0-1./(1.0+exp((x[0]-1.08)/par[10])))*(par[11]*exp(-0.5*pow(((x[0]-par[12])/par[13]),2.0)))*(1./(1.0+exp((x[0]-1.14)/par[14])))
       +(1.0-1./(1.0+exp((x[0]-1.14)/par[14])))*
        (par[15]+par[16]*x[0]+par[17]*pow(x[0],2.0)+par[18]*pow(x[0],3.0)+par[19]*pow(x[0],4.0));
   }else{
     return 1.0;
   }
}



Double_t param_MMnmiss_mod[20]={
0.000281332,
0.490481,
2.6641,
-15.438,
13.4367,
188.441,
-689.654,
998.37,
-662.618,
166.036,
0.018,
3.59474,
0.712467,
0.432273,
0.01,
260.568,
-803.045,
942.772,
-494.219,
96.9351
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


Double_t func_IMnpipmul_s(Double_t *x,Double_t *par)
{
  if(1.06<=x[0] && x[0]<1.92){
    double ret=  (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))*(1./(1.0+exp((x[0]-1.25)/par[11])))+
           (par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))*(1.0-1./(1.0+exp((x[0]-1.25)/par[11])));
    if(ret<0) return 0;
    else return ret;
  }else{
    return 1.;
  }
}

Double_t param_IMnpip_s[12]={
-327.854,
748.394,
-561.008,
138.25,
-3277.16,
12041.2,
-18359.5,
14891.1,
-6782.67,
1645.87,
-166.272,
0.01
};



Double_t func_IMnpip_wK0_mod(Double_t *x,Double_t *par)
{
  //pol3(1.08-1.14)+pol4(1.14-1.25)+pol5(1.25-1.93)
  if(1.08<x[0] && x[0]<1.63){
    double ret = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))
    *(1./(1.0+exp((x[0]-1.13)/par[4])))
    +(1.0-1./(1.0+exp((x[0]-1.13)/par[4])))
    *(par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)+par[9]*pow(x[0],4.0))
    *(1./(1.0+exp((x[0]-1.25)/par[10])))
    +(1.0-1./(1.0+exp((x[0]-1.25)/par[10])))
    *(par[11]+par[12]*x[0]+par[13]*pow(x[0],2.0)+par[14]*pow(x[0],3.0)+par[15]*pow(x[0],4.0));
  
    if(ret<0) return 0;
    else return ret;
  }else{
    return 0.0;
  }
};

Double_t param_IMnpip_wK0_mod[16]={
-88509,
242125,
-220810,
67133.3,
0.01,
-3776.85,
6211.08,
27.1978,
-4245.76,
1746.82,
0.01,
32.9473,
-35.4828,
-0.449256,
10.4611,
-2.73038
};




Double_t func_IMnpim_mod(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 0.0;
  }else if(1.07 <= x[0]){
    return (par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0)))
          *(1./(1.0+exp((x[0]-1.10)/par[3])))
          +(1.0 - 1./(1.0+exp((x[0]-1.10)/par[3])))
          //*(exp(par[4]+par[5]*x[0]));
          *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)
           +par[9]*pow(x[0],5.0));//+par[10]*pow(x[0],6.0)+par[11]*pow(x[0],7.0));
  }else{
   return 0;
  }
}


Double_t param_IMnpim_mod[10]={
1.99912,
1.11494,
0.0176427,
0.02,
133.165,
-257.346,
113.524,
67.3744,
-67.331,
14.6219
};


Double_t func_IMnpim_wK0_mod(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<2.0) {
   //1.07-1.11  gaus    1.11-1.24 pol5            1.24-2.0 exp
    return (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
         *(1./(1.0+exp((x[0]-1.11)/par[3])))
         +(1.0 - (1./(1.0+exp((x[0]-1.11)/par[3]))))
         //*(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0))
         *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0))
         *(1./(1.0+exp((x[0]-1.24)/par[8])))
         +(1.0 - (1./(1.0+exp((x[0]-1.24)/par[8]))))
         *exp(par[9]+par[10]*x[0]);
  } else {
    return 0.;
  }
}


Double_t param_IMnpim_wK0_mod[11]={
3.80053,
1.10253,
0.00869531,
0.02,
4266.42,
-10484.8,
8583.89,
-2340.8,
0.01,
8.57184,
-7.00872
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
126.564,
-1837.3,
9351.47,
-16281.6,
0.01,
11.039,
-61.0051,
118.523,
-49.6603,
-125.447,
165.959,
-59.1096
};


Double_t param_nmom_wK0_mod[12]={
2674.21,
-45316.8,
259551,
-496799,
0.01,
107.933,
-617.738,
1281.31,
-646.596,
-1456.37,
2235.93,
-909.037
};




Double_t func_IMpippim_mod(Double_t *x,Double_t *par)
{
  if(0.27<x[0]){
  return (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
        *(1./(1.0+exp((x[0]-0.31)/par[3])))            
      +(1.0 - 1./(1.0+exp((x[0]-0.31)/par[3])))
       *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))
       *(1./(1.0+exp((x[0]-0.73)/par[11])))
      +(1.0 - 1./(1.0+exp((x[0]-0.73)/par[11])))
       *(par[12]*exp(-0.5*pow(((x[0]-par[13])/par[14]),2.0)))
       ;
   }else{
    return 0.0;
   }
};


Double_t param_IMpippim_mod[15]={
0.431559,
0.297098,
0.0124366,
0.02,
2.1313,
-5.02899,
-6.19126,
17.1507,
50.3553,
26.175,
-154.199,
0.03,
4.39784,
0.671838,
0.088268
};




Double_t func_q_mod(Double_t *x,Double_t *par)
{
  if(0<x[0]) {  //&& x[0]<1.17){
    return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0));
  }else{
   return 0;
  }
}

Double_t param_q_mod[8]={
4.37725,
-5.26915,
1.8429,
12.0264,
-19.8566,
-5.04731,
27.5322,
-12.9654
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
0.190352,
6.57234,
-33.3275,
63.9037,
-40.3182,
-26.5057,
46.8204,
-16.7796
};


