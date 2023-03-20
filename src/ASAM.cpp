// Age Structured Assessment Model 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ========== Data Inputs ============================
  
  // Indices
  DATA_INTEGER(Nyear);                              // total number of years
  DATA_INTEGER(Nage);                               // number of ages
  
  // Data in likelihood
  DATA_VECTOR(harv_t);                              // observed total fishery harvest
  DATA_VECTOR(I_t);                                 // observed CPUE for each year
  DATA_MATRIX(agecomp_ta);                          // observed age composition
  DATA_MATRIX(agecomp_FI_ta);                       // observed age composition FI
  DATA_VECTOR(SS);                                  // sample size
  DATA_VECTOR(SS_FI);                               // sample size FI
  
  // Known values
  DATA_VECTOR(ages);                                // ages
  DATA_VECTOR(L_a);                                 // length-at-age
  DATA_VECTOR(W_a);                                 // weight-at-age
  DATA_VECTOR(M_vec);                               // natural mortality-at-age
  matrix<Type> M_ta(Nyear,Nage);                    // natural mortality matrix
  DATA_VECTOR(Fec_a);                               // fecundity-at-age
  DATA_INTEGER(Amin);                               // minimum age
  DATA_INTEGER(Amax);                               // maximum age
  
  
  // ========== Parameters ============================
  
  // Estimated parameters
  // PARAMETER(dummy); // for testing
  
  PARAMETER(h);                           // steepness 
  PARAMETER(log_R0);                      // unfished recruitment
  
  PARAMETER(log_q);                       // catchability
  PARAMETER(log_FI_q);                    // FI catchability
  
  PARAMETER(log_SigR);                    // recruitment variation
  PARAMETER(log_SigC);                    // observation error - catch
  PARAMETER(log_CV_I);                    // observation error - index
  PARAMETER(log_theta);                   // theta - dirichlet-multinomial
  PARAMETER(log_theta_FI);                // theta - dirichlet-multinomial FI
  
  PARAMETER_VECTOR(log_Fint);             // fishing mortality    
  PARAMETER_VECTOR(log_recruit_devs);     // recruitment deviations
  
  // selectivity
  PARAMETER(p1);                          // selectivity - peak: beginning size for plateau
  PARAMETER(p2);                          // selectivity - top: width of plateau, as logistic between peak and max size
  PARAMETER(p3);                          // selectivity - ascending width
  PARAMETER(p4);                          // selectivity - descending width
  PARAMETER(p5);                          // selectivity - initial: selectivity at first bin
  PARAMETER(p6);                          // selectivity - final: selectivity at last bin
  PARAMETER(pFI_50);                      // FI selectivity
  PARAMETER(pFI_slope);                   // FI selectivity
  
  
  // ========== Global values ============================
  
  using namespace density;
  Type JNLL = 0.0;
  vector<Type> JNLL_all(5); // 5 likelihoods (recruitment, CPUE, age comp, agecomp FI)
  JNLL_all.setZero();
  
  int a; // age
  int y; // year
  
  // ========== Transform parameters ============================
  
  Type R0 = exp(log_R0);
  
  // error terms
  Type SigR = exp(log_SigR);
  Type SigC = exp(log_SigC);
  Type CV_I = exp(log_CV_I);
  Type theta = exp(log_theta);
  Type theta_FI = exp(log_theta_FI);
  
  Type q = exp(log_q);
  Type FI_q = exp(log_FI_q);
  
  
  // ========== Derived parameters ============================
  
  vector<Type> sel_a(Nage);                                       // fishery selectivity at age
  sel_a.setZero();
  matrix<Type> sel_ta(Nyear,Nage);                                // selectivity per year and age
  vector<Type> FI_survey_sel(Nage);                               // selectivity for FI survey
  Type peak = 0;                                                  // end point where selectivity = 1.0
  vector<Type> jf1(Nage);                                         // joiner functions for asc
  vector<Type> jf2(Nage);                                         // joiner functions for dsc    
  Type t1;                                                        // t1
  Type t2;                                                        // t2
  vector<Type> asc(Nage);                                         // ascending portions
  vector<Type> dsc(Nage);                                         // descending portions
  
  vector<Type> q_t(Nyear);                                        // catchability per year
  matrix<Type> FM_ta(Nyear,Nage);                                 // instantaneous fishing mortality
  
  matrix<Type> N_ta(Nyear,Nage);                                  // Predicted abundance at age
  vector<Type> N_t(Nyear);                                        // Predicted abundance
  matrix<Type> Nbar(Nyear,Nage);                                  // Nbar for pred fishery CPUE
  
  matrix<Type> CPUE_index_ta(Nyear,Nage);                         // Predicted catch per unit of effort (abundance index)
  vector<Type> CPUE_index(Nyear);                                 // Predicted catch per unit of effort (abundance index)
  matrix<Type> CAA(Nyear,Nage);                                   // Predicted fishery catch at age
  matrix<Type> CAA_FI(Nyear,Nage);                                // Predicted fishery catch at age FI
  vector<Type> Totcatch(Nyear);                                   // Predicted total catch
  vector<Type> Totcatch_FI(Nyear);                                // Predicted total catch FI
  matrix<Type> Agecomp(Nyear,Nage);                               // Predicted age composition
  matrix<Type> Agecomp_FI(Nyear,Nage);                            // Predicted age composition for FI (bottom longline)
  matrix<Type> Harv_ta(Nyear,Nage);                               // Predicted harvest weight (kg)
  vector<Type> Harv(Nyear);                                       // Predicted harvest weight (kg)

  vector<Type> R_t(Nyear);                                        // recruitment
  matrix<Type> SSB_ta(Nyear,Nage);                                // Spawning stock biomass by year and age
  matrix<Type> VB_ta(Nyear,Nage);                                 // Vulnerable biomass by year and age
  vector<Type> SSB(Nyear);                                        // Spawning stock biomass by year
  vector<Type> VB(Nyear);                                         // Vulnerable biomass by year
  vector<Type> Depletion(Nyear);                                  // Depletion (SSB/SSB0)
  matrix<Type> TB_ta(Nyear,Nage);                                 // biomass by year and age
  vector<Type> TB(Nyear);                                         // biomass by year
  vector<Type> N0_a(Nage);                                        // Unfished numbers at age
  vector<Type> lxo(Nage);                                         // survivorship (for calculation of N0_a)
  Type SSB0 = 0.0;                                                // Unfished spawning biomass
  
  // age composition drichlet-multinomial
  matrix<Type> DM1(Nyear,Nage);
  matrix<Type> DM2(Nyear,Nage);
  matrix<Type> DM1_FI(Nyear,Nage);
  matrix<Type> DM2_FI(Nyear,Nage);
  
  
  
  // ========== Calculations ============================
  
  // Selectivity //
  peak = p1+Type(1.0)+((Type(0.99)*Amax-p1-Type(1.0))/(Type(1.0)+exp(-p2)));
  for(a=0; a<Nage; a++){
    jf1(a) = pow((1+exp(-20*((int(a)-p1)/(1+fabs(int(a)-p1))))),-1);
    jf2(a) = pow((1+exp(-20*((int(a)-peak)/(1+fabs(int(a)-peak))))),-1);
    t1 = exp(-pow(Amin-p1,2)/exp(p3));
    t2 = exp(-pow(Amax-peak,2)/exp(p4));
    asc(a) = pow(1+exp(-p5),-1)+(1-pow(1+exp(-p5),-1))*((exp(-pow(int(a)-p1,2)/exp(p3))-t1)/(1-t1));
    dsc(a) = 1+(pow(1+exp(-p6),-1)-1)*((exp(-(int(a)-peak)/exp(p4))-1)/(t2-1));
    // final selectivity at length
    sel_a(a) = asc(a)*(1-jf1(a))+jf1(a)*((1-jf2(a))+jf2(a)*dsc(a));
  }
  
  for(a=0; a<Nage; a++){
    FI_survey_sel(a) = Type(1.0) / (Type(1.0)+exp(-(int(a)-pFI_50)/pFI_slope));
  }
  FI_survey_sel = FI_survey_sel/max(FI_survey_sel);

  for(y=0; y<Nyear; y++){
    for(a=0; a<Nage; a++){
      sel_ta(y,a) = sel_a(a);
    }
  }

  // Mortalities //
  for(y=0; y<Nyear; y++){
    for(a=0; a<Nage; a++){
      FM_ta(y,a) = sel_ta(y,a)*exp(log_Fint(y)); // setting year and age specific fishing mortality as the product of selectivity and year specific fishing intensity
    }
  }
  for(y=0; y<Nyear; y++){
    for(a=0; a<Nage; a++){
      M_ta(y,a) = M_vec(a);
    }
  }



//========== Population Dynamics ============================

  // Unfished SSB //
  lxo(0) = Type(1.0);
  for(a=1; a<Nage; a++){
    lxo(a) = lxo(a-1)*exp(-M_ta(0,a-1));
  }
  lxo(Nage-1) = lxo(Nage-1) / (1-exp(-M_ta(0,Nage-1)));
  N0_a = R0*lxo;
  SSB0 = (N0_a*Fec_a).sum();

  // Initialize //
  R_t(0) = R0;

  // abundance at age and spawning biomass first year
  for(a=0; a<Nage; a++){
    if(a==0){
      N_ta(0,a) = R_t(0)*lxo(a)*exp(log_recruit_devs(a));
    }
    if((a>=1) & (a<(Nage-1))){
      N_ta(0,a) = N_ta(0,a-1)*exp(-M_ta(0,a-1)-FM_ta(0,a-1));
    }
    if(a==(Nage-1)){
      N_ta(0,a) = (N_ta(0,a-1)*exp(-M_ta(0,a-1)-FM_ta(0,a-1))) / (1-exp(-M_ta(0,a-1)-FM_ta(0,a-1)));
    }

    SSB_ta(0,a) = N_ta(0,a)*Fec_a(a);
    TB_ta(0,a) = N_ta(0,a)*W_a(a);
    VB_ta(0,a) = N_ta(0,a)*W_a(a)*sel_a(a);
    if(a>0){
      N_t(0) += N_ta(0,a);
      SSB(0) += SSB_ta(0,a);
      TB(0) += TB_ta(0,a);
      VB(0) += TB_ta(0,a);
    }
  }


  // Population loop //
  for(y=1; y<Nyear; y++){
    R_t(y) = ((4*h*R0*SSB(y-1))) / (SSB0*(1-h)+SSB(y-1)*(5*h-1)) * exp(log_recruit_devs(y));

    // age structured dynamics
    for(a=0; a<Nage; a++){
      if((y>=1) & (a==0)){
        N_ta(y,a) = R_t(y);
      }
      if((y>=1) & (a>=1) & (a<(Nage-1))){
        N_ta(y,a) = N_ta(y-1,a-1)*exp(-M_ta(y-1,a-1)-FM_ta(y-1,a-1));
      }
      if((y>=1) & (a==(Nage-1))){
        N_ta(y,a) = (N_ta(y-1,a-1)*exp(-M_ta(y-1,a-1)-FM_ta(y-1,a-1))) + (N_ta(y-1,a)*exp(-M_ta(y-1,a)-FM_ta(y-1,a)));
      }

      SSB_ta(y,a) = N_ta(y,a)*Fec_a(a);
      TB_ta(y,a) = N_ta(y,a)*W_a(a);
      VB_ta(y,a) = N_ta(y,a)*W_a(a)*sel_a(a);

      if(a>0){
        N_t(y) += N_ta(y,a);
        SSB(y) += SSB_ta(y,a);
        TB(y) += TB_ta(y,a);
        VB(y) += VB_ta(y,a);
      }
    }
    // Depletion
    Depletion(y) = SSB(y)/SSB0;
  }



  // ========== Fishery data ============================
  for(y=0; y<Nyear; y++){
    for(a=0; a<Nage; a++){
      CAA(y,a) = N_ta(y,a) * (Type(1.0)-exp(-M_ta(y,a)-FM_ta(y,a))) * (FM_ta(y,a))/(M_ta(y,a)+FM_ta(y,a));
      CAA_FI(y,a) = N_ta(y,a)*FI_q*FI_survey_sel(a);
      CPUE_index_ta(y,a) = TB_ta(y,a)*sel_ta(y,a)*q;
      Harv_ta(y,a) = CAA(y,a)*W_a(a);
    }
    Totcatch(y) = CAA.row(y).sum();
    Totcatch_FI(y) = CAA_FI.row(y).sum();
    Agecomp.row(y) = CAA.row(y)/(Totcatch(y)+Type(1e-8));
    Agecomp_FI.row(y) = CAA_FI.row(y)/(Totcatch_FI(y)+Type(1e-8));
    CPUE_index(y) = CPUE_index_ta.row(y).sum();
    Harv(y) = Harv_ta.row(y).sum();
  }



  // ========== Likelihoods ============================

  // Recruitment //
  JNLL_all(0) = Type(0.0);
  for(a=0; a<Nyear; a++){
    JNLL_all(0) -= dnorm(log_recruit_devs(a), Type(0.0), SigR, true);
  }
  
  // CPUE index //
  JNLL_all(1) = Type(0.0);
  for(y=0; y<Nyear; y++){
    JNLL_all(1) += log(CV_I*CPUE_index(y)+1e-3)+0.5*pow((I_t(y)-CPUE_index(y))/(CV_I*CPUE_index(y)+1e-3), 2);
  }

  // Age comp //
  JNLL_all(2) = Type(0.0);
  vector<Type> Beta(Nyear);
  Beta = theta*SS;
  for(y=0; y<Nyear; y++){
    // if ages observed, calculate likelihood
    if(SS(y)>0){
      for(a=0; a<Nage; a++){
        DM1(y,a) += lgamma(SS(y)*agecomp_ta(y,a)+1);
        DM2(y,a) += lgamma(SS(y)*agecomp_ta(y,a) + Beta(y)*Agecomp(y,a)) - lgamma(Beta(y)*Agecomp(y,a));
      }
      JNLL_all(2) -= lgamma(SS(y)+1)-DM1.row(y).sum() + lgamma(Beta(y)) - lgamma(SS(y)+Beta(y)) + DM2.row(y).sum();
    }
  }
  
  // Age comp FI //
  JNLL_all(3) = Type(0.0);
  vector<Type> Beta_FI(Nyear);
  Beta_FI = theta_FI*SS_FI;
  for(y=0; y<Nyear; y++){
    // if ages observed, calculate likelihood
    if(SS_FI(y)>0){
      for(a=0; a<Nage; a++){
        DM1_FI(y,a) += lgamma(SS_FI(y)*agecomp_FI_ta(y,a)+1);
        DM2_FI(y,a) += lgamma(SS_FI(y)*agecomp_FI_ta(y,a) + Beta_FI(y)*Agecomp_FI(y,a)) - lgamma(Beta_FI(y)*Agecomp_FI(y,a));
      }
      JNLL_all(3) -= lgamma(SS_FI(y)+1)-DM1_FI.row(y).sum() + lgamma(Beta_FI(y)) - lgamma(SS_FI(y)+Beta_FI(y)) + DM2_FI.row(y).sum();
    }
  }

  // Harvest //
  JNLL_all(4) = Type(0.0);
  for(y=0; y<Nyear; y++){
    JNLL_all(4) += -dnorm(log(harv_t(y)+1e-3), log(Harv(y)+1e3), SigC, true);
  }

  // Total likelihood
  // JNLL = dummy*dummy; // for testing
  JNLL = sum(JNLL_all);
  
  // std::cout << "Test NLL" << JNLL_all << std::endl;
  
  
  
  // ========== Report Section ============================
  
  ADREPORT(Depletion);
  ADREPORT(SSB);
  ADREPORT(CPUE_index);
  ADREPORT(Agecomp);
  ADREPORT(Agecomp_FI);
  ADREPORT(Harv);
  
  REPORT(Depletion);
  REPORT(sel_a);
  REPORT(FI_survey_sel);
  REPORT(SSB);
  REPORT(VB);
  REPORT(SSB0);
  REPORT(CPUE_index);
  REPORT(Agecomp);
  REPORT(Agecomp_FI);
  REPORT(Harv);
  REPORT(N_ta);
  
  REPORT(JNLL_all);
  REPORT(JNLL);
  
  return JNLL;
}