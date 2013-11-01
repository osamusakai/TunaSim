#include <admodel.h>

#include <catsim.htp>

model_data::model_data(void)
{
  rho_sd.allocate("rho_sd");
  catch_sd.allocate("catch_sd");
  sample_size.allocate("sample_size");
  effort_sd.allocate("effort_sd");
  nyrs.allocate("nyrs");
  nages.allocate("nages");
  num_fish.allocate("num_fish");
  nlint.allocate("nlint");
  effort.allocate(1,num_fish,1,nyrs,"effort");
  qchange.allocate(1,num_fish,1,nyrs,"qchange");
  obs_effort.allocate(1,num_fish,1,nyrs);
  M.allocate(1,nages,"M");
  vbcoff.allocate(1,3,"vbcoff");
  sd.allocate(1,nages,"sd");
  mean_length_sd.allocate(1,nages,"mean_length_sd");
  pop1.allocate(1,nages,"pop1");
  init_recruit.allocate(2,nyrs,"init_recruit");
  sel.allocate(1,num_fish,1,nages,"sel");
  q.allocate(1,num_fish,"q");
  xmid.allocate(1,nlint);
  J.allocate(1,nages);
 J.fill_seqadd(1,1);
  effort_devs.allocate(1,num_fish,1,nyrs);
  catch_devs.allocate(1,num_fish,1,nyrs);
  mean_length_devs.allocate(1,num_fish,1,nyrs,1,nages);
  cohort_rho_devs.allocate(1-nages,nyrs);
  q1.allocate(1,num_fish,1,nyrs,1,nages,1,nlint);
  if (ad_comm::argc > 1)
  {
    int on=0;
    if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-iseed"))>-1)
    {
      {
        iseed = atoi(ad_comm::argv[on+1]);
        random_number_generator rng(iseed);
        dvector tmp(1,3);
        tmp.fill_randu(rng);
        tmp          = 6400*(tmp-.5);
        eff_seed     = int(tmp(1));
        len_seed     = int(tmp(2));
        len_dev_seed = int(tmp(3));
        cout<<  "Currently using "   << endl<<
        adstring(ad_comm::argv[on+1])<< endl<<
        eff_seed     << endl<<
        len_seed     << endl<<
        len_dev_seed << endl<<
        " as random number seeds for sims"<<endl;
      }
    }
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 ad_comm(argc,argv), model_data() , function_minimizer(sz)
{
  initializationfunction();
  rbio.allocate(1,nyrs,"rbio");
  #ifndef NO_AD_INITIALIZE
    rbio.initialize();
  #endif
  obs_freq.allocate(1,num_fish,1,nyrs,1,nlint,"obs_freq");
  #ifndef NO_AD_INITIALIZE
    obs_freq.initialize();
  #endif
  pred_freq.allocate(1,num_fish,1,nyrs,1,nlint,"pred_freq");
  #ifndef NO_AD_INITIALIZE
    pred_freq.initialize();
  #endif
  F.allocate(1,num_fish,1,nyrs,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,nyrs,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,nyrs,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,nyrs,1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  C.allocate(1,num_fish,1,nyrs,1,nages,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  P.allocate(1,num_fish,1,nyrs,1,nages,"P");
  #ifndef NO_AD_INITIALIZE
    P.initialize();
  #endif
  obs_catch.allocate(1,num_fish,1,nyrs,"obs_catch");
  #ifndef NO_AD_INITIALIZE
    obs_catch.initialize();
  #endif
  mu.allocate(1,nyrs,1,nages,"mu");
  #ifndef NO_AD_INITIALIZE
    mu.initialize();
  #endif
  f.allocate("f");
}

void model_parameters::preliminary_calculations(void)
{
  get_mortality_and_survival_rates();
  get_numbers_at_age();
  get_catch_at_age();
  get_mean_lengths();
  get_catch_at_length();
  get_biomass_levels();  
  write_simulated_data();
  write_maundered_data();
  exit(1);
}

void model_parameters::userfunction(void)
{
}

void model_parameters::get_mortality_and_survival_rates(void)
{
  for (int fi=1;fi<=num_fish;fi++)
  {
    F(fi)=outer_prod(q(fi)*effort(fi),sel(fi));
    Z+=F(fi);
  }
  dvar_matrix TZ=trans(Z);
  for (int j=1;j<=nages;j++)
    TZ(j)+=M(j);
  Z=trans(TZ);
  S=mfexp(-1.0*Z);
}

void model_parameters::get_numbers_at_age(void)
{
  int i,j;
  for (i=2;i<=nyrs;i++)
    N(i,1)=init_recruit(i);
 
  for (j=1;j<=nages;j++)
    N(1,j)=pop1(j);
 
  for (i=1;i<nyrs;i++)
  {
    for (j=1;j<nages;j++)
      N(i+1,j+1)=N(i,j)*S(i,j);
    N(i+1,nages)+=N(i,nages)*S(i,nages);
  }
}

void model_parameters::get_catch_at_age(void)
{
  for (int fi=1;fi<=num_fish;fi++)
    C(fi)=elem_prod(elem_div(F(fi),Z),elem_prod(1.-S,N));
}

void model_parameters::get_mean_lengths(void)
{
  random_number_generator rng(len_dev_seed);
  cohort_rho_devs.fill_randn(rng);
  cohort_rho_devs*=rho_sd;
  double rho=exp(-vbcoff(3)); 
  // Mean length at age
  for (int i=1;i<=nyrs;i++)
    for (int j=1;j<=nages;j++)
    {
      rho = exp(-vbcoff(3) + cohort_rho_devs(i-j+1));
      mu(i,j) = vbcoff(1) + (vbcoff(2)-vbcoff(1))* (1.-pow(rho,J(j)-1.0))/(1.-pow(rho,nages-1.0));
    }
  //cout << mu<<endl;
  //cout << pow(mu(1),3)<<endl;
  for (int fi=1;fi<=num_fish;fi++)
  {
    mean_length_devs(fi).fill_randn(rng);
    for (int i=1;i<=nyrs;i++)
      mean_length_devs(fi,i)=elem_prod(mean_length_devs(fi,i), mean_length_sd);
  }
  // Add cohort specific changes? or annual events (fast growth in all ages w/in one year)
}

void model_parameters::get_catch_at_length(void)
{
  int i,fi;
  double cutoff=1.0;
  random_number_generator rng1(eff_seed);
  effort_devs.fill_randn(rng1);
  obs_effort=elem_prod(effort,mfexp(effort_sd*effort_devs));
  catch_devs.fill_randn(rng1);
  for (fi=1;fi<=num_fish;fi++)
    for (i=1;i<=nyrs;i++)
      obs_catch(fi,i) = sum(C(fi,i))*mfexp(catch_sd*catch_devs(fi,i));
  for (fi=1;fi<=num_fish;fi++)
    obs_effort(fi)=elem_div(obs_effort(fi),qchange(fi));
  //cout<< "True effort " <<endl<<effort<<endl<<" obs_eff: "<<endl<<obs_effort<<endl;
  //Set up length bins (Too clean??)
  l1=value(mean(trans(mu)(1))-3.0*sd(1));
  // l1=value(mu(1)-3.0*sd(1));
  //cout <<"MEAN Age 1 "<<l1<<endl;
  //ln=value(mu(nages)+3.0*sd(1));
  ln=value(mean(trans(mu)(nages))+3.0*sd(1));
  //cout <<"MEAN Age n "<<ln<<endl;
  delta=(ln-l1)/(nlint-1.0);
  xmid.fill_seqadd(l1+0.5*delta,delta);
  cout<<endl<<xmid<<endl;
  //cout<<endl<<endl<<"MeanLength "<<endl;
  for (fi=1;fi<=num_fish;fi++)
  {
    for (i=1;i<=nyrs;i++)
    {
      dvar_vector ml=mean_length_devs(fi,i)+mu(i);
      // cout<<ml<<endl;
      for (int j=1;j<=nages;j++)
      {
        // add a 5% outliers with 3sd's
        // q1 4d array (length last index), xmid vector of length bin mid-point
        q1(fi,i,j)=
        value(0.95/sd(j)*exp(-square(xmid-ml(j))/(2.0*square(sd(j))))) + value(0.05/(3.0*sd(j)) *exp(-square(xmid-ml(j))/(2.0*9.0*square(sd(j))))); 
        //=(0.95/sd*exp(-(xmid-ml)^2/(2.0*(sd^2)))) + (0.05/(3.0*sd) *exp(-(xmid-ml)^2/(2.0*9.0*(sd^2)))); 
        q1(fi,i,j)/=sum(q1(fi,i,j));
      }
    }
  }
  random_number_generator rng(len_seed);
  obs_freq.initialize();
  for (fi=1;fi<=num_fish;fi++)
  {
    for (i=1;i<=nyrs;i++)
    {
      P(fi,i)=C(fi,i)/sum(C(fi,i));
    // Vector * matrix to get annual, fishery specific l freqs
      pred_freq(fi,i)=P(fi,i)*q1(fi,i);
    // Potentially add alternative harvest sampling process (other than multinomial)
      dvector tmp(1,sample_size);
    // Fills tmp with index values corresponding to len intervals in pred_freq
      tmp.fill_multinomial(rng,value(pred_freq(fi,i)));
    // Build the histograms
      for (int ii=tmp.indexmin();ii<=tmp.indexmax();ii++)
        obs_freq(fi,i,int(tmp(ii))) +=1; 
    }  
  }
}

void model_parameters::get_biomass_levels(void)
{
  for (int i=1;i<=nyrs;i++)
    rbio(i)=N(i)*pow(mu(i),3);
}

void model_parameters::write_maundered_data(void)
{
  int i,fi;
  int nregions=1;
  // count number of fishing situations
  int num_data_sets=0;
  double cutoff=1.0;
  for (fi=1;fi<=num_fish;fi++)
    for (i=1;i<=nyrs;i++)
      if ( sum(C(fi,i)) > cutoff ) num_data_sets++;
  ofstream ofs("maunder.dat");
  ofs << "#init_matrix qindexdata(1,Ngears,StartYear,EndYear) "<<endl;
  for (i=1;i<=nyrs;i++) ofs << " 1 ";ofs<< endl; 
  for (i=1;i<=nyrs;i++) ofs << " 1 ";ofs<< endl; 
  ofs << "#init_int NRindexdata"<<endl<<"# 1"<<endl
      << "#init_matrix Rindexdata(1,NRindexdata,1,2)" <<endl<<"#1 1"<<endl;
  for (i=1;i<=nyrs+1;i++) ofs << " 1 ";ofs<< endl; 
  ofs << "#init_int Nsurveydata"<<endl<<" 1"<<endl
      << "#init_matrix surveydata(1,Nsurveydata,1,2)" <<endl<<"#survey year value cv"<<endl
      << "1 10 10 0.5"<<endl;
  ofs << "#init_int NCPUEdata"<<endl;
  ofs << "1 "<<endl;
  ofs << "#init_matrix CPUEdata(1,Ncpuedata,1,4)"<<endl;
  ofs << "#fishery year value cv"<<endl;
  ofs << "1 1 447.0393236 0.2"<<endl;
  ofs << "# init_vector Ctype(1,Ngears) // 0 = numbers, 1 = weight"<<endl; 
  ofs << " 0 0 " <<endl;
  ofs << "# init_int NC " << endl;
  ofs << " 80 "<<endl;
  ofs << "# init_matrix C(1,NC,1,5) // gear, time, landings, discard rate, cv "<<endl;
  for (i=1;i<=nyrs;i++)
    for (fi=1;fi<=num_fish;fi++)
      ofs << setfixed()<< setprecision(5)<<"   " << fi <<" "<< i << "   " << obs_catch(fi,i) << "   1 0.1 " <<endl;
  ofs << "# init_matrix E(1,Ngears, StartTime, EndTime)" <<endl;
  ofs <<  obs_effort(1) <<  endl;
  ofs <<  obs_effort(2) <<  endl;
  ofs << "# init_int NCAL"<<endl;
  ofs << "80             "<<endl;
  ofs << " # init_matrix CALin(1,NCAL,1,Nlengths+3) // gear, time, sample size, L1, L2, ...  "<<endl;
  for (i=1;i<=nyrs;i++)
    for (fi=1;fi<=num_fish;fi++)
    {
      ofs << fi << " " << i << " " << sum(obs_freq(fi,i)) << " ";
      ofs << setfixed()  << setprecision(8) << obs_freq(fi,i)/sum(obs_freq(fi,i)) << endl;
    }
}

void model_parameters::write_simulated_data(void)
{
  int i,fi;
  int nregions=1;
  // count number of fishing situations
  int num_data_sets=0;
  double cutoff=1.0;
  for (fi=1;fi<=num_fish;fi++)
    for (i=1;i<=nyrs;i++)
    {
      // cout<<obs_catch(fi,i)<<" real catch "<<sum(C(fi,i))<<endl;
      if ( sum(C(fi,i)) > cutoff ) num_data_sets++;
    }
  ofstream ofs("catsim.out");
  ofs << "# Number of number of   use generic number of     year 1" << endl;
  ofs << "#  regions   fisheries  diffusion   tag groups" << endl;
  ofs << "      " << nregions << "          " << num_fish <<  "            0 " 
      << "         0       1960      0    0    0  0  0" << endl;
  ofstream * pofs = &ofs;
  ofs << "# Regions in which fishery is located" << endl;
  for (fi=1;fi<=num_fish;fi++)
    ofs << "  1 "<< endl;
  ofs << "# Regions Area" << endl;
  ofs << "  1"  << endl;
  ofs << "# Datasets  Intervals  First  Width  Factor" << endl;
  ofs << "   " << num_data_sets << "         " <<nlint<< "     " <<l1<<"    "  << delta <<  "  1 "  << endl;
  for (i=1;i<=nyrs;i++)
  {
    for (fi=1;fi<=num_fish;fi++)
    {
      double sumcatch=sum(value(C(fi,i)));
      if (sumcatch>cutoff)
      {
        ofs << "# Year   Month  Week" << endl;
        ofs << "  " << 1959+i <<  "   1   1 " << endl;
        ofs << "# Fishery Catch Effort" << endl;
        ofs << setfixed()  << setprecision(2) << "   " << fi << "   " << obs_catch(fi,i) << "   " <<  obs_effort(fi,i) <<  endl;
        ofs << setfixed()  << setprecision(1) << obs_freq(fi,i) << endl;
      }
    }
  }
  dvar_vector relbio = rbio/max(rbio); 
  ofstream ofs1("trec.dat");
  ofstream ofs2("tbiom.dat");
  ofstream ofs3("trbiom.dat");
  ofstream ofs4("ttot.dat");
  ofstream ofs5("tinit.dat");
  ofstream ofs6("tend.dat");
  ofstream ofs7("tF1.dat");
  ofstream ofs8("tF2.dat");
  for (i=1;i<=nyrs;i++)
  {
    ofs1 << i<< " "<< setscientific()<< setprecision(5) << N(i,1)   << endl;
    ofs2 << i<< " "<< setscientific()<< setprecision(5) << rbio(i)  << endl;
    ofs3 << i<< " "<< setscientific()<< setprecision(5) << relbio(i)<< endl;
    ofs4 << i<< " "<< setscientific()<< setprecision(5) << sum(N(i))<< endl;
  }
  ofs5 << setscientific()<< setprecision(5) << N(1)     << endl;
  ofs6 << setscientific()<< setprecision(5) << N(nyrs)  << endl;
  ofs7 << setscientific()<< setprecision(5) << q(1)*effort(1,1)*sel(1)  <<" "<<  q(2)*effort(2,1)*sel(2)  << endl;
  ofs8 << setscientific()<< setprecision(5) << q(1)*effort(1,nyrs)*sel(1)  <<" "<<  q(2)*effort(2,nyrs)*sel(2)  << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(void){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  arrmblsize=400000;
  /*
  Additional simulation scenarios
    Run model w/ trends on q when sims aren't DONE
    Model with trend in q not analzyed but in data DONE
    Flat or increase in trajectory DONE
    Error in catch data DONE
    Add holes in data
    trend in misreporting w/ trend in Q's
    Changes in growth should be cumulative in time DONE
  Harder ones:
    Multi area specification wrong
    Changes in growth should be larger 
  */
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}
