08003765432
082
DATA_SECTION
  // add time trends to catchability to catchsim2.tpl
  int iseed
  int eff_seed
  int len_seed
  int len_dev_seed
  init_number rho_sd   
  init_number catch_sd   
  init_number sample_size
  init_number effort_sd
  init_int nyrs  
  init_int nages
  init_int num_fish
  init_int nlint;
  init_matrix effort(1,num_fish,1,nyrs)
  init_matrix qchange(1,num_fish,1,nyrs)
  matrix obs_effort(1,num_fish,1,nyrs)
  init_vector M(1,nages)
  init_vector vbcoff(1,3)
  init_vector sd(1,nages)
  init_vector mean_length_sd(1,nages)
  init_vector pop1(1,nages);
  init_vector init_recruit(2,nyrs)
  init_matrix sel(1,num_fish,1,nages)
  init_vector q(1,num_fish)
  vector xmid(1,nlint)
  vector J(1,nages)
 !! J.fill_seqadd(1,1);
  int l1;
  int ln;
  matrix effort_devs(1,num_fish,1,nyrs)
  matrix catch_devs(1,num_fish,1,nyrs)
  3darray mean_length_devs(1,num_fish,1,nyrs,1,nages)
  vector cohort_rho_devs(1-nages,nyrs)
  4darray q1(1,num_fish,1,nyrs,1,nages,1,nlint)
  number delta
 LOCAL_CALCS
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
 END_CALCS
 
PARAMETER_SECTION
  vector rbio(1,nyrs)
  3darray obs_freq(1,num_fish,1,nyrs,1,nlint)
  3darray pred_freq(1,num_fish,1,nyrs,1,nlint)
  3darray F(1,num_fish,1,nyrs,1,nages)
  matrix Z(1,nyrs,1,nages)
  matrix S(1,nyrs,1,nages)
  matrix N(1,nyrs,1,nages)
  3darray C(1,num_fish,1,nyrs,1,nages)
  3darray P(1,num_fish,1,nyrs,1,nages)
  matrix obs_catch(1,num_fish,1,nyrs)
  matrix mu(1,nyrs,1,nages)
  objective_function_value f
    
PRELIMINARY_CALCS_SECTION
  get_mortality_and_survival_rates();
  get_numbers_at_age();
  get_catch_at_age();
  get_mean_lengths();
  get_catch_at_length();
  get_biomass_levels();  
  write_simulated_data();

  write_maundered_data();
  exit(1);
PROCEDURE_SECTION
FUNCTION get_mortality_and_survival_rates
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

FUNCTION get_numbers_at_age
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

FUNCTION get_catch_at_age
  for (int fi=1;fi<=num_fish;fi++)
    C(fi)=elem_prod(elem_div(F(fi),Z),elem_prod(1.-S,N));

FUNCTION get_mean_lengths
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

FUNCTION get_catch_at_length
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

FUNCTION get_biomass_levels  
  for (int i=1;i<=nyrs;i++)
    rbio(i)=N(i)*pow(mu(i),3);

FUNCTION write_maundered_data
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
FUNCTION write_simulated_data
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

TOP_OF_MAIN_SECTION
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
