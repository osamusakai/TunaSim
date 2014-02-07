DATA_SECTION
  // add time trends to catchability to catchsim2.tpl
  // modify the output format for the recent mfcl (.frq file). 
  // add equation for recruitment calculation using Beverton-Holt stock-recruitment relationship
  // add maturity schedule(knife edge)
  // add SS3-like selectivity options
  int iseed
  int eff_seed
  int len_seed
  int len_dev_seed
  int rec_dev_seed   // OS add
  init_number rho_sd   
  init_number catch_sd   
  init_number sample_size
  init_number effort_sd
  init_int catch_mth // OS temporaly add 
  init_int catch_wek // OS temporaly add
  init_int nyrs  
  init_int nages
  init_vector mage(1,2)  // OS added as the mim maturity age
  init_int num_fish
  init_int nlint;
  init_matrix effort(1,num_fish,1,nyrs)
  init_matrix qchange(1,num_fish,1,nyrs)
  matrix obs_effort(1,num_fish,1,nyrs)
  init_vector M(1,nages)
  init_vector vbcoff(1,3)
  init_vector sd(1,nages)
  init_vector mean_length_sd(1,nages)
  init_vector pop1(1,nages)
  init_int recruit_sw // OS add: recruitment input switch
  init_vector init_recruit(2,nyrs)
  init_number lnR0  // OS add (num)
  init_number B0  // OS add (Mt)
  init_number h   // OS add
  init_number v   // OS add
  init_number sigmaR  // OS add
  //init_matrix sel(1,num_fish,1,nages)
  init_vector sel_sw(1,num_fish) // OS add: sel input switch
  init_matrix selpl(1,num_fish,1,8)  // OS add
  init_vector q(1,num_fish)
  init_vector lwcoff(1,2) // OS add: W=aL^b
  vector xmid(1,nlint)
  vector J(1,nages)
 !! J.fill_seqadd(1,1);
  int l1;
  int ln;
  matrix effort_devs(1,num_fish,1,nyrs)
  matrix catch_devs(1,num_fish,1,nyrs)
  3darray mean_length_devs(1,num_fish,1,nyrs,1,nages)
  vector cohort_rho_devs(1-nages,nyrs)
  vector rec_devs(1,nyrs)      // OS add for BH
  4darray q1(1,num_fish,1,nyrs,1,nages,1,nlint)
  number alpha // OS add for BH
  number beta  // OS add for BH
  number delta 
  //number nlint // OS modified
  matrix selplt(1,num_fish,1,8)  // OS add for sel
  vector age(1,nages) // OS added
 LOCAL_CALCS
  if (ad_comm::argc > 1)  {
    int on=0;
    if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-iseed"))>-1)
    {
      {
        iseed = atoi(ad_comm::argv[on+1]);
        random_number_generator rng(iseed);
        dvector tmp(1,4);
        tmp.fill_randu(rng);
        tmp          = 6400*(tmp-.5);
        eff_seed     = int(tmp(1));
        len_seed     = int(tmp(2));
        len_dev_seed = int(tmp(3));
        rec_dev_seed = int(tmp(4));    // OS add for BH
        cout<<  "Currently using "   << endl<<
        adstring(ad_comm::argv[on+1])<< endl<<
        eff_seed     << endl<<
        len_seed     << endl<<
        len_dev_seed << endl<<
        rec_dev_seed << endl<<        // OS add for BH
        " as random number seeds for sims"<<endl;
      }
    }
  }

  for (int i=1;i<=nages;i++)
    {
    age(i) = i-1;
    }  
 END_CALCS
 
PARAMETER_SECTION
  matrix sel(1,num_fish,1,nages) // OS added
  vector mature(1,nages)       // OS added
  vector SSB(1,nyrs)           // OS add for BH (Mt)
  vector rbio(1,nyrs)
  3darray obs_freq(1,num_fish,1,nyrs,1,nlint)
  3darray pred_freq(1,num_fish,1,nyrs,1,nlint)
  3darray F(1,num_fish,1,nyrs,1,nages)
  matrix Z(1,nyrs,1,nages)
  matrix S(1,nyrs,1,nages)
  matrix N(1,nyrs,1,nages)
  3darray C(1,num_fish,1,nyrs,1,nages)
  3darray Cbio(1,num_fish,1,nyrs,1,nages)  // OS added
  3darray P(1,num_fish,1,nyrs,1,nages)
  matrix obs_catch(1,num_fish,1,nyrs)
  matrix obs_catbio(1,num_fish,1,nyrs)  // OS added
  matrix mu(1,nyrs,1,nages)
  vector flselpeak1(1,num_fish);
  vector flselpeak2(1,num_fish);
  objective_function_value f
    
PRELIMINARY_CALCS_SECTION
  get_maturity();
  get_selectivity();
  get_mean_lengths();
  get_mortality_and_survival_rates();
  get_numbers_at_age();
  get_catch_at_age();
  get_catch_at_length();
  get_biomass_levels();  
  write_simulated_data();
  write_maundered_data();
  write_simulated_data_for_SS3();
  exit(1);

PROCEDURE_SECTION
FUNCTION get_maturity
  for (int i=1;i<=nages;i++) {
    mature(i) = 1/(1+exp(mage(2)*(age(i)-mage(1))));
    }

FUNCTION get_selectivity
  for (int fi=1;fi<=num_fish;fi++) {
    flselpeak1(fi) = vbcoff(1) + (vbcoff(2)-vbcoff(1))* (1-pow(exp(-vbcoff(3)),selpl(fi,1)-1))/(1-pow( exp(-vbcoff(3)) ,nages-1));
    flselpeak2(fi) = vbcoff(1) + (vbcoff(2)-vbcoff(1))* (1-pow(exp(-vbcoff(3)),selpl(fi,1)+selpl(fi,8)-1))/(1-pow( exp(-vbcoff(3)) ,nages-1));
   // logistic (N_pal=2)
    if(sel_sw(fi)==1) {
      for (int i=1;i<=nages;i++) {
        sel(fi,i) = 1/(1+mfexp((-1)*log(19)*(age(i)-selpl(fi,1))/selpl(fi,2)));
      }
    }
   // double normal (N_pal=4)
    else if(sel_sw(fi)==22) {
      for (int i=1;i<=nages;i++) {
        if(age(i) < selpl(fi,1)) {
          sel(fi,i) = (1/sqrt(2*3.1415*pow(exp(selpl(fi,3)/2),2))*exp(((-1)*pow(age(i)-selpl(fi,1),2))/(2*pow(exp(selpl(fi,3)/2),2))))/(1/sqrt(2*3.1415*pow(exp(selpl(fi,3)/2),2)));
        }
        else if(age(i) <= selplt(fi,1)) {
          sel(fi,i) = 1;
        }
        else {
          selplt(fi,1)= selpl(fi,1)+(0.99*nages-selpl(fi,1))*1/(1+exp((-1)*selpl(fi,2)));
          sel(fi,i) = (1/sqrt(2*3.1415*pow(exp(selpl(fi,4)/2),2))*exp(((-1)*pow(age(i)-selplt(fi,1),2))/(2*pow(exp(selpl(fi,3)/2),2))))/(1/sqrt(2*3.1415*pow(exp(selpl(fi,4)/2),2)));
        }
      }
    }
   // double logistic (N_pal=8)
    else if(sel_sw(fi)==8) {
      selplt(fi,1)= 1+(1/(1+exp(selpl(fi,3))))*(selpl(fi,1)-1); // 
      selplt(fi,2)= 1/(1+exp((-1)*selpl(fi,4)*(1-selplt(fi,1))))*0.9999; 
      selplt(fi,3)= 1/(1+exp((-1)*selpl(fi,4)*(selpl(fi,1)-selplt(fi,1))))*1.00001;
      selplt(fi,4)= log10(0.5)/log10((0.5-selplt(fi,2))/(selplt(fi,3)-selplt(fi,2)));
      selplt(fi,5)= (selpl(fi,1)+selpl(fi,8))+(1/(1+exp((-1)*selpl(fi,6))))*(nages-(selpl(fi,1)+selpl(fi,8)));
      selplt(fi,6)= 1/(1+exp((-1)*selpl(fi,7)*(selpl(fi,1)+selpl(fi,8)-selplt(fi,5))))*0.9999;
      selplt(fi,7)= 1/(1+exp((-1)*selpl(fi,7)*(nages-selplt(fi,5))))*1.00001;
      selplt(fi,8)= log10(0.5)/log10((0.5-selplt(fi,6))/(selplt(fi,7)-selplt(fi,6)));
      
      for (int i=1;i<=nages;i++) {
        if(age(i) < selpl(fi,1)) {
          sel(fi,i) = selpl(fi,2)+(1-selpl(fi,2))*pow(((1/(1+exp((-1)*selpl(fi,4)*(age(i)-selplt(fi,1)))))-selplt(fi,2))/(selplt(fi,3)-selplt(fi,2)),selplt(fi,4));
        }
        else if(age(i) <= (selpl(fi,1)+selpl(fi,8))) {
          sel(fi,i) = age(i)*0+1;
        }
        else {
          sel(fi,i) = 1+((1/(1+exp((-1)*selpl(fi,5))))-1)*pow(((1/(1+exp((-1)*selpl(fi,7)*(age(i)-selplt(fi,5)))))-selplt(fi,6))/(selplt(fi,7)-selplt(fi,6)), selplt(fi,8));
        }
      }
    }
  }
  ///

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
  // Rec is direct input (base)
  if (recruit_sw==1)
     {
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
  // Rec is based on BH (OS add)
  if (recruit_sw==2)
     {
      alpha = (4*h*exp(lnR0))/(5*h-1);  // R0(num) <- exp(lnR0)
      beta = B0*(1-h)/(5*h-1);          // B0(Mt)

      for (j=1;j<=nages;j++) {
          N(1,j)=pop1(j);
          SSB(1) += pop1(j)*mature(j)*pow(mu(1,j),lwcoff(2))*lwcoff(1)/1000;
          }
      //for (j=mage;j<=nages;j++)
      //    SSB(1) += pop1(j)*pow(mu(1,j),lwcoff(2))*lwcoff(1)/1000;

      for (i=1;i<=nyrs-1;i++)
          {
          random_number_generator rng(rec_dev_seed);
          rec_devs.fill_randn(rng);
          rec_devs*=sigmaR;
          }
      for (i=1;i<nyrs;i++)
           {
            //SSB(i)= 0;
            for (j=1;j<nages;j++)
                 {
                  N(i+1,j+1)=N(i,j)*S(i,j);
                 }
            N(i+1,nages)+=N(i,nages)*S(i,nages);
            for (j=1;j<=nages;j++)
                 SSB(i+1) +=N(i+1,j)*mature(j)*pow(mu(i+1,j),lwcoff(2))*lwcoff(1)/1000;  // SSB(Mt)
            //for (j=mage;j<=nages;j++)
            //     SSB(i+1) +=N(i,j)*S(i,j)*pow(mu(i,j),lwcoff(2))*lwcoff(1)/1000;  // SSB(Mt)

            N(i+1,1)=(alpha*SSB(i+1))/(beta+SSB(i+1))*exp(rec_devs(i)+pow((-0.5)*sigmaR,2))*(1-exp(log(0.5)*SSB(i+1)/(v*B0)));
           }
     }


FUNCTION get_catch_at_age
  for (int fi=1;fi<=num_fish;fi++)
    C(fi)=elem_prod(elem_div(F(fi),Z),elem_prod(1.-S,N));

  for (int fi=1;fi<=num_fish;fi++) {
    for (int i=1;i<=nyrs;i++)
      for (int j=1;j<=nages;j++)
        Cbio(fi,i,j) = C(fi,i,j)*pow(mu(i,j),lwcoff(2))*lwcoff(1);  //Catch biomass; L(cm)-W(kg) relationsip like SBT (OS modify) ???j???
    }

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
      mu(i,j) = vbcoff(1) + (vbcoff(2)-vbcoff(1))* (1.-pow(rho,age(j)))/(1.-pow(rho,nages-1.0));
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
    {
    for (i=1;i<=nyrs;i++)
      {
      obs_catch(fi,i) = sum(C(fi,i))*mfexp(catch_sd*catch_devs(fi,i));
      obs_catbio(fi,i) = sum(Cbio(fi,i))*mfexp(catch_sd*catch_devs(fi,i)); // OS add for catch weight (kg)
      }
    }

  for (fi=1;fi<=num_fish;fi++)
    obs_effort(fi)=elem_div(obs_effort(fi),qchange(fi));

  //cout<< "True effort " <<endl<<effort<<endl<<" obs_eff: "<<endl<<obs_effort<<endl;

  //Set up length bins (Too clean??)
  l1=ceil(value(mean(trans(mu)(1))-3.0*sd(1))); // OS add "ceil"
  // l1=value(mu(1)-3.0*sd(1));
  //cout <<"MEAN Age 1 "<<l1<<endl;
  //ln=value(mu(nages)+3.0*sd(1));
  ln=ceil(value(mean(trans(mu)(nages))+3.0*sd(1))); // OS add "ceil"
  //cout <<"MEAN Age n "<<ln<<endl;
  delta=ceil((ln-l1)/(nlint-1.0));  // OS add "ceil"
  xmid.fill_seqadd(l1+0.5*delta,delta);  cout<<endl<<xmid<<endl;
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
    rbio(i)=N(i)*pow(mu(i),lwcoff(2))*lwcoff(1); //L(cm)-W(kg) relationsip like SBT (OS modify)

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
  int i,j,fi;
  // ignoring spatial structure ,  OS add "j"   
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
  ofstream ofs("tunasim_out_forMFCL.frq");
  ofs << "# Simulated data for MULTIFAN-CL analysis" << endl;
  ofs << "#" << endl;
  ofs << "# Number of number of   use generic number of     year1  a b c d e" << endl;
  ofs << "#  regions   fisheries  diffusion   tag groups" << endl;
  ofs << "      " << nregions << "          " << num_fish <<  "            1 " 
      << "         0         1960  0  0  1  0  4" << endl;    // OS modified
  ofstream * pofs = &ofs;
  ofs << "# Relative region sizes" << endl; // OS added
  ofs << "  1"  << endl;                    // OS added
  ofs << "# Regions in which fishery is located" << endl;
  for (fi=1;fi<=num_fish;fi++)
    ofs << "  1 "<< endl;
  ofs << "# Incidence matrix (not vestige of a movement matrix in this simulation)" << endl;       // OS added
  ofs << "#  1"  << endl;                    // OS added
  // ofs << "# Regions Area" << endl;        // OS deleted
  // ofs << "  1"  << endl;                  // OS deleted
  ofs << "# Data flags (1st line: 0=catch in number; 1=catch in weight)" << endl;
  for (i=1;i<=5;i++)               // OS added
  {
    for (fi=1;fi<=num_fish;fi++)             // OS added
        ofs << "  0 "<< endl;                // OS added
  }
  ofs << "# num movement" << endl;           // OS added
  ofs << "#  periods" << endl;               // OS added
  ofs << "  1"  << endl;                     // OS added
  ofs << "#  move weeks (dummy number)" << endl;            // OS added
  ofs << "  1"  << endl;                     // OS added
  ofs << "#" << endl;           // OS added
  // ofs << "# Datasets  / LFIntervals  LFFirst  LFWidth  LFFactor" << endl;   // OS deleted
  // ofs << "   " << num_data_sets << "         " <<nlint<< "     " <<l1<<"    "  << delta <<  "  1 "  << endl;   // OS deleted
  ofs << "# Datasets/ LFIntervals LFFirst LFWidth LFFactor/ WFIntervals WFFirst WFWidth WFFactor (all WF flag was set to 0 in this simulation)" << endl;   // OS added
  ofs << "   " << num_data_sets << "         " <<nlint<< "       " <<l1<<"      "  << delta <<  "      1  0  0  0  0 "  << endl;                                          // OS added
  ofs << "# Year Month Week Fishery Catch Effort bin####" << endl;
  for (i=1;i<=nyrs;i++)
  {
    for (fi=1;fi<=num_fish;fi++)
    {
      double sumcatch=sum(value(C(fi,i)));
      if (sumcatch>cutoff)
      {
  //    ofs << 1959+i <<  "     2    1"   // OS temporaly modify
        ofs << 1959+i <<  "     " << catch_mth << "    " << catch_wek 
        << setfixed() << setprecision(2) << "    " << fi << "     " << obs_catch(fi,i) << "  " <<  obs_effort(fi,i) 
        << setfixed()  << setprecision(1) << "  " << obs_freq(fi,i) << endl;   // OS modified
        // ofs << "# Year   Month  Week" << endl;
        // ofs << "  " << 1959+i <<  "   1   1 " << endl;
        // ofs << "# Fishery Catch Effort" << endl;
        // ofs << setfixed()  << setprecision(2) << "   " << fi << "   " << obs_catch(fi,i) << "   " <<  obs_effort(fi,i) <<  endl;
        // ofs << setfixed()  << setprecision(1) << obs_freq(fi,i) << endl;
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
  ofstream ofs9("tpop.dat");    // OS add as the true age structure of population
  ofstream ofs10("tC.dat");     // OS add as the true catch number of fish
  ofstream ofs11("tCbio.dat");     // OS add as the true catch biomass of fish
  ofstream ofs12("tmu.dat");     // OS add as the true length of fish
  ofstream ofs13("twe.dat");     // OS add as the true weight of fish
  ofstream ofs14("tsel.dat");     // OS add as the true catch number of fish
  for (i=1;i<=nyrs;i++)
  {
    ofs1 << i<< " "<< setscientific()<< setprecision(5) << N(i,1)   << endl;
    ofs2 << i<< " "<< setscientific()<< setprecision(5) << "SSB(t): " << SSB(i) << "   Biomass(kg):  " << rbio(i)  << endl;
    ofs3 << i<< " "<< setscientific()<< setprecision(5) << relbio(i)<< endl;
    ofs4 << i<< " "<< setscientific()<< setprecision(5) << sum(N(i))<< endl;
  }
  ofs5 << setscientific()<< setprecision(5) << N(1)     << endl;
  ofs6 << setscientific()<< setprecision(5) << N(nyrs)  << endl;
  ofs7 << setscientific()<< setprecision(5) << q(1)*effort(1,1)*sel(1)  <<" "<<  q(2)*effort(2,1)*sel(2)  << endl;
  ofs8 << setscientific()<< setprecision(5) << q(1)*effort(1,nyrs)*sel(1)  <<" "<<  q(2)*effort(2,nyrs)*sel(2)  << endl;
  
  // OS add as the true age structure of population
  for (i=1;i<=nyrs;i++)
  {
    ofs9 << i << " "<<  setfixed()  << setprecision(2) << N(i) << endl;   // OS modified
    for (fi=1;fi<=num_fish;fi++)
    {
        ofs10 << 1959+i
        << setfixed() << setprecision(2) << "   " << fi << "    " << sum(value(C(fi,i))) << "   " <<  C(fi,i) << endl;   // OS modified
        ofs11 << 1959+i
        << setfixed() << setprecision(2) << "   " << fi << "    " << sum(value(Cbio(fi,i))) << "   " <<  Cbio(fi,i) << endl;   // OS modified
    }
    ofs12 << 1959+i
    << setfixed() << setprecision(2) << "   " <<  mu(i) << endl;   // OS modified
    ofs13 << 1959+i
    << setfixed() << setprecision(2) << "   " <<  pow(mu(i),lwcoff(2))*lwcoff(1) << endl; // OS modified
  }
  ofs14 << "Para:  l1: " << l1 << "  ln: "  << ln << "  delta: "  << delta << endl;
  ofs14 << "J:   " << J << endl;
  ofs14 << "mature:  " << mature << endl;
  ofs14 << "Age:   " << age << endl;
  for (fi=1;fi<=num_fish;fi++)
  {
    ofs14 << "Fleet:" << fi << "   " << sel(fi) << endl;
    // ofs11 << "Seldlt:" << fi << "   " << selplt(fi) << endl;
    // ofs11 << "Seldl:" << fi << "   " << selpl(fi) << endl;
  }

FUNCTION write_simulated_data_for_SS3
  int i,j,fi,fj;
  // ignoring spatial structure ,  OS add "j"   
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
  ofstream ofs("tunasim_out_forSS3.dat");
  ofs << "#V3.24f" << endl;
  ofs << "#Simulated data for SS3(ver3.24) analysis" << endl;
  ofs << "#" << endl;
  ofs << "###### Model Dimensions ##################" << endl;
  ofs << 1959+1 << "     # start year" << endl;
  ofs << 1959+nyrs << "     # end year" << endl;
  ofs << "1" << "     # number of seasons par year" << endl;
  ofs << " 12" << "     # number of monthes in each seasons" << endl;
  ofs << "1" << "     # spowning season" << endl;
  ofs << num_fish << "     # number of fishing fleet" << endl;
  ofs << num_fish << "     # number of surveys (cpue etc.)" << endl;
  ofs << nregions << "     # number of areas" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    ofs << "FISHERY" << fi << "%";
    }
  for (fi=1;fi<=num_fish;fi++) {
    if(fi < num_fish) {
       ofs << "SURVEY" << fi << "%";
       }
    else {
       ofs << "SURVEY" << fi << endl;
       }
    }
  for (fi=1;fi<=num_fish;fi++) {
    ofs << "0.50  ";
    }
  for (fi=1;fi<=num_fish;fi++) {
    if(fi < num_fish) {
       ofs << "0.50  ";
       }
    else{
       ofs << "0.50  " << "   # surveytiming_in_season" << endl;
       }
    }
  for (fi=1;fi<=num_fish;fi++) {
    ofs << "1  ";
    }
  for (fi=1;fi<=num_fish;fi++) {
    if(fi < num_fish) {
       ofs << "1  ";
       }
    else {
       ofs << "1  " << "   # area_assignments_for_each_fishery_and_survey" << endl;
       }
    }
  for (fi=1;fi<=num_fish;fi++) {
    if(fi < num_fish) {
       ofs << "1  ";
       }
    else {
       ofs << "1" << "     # units of catch for each fleet:  1=bio; 2=num" << endl;
       }
    }
  for (fi=1;fi<=num_fish;fi++) {
    if(fi < num_fish) {
       ofs << "0.05  ";
       }
    else {
       ofs << "0.05" << "     #_se of log(catch) for each fleet only used for init_eq_catch and for Fmethod 2 and 3; use -1 for discard only fleets" << endl;
       }
    }
  ofs << "1" << "     # number of genders" << endl;
  ofs << nages << "     # number of ages" << endl;
  //
  ofs << "#" << endl;
  ofs << "###### Catch ##################" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    if(fi < num_fish) {
       ofs << "0  ";
       }
    else {
       ofs << "0" << "     # init_equil_catch_for_each_fishery" << endl;
       }
    }
  ofs << nyrs << "     # number of lines for catch to read" << endl;
  ofs << "# catch by  each fleet by year and seasons" << endl;
  for (i=1;i<=nyrs;i++) {
    for (fi=1;fi<=num_fish;fi++) {
          ofs << obs_catbio(fi,i)/1000 << "   ";  // catch biomass(t) by each fisheries
          }
    ofs << 1959+i << "   1" << endl;  // year and season
    }
  //
  ofs << "#" << endl;
  ofs << "###### Abundance Index ##################" << endl;
  ofs << num_data_sets << "     # N_cpue_and_surveyabundance_observations" << endl;
  ofs << "# Units:  0=numbers; 1=biomass; 2=F" << endl;
  ofs << "# Errtype:  -1=normal; 0=lognormal; >0=T" << endl;
  ofs << "# fleet units errtype" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    ofs << fi << "   0   0     #FISHERY" << fi << endl;;
    }
  for (fi=1;fi<=num_fish;fi++) {
    ofs << fi+num_fish << "   0   0     #CPUE" << fi << endl;;
    }
  ofs << "# year season index obs_err" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    for (i=1;i<=nyrs;i++) {
      double sumcatch2=sum(value(C(fi,i)));
      if (sumcatch2>cutoff) {
         ofs << 1959+i <<  "   1   " << fi+num_fish << "   " << obs_catch(fi,i)/obs_effort(fi,i)/100000 << "   "<< 0.04*pow(10,fi-1) << endl;
         }
      }
    }
  //
  ofs << "#" << endl;
  ofs << "###### Discard ##################" << endl;
  ofs << "0     # number of fleets with discard" << endl;
  ofs << "0     # number of discard observation" << endl;
  //
  ofs << "#" << endl;
  ofs << "###### Mean-body weight ##################" << endl;
  ofs << "0     # number of mean bodyweight observed" << endl;
  ofs << "30     # DF for mean bodyweight T-distribution like" << endl;
  //
  ofs << "#" << endl;
  ofs << "###### Population length Bins ##################" << endl;
  //ofs << "1     # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector" << endl;
  ofs << "2     # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector" << endl; 
  ofs << delta << "     # binwidth for population size comp" << endl;
  ofs << l1 << "     # minimum size in the population (lower edge of first bin and size at age 0.00)" << endl;
  //ofs << l1-(delta/2) << "     # minimum size in the population (lower edge of first bin and size at age 0.00)" << endl;
  ofs << l1+(delta*nlint)-(delta/2) << "     # maximum size in the population (lower edge of last bin)" << endl;
  ofs << "# length composition #" << endl;
  ofs << "0.0001     # comp_tail_compression" << endl;
  ofs << "0.0001     # add_to_comp" << endl;
  ofs << "0      #_combine males into females at or below this bin number" << endl;
  ofs << nlint << "     # N_LengthBins" << endl;
  for (fj=1;fj<=nlint;fj++) {
    if(fj < nlint){
      ofs << " " << l1+(delta*(fj-1));
      }
    else {
      ofs << " " << l1+(delta*(fj-1)) << endl;
      }
    }
  ofs << num_data_sets*2 << "     # N_Length_obs" << endl;
  ofs << "#Yr  Season  Flt/Svy  Gender  Part  Nsamp  datavector(female+male)" << endl;
  for (fi=1;fi<=num_fish;fi++)
    {
    for (i=1;i<=nyrs;i++)
      {
      double sumcatch=sum(value(C(fi,i)));
      if (sumcatch>cutoff)
        {
        ofs << 1959+i <<  "     1     " << fi << "     0     0     "
        << setfixed()  << setprecision(1) << sum(obs_freq(fi,i)) << "     " << obs_freq(fi,i) << endl;
        }
      }
    }
  for (fi=1;fi<=num_fish;fi++)
    {
    for (i=1;i<=nyrs;i++)
      {
      double sumcatch=sum(value(C(fi,i)));
      if (sumcatch>cutoff)
        {
        ofs << 1959+i <<  "     1     " << num_fish+fi << "     0     0     "
        << setfixed()  << setprecision(1) << sum(obs_freq(fi,i)) << "     " << obs_freq(fi,i) << endl;
        }
      }
    }
  //
  ofs << "#" << endl;
  ofs << "###### Age Composition  ##################" << endl;
  ofs << "0     # N_age_bins" << endl;
  ofs << "0     # N_ageerror_definitions" << endl;
  ofs << "0     # N_Agecomp_obs" << endl;
  ofs << "2     # Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths" << endl; //?? 2?
  ofs << "0     # combine males into females at or below this bin number" << endl;
  ofs << "#Yr Seas Flt/Svy Gender Part Ageerr Lbin_lo Lbin_hi Nsamp datavector(female+male)" << endl;
  //
  ofs << "#" << endl;
  ofs << "###### Mean lenght or bodyweight at Age  ##################" << endl;
  ofs << "0     # N_MeanSize-at-Age_obs" << endl;
  //
  ofs << "#" << endl;
  ofs << "###### Environmental data  ##################" << endl;
  ofs << "0     # N_environ_variables" << endl;
  ofs << "0     # N_environ_obs"  << endl;
  //
  ofs << "#" << endl;
  ofs << "###### Generalized size composition data  ##################" << endl;
  ofs << "0     # N sizefreq methods to read" << endl; 
  //
  ofs << "#" << endl;
  ofs << "###### Tag recapture data  ##################" << endl;
  ofs << "0     # no tag data" << endl; 
  //
  ofs << "#" << endl;
  ofs << "###### Stock composition data  ##################" << endl;
  ofs << "0     # do morphcomp data" << endl; 
  ofs << "#" << endl;
  ofs << "999" << endl;
  //
  //
  ofstream ofs1("tunasim_out_forSS3.ctl");
  ofs1 << "#V3.24f" << endl;
  ofs1 << "#Simulated data for SS3(ver3.24) analysis" << endl;
  // 
  ofs1 << "#" << endl;
  ofs1 << "###### TERMINOLOGY ##################" << endl;
  ofs1 << "1    # number of growth patterns" << endl;
  ofs1 << "1    # number of submorphs within growth pattern" << endl;
  ofs1 << "0    # number of block patterns" << endl;
  // 
  ofs1 << "#" << endl;
  ofs1 << "###### BIOLOGY ##################" << endl;
  ofs1 << "0.5  # fraction female" << endl; 
  ofs1 << "#" << endl;
  ofs1 << "# M #" << endl;
  ofs1 << "0    #_natural mortality (naM) options (0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "# LAA #" << endl;
  ofs1 << "1    #_growth model (1=vonBert w/ L1&L2; 2=Richards w/ L1&L2; 3=age_speciific_K; 4=not_implemented)" << endl;
  ofs1 << " 0    #_growth age-min (1st size at age parameter)" << endl;
  ofs1 << " " << nages-1 << "     # growth age-max (2nd size at age parameter; 999 to use as Linf)" << endl;
  ofs1 << "0    #_sd add to length at age (LAA). Use 0.1 for SS2 ver1.x compatibility. 0 was Recommended." << endl;
  ofs1 << "0    #_cv growth pattern (0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A))" << endl;
  ofs1 << "#" << endl;
  ofs1 << "# Maturity #" << endl;
  ofs1 << "2    #_maturity_option (1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=read fec and wt from wtatage.ss)" << endl;
  ofs1 << "3    # first_mature_age" << endl;
  ofs1 << "1    # fecundity option ((1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W)" << endl;
  ofs1 << "0    # hermaphroditism option (0=none; 1=age-specific fxn)" << endl;
  ofs1 << "1    #_parameter offset method (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)" << endl;
  ofs1 << "2    #_time-varying adjustment constraint (env/block/dev_adj) (1=standard; 2=logistic transform keeps in base parm bounds; 3=standard w/ no bound check)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "#_Parameters #" << endl;
  ofs1 << "#_LOW    HIGH   INIT    PRIOR    PR    SD   PHASE   env   use   dev    dev    dev     Block  Block" << endl;
  ofs1 << "# ---    ----   ----    -----    -type --   -----   -var  -dev  -minyr -maxyr -stddev -----  -Fxn" << endl;
  ofs1 << "  0.01   0.50   0.20    0.20     -1    0.05   3    0     0     0      0      0       0      0   # [MGparam[1]] natM (PRIOR and SD is ignored (PR_type=-1))" << endl; 
  ofs1 << "#" << endl;
  ofs1 << "  0.50   50.0   20.0    25.0     0     10     2     0     0     0      0      0       0      0   # [MGparam[2]] body length(cm) of LAA(min)" << endl;
  ofs1 << "  90.0   200.0  170.0   180.0    0     10     4     0     0     0      0      0       0      0   # [MGparam[3]] body length(cm) of LAA(max)" << endl;
  ofs1 << "  0.05   0.50   0.17    0.18     0     0.8    4     0     0     0      0      0       0      0   # [MGparam[4]] growth coefficient of VonBert(K)" << endl;
  ofs1 << "  0.005  0.50   0.10    0.10     -1    0.8    -3    0     0     0      0      0       0      0   # [MGparam[5]] cv of LAA (young)" << endl; 
  ofs1 << "  0.005  0.25   0.10    0.10     -1    0.8    -3    0     0     0      0      0       0      0   # [MGparam[6]] cv of LAA (old)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "  -3.0   3.00   0.00002 0.00002  -1    0.8    -3    0     0     0      0      0       0      0   # [MGparam[7]] coefficient in L(cm)-W(kg)" << endl;
  ofs1 << "  -4.0   4.00   3.00    3.00     -1    0.8    -3    0     0     0      0      0       0      0   # [MGparam[8]] exponent in L(cm)-W(kg)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "  0.00   " << nages-1 << "     " << mage(1) << "       " << mage(1) << "        -1    0.8    3     0     0     0      0      0       0      0   # [MGparam[9]] maturity logistic inflection(cm or age) (Mat50%)" << endl;
  ofs1 << "  -50    -0.01  -0.50   -0.50    -1    0.8    3     0     0     0      0      0       0      0   # [MGparam[10]] maturity logistic slope (must have a negative value)" << endl;
  ofs1 << "  -3.0   3.00   1.00    1.00     -1    0.8    -3    0     0     0      0      0       0      0   # [MGparam[11]] fecundity parameter (Eggs/kg_inter)" << endl;
  ofs1 << "  -3.0   3.00   0.00    0.00     -1    0.8    -3    0     0     0      0      0       0      0   # [MGparam[12]] fecundity parameter (Eggs/kg_slope)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "  0.00   0.00   0.00    0.00     -1    0      -4    0     0     0      0      0       0      0   # [MGparam[13]] hermaphroditism parameter. Dont omit." << endl;
  ofs1 << "  0.00   0.00   0.00    0.00     -1    0      -4    0     0     0      0      0       0      0   # [MGparam[14]] hermaphroditism parameter. Dont omit." << endl;
  ofs1 << "  0.00   0.00   0.00    0.00     -1    0      -4    0     0     0      0      0       0      0   # [MGparam[15]] hermaphroditism parameter. Dont omit." << endl;
  ofs1 << "  0.00   0.00   0.00    0.00     -1    0      -4    0     0     0      0      0       0      0   # [MGparam[16]] cohort growth  dev. Dont omit." << endl;
  ofs1 << "#" << endl;
  ofs1 << "# seasonality for selected biology parmeters #" << endl;
  ofs1 << "  0   0   0   0   0   0   0   0   0   0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K (10 integers, Cannot omit.)" << endl;
  //
  ofs1 << "#" << endl;
  ofs1 << "## Spawner-Recruitment ##########################" << endl;
  ofs1 << "3   #_spawner recruitment specification functions (2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "#_Parameters #" << endl;
  ofs1 << "#_LOW   HI   INIT   PRIOR  PR_type  SD     PHASE" << endl;
  ofs1 << "  3.0   15   14     14     1        10.0   1    # [SR_parm[1]] SR_Log(R0)" << endl; 
  ofs1 << "  0.2   1    0.75   0.75   1        0.10   4   # [SR_parm[2]] SR_BH_steepness" << endl; 
  ofs1 << "  0.0   2    0.25   0.25   1        0.80   -4   # [SR_parm[3]] SR_sigmaR" << endl;
  ofs1 << "  -5.0  5    0.00   0.10   -1       1.00   -3   # [SR_parm[4]] SR_environmental linkage coefficient" << endl;
  ofs1 << "  -5.0  5    0.00   0.00   -1       1.00   -4   # [SR_parm[5]] SR_virgin recruitment(R1) offset" << endl;
  ofs1 << "  0.0   0    0.00   0.00   -1       0.00   -99  # [SR_parm[6]] SR_autocorrelation" << endl;
  ofs1 << "#" << endl;
  ofs1 << "# Environmental effect #" << endl;
  ofs1 << "0   #_environmental variable for adjustment of SR (SR_env_link)" << endl;
  ofs1 << "0   #_the factor of environmental variable (SR_env_target: 0=none;1=devs;_2=R0;_3=steepness)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "# Others #" << endl;
  ofs1 << "1   #_way for recruitment deviation (do_recdev: 0=none; 1=devvector; 2=simple deviations)" << endl;
  ofs1 << 1959+1 << " # first year of main recr_devs; early devs can preceed this era" << endl;
  ofs1 << 1959+ nyrs << " # last year of main recr_devs; forecast devs start in following year" << endl;
  ofs1 << "2   #_recdev phase" << endl;
  ofs1 << "1   # (0/1) to read 13 advanced options" << endl;
  ofs1 << " 0  #_recdev_early_start (0=none; neg value makes relative to recdev_start)" << endl;
  ofs1 << " -4 #_recdev_early_phase" << endl;
  ofs1 << " 0  #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)" << endl;
  ofs1 << " 1  #_lambda for Fcast_recr_like occurring before endyr+1" << endl;
  ofs1 << " 1900 #_last_early_yr_nobias_adj_in_MPD" << endl;
  ofs1 << " 1900 #_first_yr_fullbias_adj_in_MPD" << endl;
  ofs1 << " " << 1959+ nyrs << " #_last_yr_fullbias_adj_in_MPD" << endl;
  ofs1 << " " << 1959+ nyrs+ 1 << " #_first_recent_yr_nobias_adj_in_MPD" << endl;
  ofs1 << " 1  #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)" << endl;
  ofs1 << " 0  #_period of cycles in recruitment (N parms read below)" << endl;
  ofs1 << " -5 #min rec_dev" << endl;
  ofs1 << " 5  #max rec_dev" << endl;
  ofs1 << " 0  #_read_recdevs" << endl;
  ofs1 << "#_end of advanced SR options" << endl;
  //
  ofs1 << "#" << endl;
  ofs1 << "## Fishing Mortality Info ##########################" << endl;
  ofs1 << "0.3    #_F ballpark for tuning early phases" << endl;
  ofs1 << (-1)*(1959+nyrs) << " #_negative value disable F ballpark (F ballpark year)" << endl;
  ofs1 << "3      #_F_Method (1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended))" << endl;
  ofs1 << "4      #_maximum F or harvest rate, depends on F_Method (example is 2.9)" << endl;
  ofs1 << " 4      #_number of tuning itreations in hybrid method (recommend 3 to 7)" << endl;
  ofs1 << "#" << endl;
  ofs1 << "# Initial Fishing Mortality #" << endl;
  ofs1 << "#_LOW   HI   INIT   PRIOR  PR_type  SD     PHASE" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 << " 0.0  1    0.00   0.01   0        99     -1    # [init_F[" << fi << "]] FISHERY" << fi << endl;
    }
  ofs1 << "#" << endl;
  ofs1 << "# Catchability #" << endl;
  ofs1 << "#_Q_setup (Power, Environmnetal link, Extra sd, Types)" << endl;
  ofs1 << "#__Q_type options:  <0=mirror, 0=float_nobiasadj, 1=float_biasadj, 2=parm_nobiasadj, 3=parm_w_random_dev, 4=parm_w_randwalk, 5=mean_unbiased_float_assign_to_parm" << endl;
  ofs1 << "#__for_env-var:_enter_index_of_the_env-var_to_be_linked" << endl;
  ofs1 << "#_Den-dep  env-var  extra_se  Q_type" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 << " 0        0        0         0    # FISHERY" << fi << endl;
    }
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 << " 0        0        0         0    # CPUE" << fi << endl;
    }
  ofs1 << "# LOW   HI   INIT   PRIOR   PR_type   SD    PHASE" << endl;
  ofs1 << "#" << endl;
  ofs1 << "# Size Selectivity and Discard #" << endl;
  ofs1 << "#_discard_options (0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead)" << endl;
  ofs1 << "#_Pattern Discard Male Special" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 <<  sel_sw(fi) << "     0     0     0    # FISHERY" << fi << endl;
    }
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 <<  sel_sw(fi) << "     0     0     0    # CPUE" << fi << endl;
    }
  ofs1 << "#" << endl;
  ofs1 << "#_age_selex_types" << endl;
  ofs1 << "#_Pattern ___ Male Special" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 << " 10    0     0     0    # FISHERY" << fi << endl;
    }
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 << " 10    0     0     0    # CPUE" << fi << endl;
    }
  ofs1 << "#" << endl;
  ofs1 << "#_Parameters (length Sel)#" << endl;
  ofs1 << "#_LOW    HIGH   INIT    PRIOR    PR    SD   PHASE   env   use   dev    dev    dev     Block  Block" << endl;
  ofs1 << "# ---    ----   ----    -----    -type --   -----   -var  -dev  -minyr -maxyr -stddev -----  -Fxn" << endl;
  ofs1 << "# Fisheries #" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    //           ## LO ##                ## HI ##            ##  INIT  ##                              ##  PRIOR  ##                    # PR #  #SD#  PHASE   env  use  dev   dev   dev    Block  Block#
    //           ##    ##                ##    ##            ##        ##                              ##         ##                    #type#  #  #          var  dev  minyr maxyr stddev        Fxn #
    ofs1 << "  " << l1-delta << "    " << ln-delta << "  " << flselpeak1(fi)                << "    " << flselpeak1(fi)                << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " peak" <<  endl;
    if(sel_sw(fi)==1) {
      ofs1 << "  " << "-100"   << "  "   << "100"  << "  " << selpl(fi,2)*5                << "     " << selpl(fi,2)*5                 << "    0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " wideth" <<  endl;
    }
    else if(sel_sw(fi)==8) {
      ofs1 << "  " << "0   "   << "  "   << "0.8"    << "  " << selpl(fi,2)                   << "    " << selpl(fi,2)                   << "   0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " initMin" <<  endl;
      ofs1 << "  " << "-20 "   << "  "   << "10 "    << "  " << selpl(fi,3)                   << "    " << selpl(fi,3)                   << "   0   1000  -2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " infl1L" <<  endl;
      ofs1 << "  " << "-20 "   << "  "   << "20 "    << "  " << selpl(fi,4)                   << "    " << selpl(fi,4)                   << "   0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " slope1L" <<  endl;
      ofs1 << "  " << "-10 "   << "  "   << "10 "    << "  " << selpl(fi,5)                   << "    " << selpl(fi,5)                   << "   0   1000  -4      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " finalMax" <<  endl;
      ofs1 << "  " << "-20 "   << "  "   << "10 "    << "  " << selpl(fi,6)                   << "    " << selpl(fi,6)                   << "   0   1000  -3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " infl2R" <<  endl;
      ofs1 << "  " << "0   "   << "  "   << "15 "    << "  " << selpl(fi,7)                   << "    " << selpl(fi,7)                   << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " slope2R" <<  endl;
      ofs1 << "  " << "-233"   << "  "   << "233"    << "  " << flselpeak2(fi)-flselpeak1(fi) << "    " << flselpeak2(fi)-flselpeak1(fi) << "   0   1000   4      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " wideth" <<  endl;
    }
    else if(sel_sw(fi)==22) {
      ofs1 << "  " << "-10 "   << "  "   << "10 "    << "  " << selpl(fi,2)                   << "  " << selpl(fi,2)                   << "   0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " peak2" <<  endl;
      ofs1 << "  " << "0   "   << "  "   << "10 "    << "  " << selpl(fi,3)                   << "  " << selpl(fi,3)                   << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " upslope" <<  endl;
      ofs1 << "  " << "0   "   << "  "   << "10 "    << "  " << selpl(fi,4)                   << "  " << selpl(fi,4)                   << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " downslope" <<  endl;
    }
    ofs1 << "#" << endl;
  }
  ofs1 << "# CPUE #" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    //           ## LO ##                ## HI ##            ##  INIT  ##                              ##  PRIOR  ##                    # PR #  #SD#  PHASE   env  use  dev   dev   dev    Block  Block#
    //           ##    ##                ##    ##            ##        ##                              ##         ##                    #type#  #  #          var  dev  minyr maxyr stddev        Fxn #
    ofs1 << "  " << l1-delta << "    " << ln-delta << "  " << flselpeak1(fi)                << "    " << flselpeak1(fi)                << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " peak" <<  endl;
    if(sel_sw(fi)==1) {
      ofs1 << "  " << "-100"   << "  "   << "100"  << "  " << selpl(fi,2)*5                << "     " << selpl(fi,2)*5                 << "    0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " wideth" <<  endl;
    }
    else if(sel_sw(fi)==8) {
      ofs1 << "  " << "0   "   << "  "   << "0.8"    << "  " << selpl(fi,2)                   << "    " << selpl(fi,2)                   << "   0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " initMin" <<  endl;
      ofs1 << "  " << "-20 "   << "  "   << "10 "    << "  " << selpl(fi,3)                   << "    " << selpl(fi,3)                   << "   0   1000  -2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " infl1L" <<  endl;
      ofs1 << "  " << "-20 "   << "  "   << "20 "    << "  " << selpl(fi,4)                   << "    " << selpl(fi,4)                   << "   0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " slope1L" <<  endl;
      ofs1 << "  " << "-10 "   << "  "   << "10 "    << "  " << selpl(fi,5)                   << "    " << selpl(fi,5)                   << "   0   1000  -4      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " finalMax" <<  endl;
      ofs1 << "  " << "-20 "   << "  "   << "10 "    << "  " << selpl(fi,6)                   << "    " << selpl(fi,6)                   << "   0   1000  -3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " infl2R" <<  endl;
      ofs1 << "  " << "0   "   << "  "   << "15 "    << "  " << selpl(fi,7)                   << "    " << selpl(fi,7)                   << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " slope2R" <<  endl;
      ofs1 << "  " << "-233"   << "  "   << "233"    << "  " << flselpeak2(fi)-flselpeak1(fi) << "    " << flselpeak2(fi)-flselpeak1(fi) << "   0   1000   4      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " wideth" <<  endl;
    }
    else if(sel_sw(fi)==22) {
      ofs1 << "  " << "-10 "   << "  "   << "10 "    << "  " << selpl(fi,2)                   << "  " << selpl(fi,2)                   << "   0   1000   2      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " peak2" <<  endl;
      ofs1 << "  " << "0   "   << "  "   << "10 "    << "  " << selpl(fi,3)                   << "  " << selpl(fi,3)                   << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " upslope" <<  endl;
      ofs1 << "  " << "0   "   << "  "   << "10 "    << "  " << selpl(fi,4)                   << "  " << selpl(fi,4)                   << "   0   1000   3      0    0    0     0     0.5    0      0    # [selparm[*]]#" << fi << " downslope" <<  endl;
    }
    ofs1 << "#" << endl;
  }
  ofs1 << "#" << endl;
  ofs1 << "#_Parameters (Age Sel)#" << endl;
  ofs1 << "# Fisheries #" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    //            ## LO ##         ## HI ##        ##  INIT  ##       ##  PRIOR  ##    PR    SD    PHASE  env   use  dev    dev   dev    Block  Block
    //            ##    ##         ##    ##        ##        ##       ##         ##    type               var   dev  minyr  maxyr stddev        Fxn 
    ofs1 << "# " << "0.00" << " " << nages << "  " << "0.00" << "  "  << "0.00" << "   -1    99.9   -1    0     0    0      0     0      0      0    # AgeSel_1P_1_#" << fi << endl;
    ofs1 << "# " << "0.00" << " " << nages << "  " << nages << "    " << nages  << "   -1    99.9   -1    0     0    0      0     0      0      0    # AgeSel_1P_2_#" << fi << endl;
    ofs1 << "#" << endl;
    }
    ofs1 << "# CPUE #" << endl;
  for (fi=1;fi<=num_fish;fi++) {
    ofs1 << "# " << "0.00" << " " << nages << "  " << "0.00" << "  "  << "0.00" << "   -1    99.9   -1    0     0    0      0     0      0      0    # AgeSel_1P_1_#" << fi << endl;
    ofs1 << "# " << "0.00" << " " << nages << "  " << nages << "    " << nages  << "   -1    99.9   -1    0     0    0      0     0      0      0    # AgeSel_1P_2_#" << fi << endl;
    ofs1 << "#" << endl;
    } 
  ofs1 << "#" << endl;
  ofs1 << "# Tag loss and Tag reporting parameters  ##########################" << endl;
  ofs1 << "0   # TG_custom:  0=no read; 1=read if tags exist" << endl;
  //
  ofs1 << "#" << endl;
  ofs1 << "# Variance Adjustment Factors ##########################" << endl;
  ofs1 << "0 #_Variance_adjustments_to_input_values" << endl;
  ofs1 << "#_fleet and cpue: 1 2 3 4" << endl;
  ofs1 << "# 0    0    0    0   #_add_to_survey_CV" << endl;
  ofs1 << "# 0    0    0    0   #_add_to_discard_stddev" << endl;
  ofs1 << "# 0    0    0    0   #_add_to_bodywt_CV" << endl;
  ofs1 << "# 1    1    1    1   #_mult_by_lencomp_N" << endl;
  ofs1 << "# 1    1    1    1   #_mult_by_agecomp_N" << endl;
  ofs1 << "# 1    1    1    1   #_mult_by_size-at-age_N" << endl; 
  //
  ofs1 << "#" << endl;
  ofs1 << "# Lambdas ##########################" << endl;
  ofs1 << "4   #_maxlambdaphase" << endl;
  ofs1 << "2   #_sd_offset" << endl;
  ofs1 << "#" << endl;
  ofs1 << "2   # number of changes to make to default Lambdas (default value is 1.0)" << endl;
  ofs1 << "# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch;" << endl; 
  ofs1 << "# 9=init_equ_catch; 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin" << endl;
  ofs1 << "#like_comp fleet/survey  phase  value  sizefreq_method" << endl;
  ofs1 << " 1   3   3   1   1" << endl;
  ofs1 << " 1   4   3   1   1" << endl;
  ofs1 << " # 4   1   3   1   1" << endl;
  ofs1 << " # 4   2   3   1   1" << endl;
  ofs1 << " # 4   3   3   1   1" << endl;
  ofs1 << " # 4   4   3   1   1" << endl;
  //
  ofs1 << "#" << endl;
  ofs1 << "# Controls for variance of derived quantities ##########################" << endl;
  ofs1 << "0 # (0/1) read specs for more stddev reporting" << endl; 
  // 
  ofs1 << "999" << endl;


TOP_OF_MAIN_SECTION
  arrmblsize=500000;

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
