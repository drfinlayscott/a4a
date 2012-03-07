DATA_SECTION
  init_int Syr;
  !!cout << Syr << endl;
  init_int Fyr;
  !!cout << Fyr << endl;
  init_vector Land(Syr,Fyr);
  !!cout << Land << endl;
  init_int NBiom;
  !!cout << NBiom << endl;
  init_matrix CatDynBiom(1,NBiom,1,3);
  !!cout << CatDynBiom << endl;
  vector YrBiom(1,NBiom);
  vector MLBiom(1,NBiom);
  vector SDBiom(1,NBiom);
INITIALIZATION_SECTION
  logB0  11.08;
  logK 11.23;
  logr 0.405;
  logp 0.693;
PARAMETER_SECTION
  init_bounded_number logB0(5.5,15.0,2);
  init_bounded_number logK(9.5,15.0,2);
  init_bounded_number logr(0.1,2.5,3);
  init_bounded_number logp(-1.0,2.0,1);
  vector Biom(Syr,Fyr+1);
  vector PredBiom(1,NBiom);
  vector SquDiff(1,NBiom);
  number Prod;
  sdreport_number B0;
  sdreport_number K;
  sdreport_number r;
  sdreport_number p;
  //sdreport_vector Biom_sd(Syr,Fyr+1);
  objective_function_value ff;
PRELIMINARY_CALCS_SECTION
  YrBiom=column(CatDynBiom,1);
  cout<<"YrBiom="<<YrBiom<<endl;
  MLBiom=column(CatDynBiom,2);
  cout<<"MLBiom="<<MLBiom<<endl;
  SDBiom=column(CatDynBiom,3);
  cout<<"SDBiom="<<SDBiom<<endl;
PROCEDURE_SECTION
  B0=mfexp(logB0);
  cout<<"B0="<<B0<<endl;
  K=mfexp(logK);
  cout<<"K="<<K<<endl;
  r=mfexp(logr);
  cout<<"r="<<r<<endl;
  p=mfexp(logp);
  cout<<"p="<<p<<endl;
  int yr;
  int yrBiom;
  Biom(Syr)=B0;
  for(yr=Syr+1;yr<=Fyr+1;yr++)
    {
     Prod=Biom(yr-1)+r*Biom(yr-1)*(1-pow(Biom(yr-1)/K,p));
     if(Prod<=Land(yr-1)) Prod=Land(yr-1)+1000;
     Biom(yr)=Prod-Land(yr-1);
    }
  cout<<"Biom="<<Biom<<endl;
  dvariable loglikBiom;
  loglikBiom=0; 
  yrBiom=1;
  for(yr=Syr;yr<=Fyr;yr++)
    {
      if(yr==YrBiom(yrBiom))
        {
         PredBiom(yrBiom)=(Biom(yr)+Biom(yr+1))/2.;
         SquDiff(yrBiom)=square(MLBiom(yrBiom)-PredBiom(yrBiom));
         loglikBiom = loglikBiom - 0.5*(log(2*3.1416*square(SDBiom(yrBiom)))+SquDiff(yrBiom)/square(SDBiom(yrBiom)));
         yrBiom=yrBiom+1;
        }
    }
  cout<<"yrBiom="<<yrBiom<<endl;
  cout<<"Year="<<yr<<endl;
  cout<<"PredBiom="<<PredBiom<<endl;
  cout<<"SquDiff="<<SquDiff<<endl;
  cout<<"loglikBiom"<<loglikBiom<<endl;
  ff = -(loglikBiom);
  //Biom_sd=Biom;
