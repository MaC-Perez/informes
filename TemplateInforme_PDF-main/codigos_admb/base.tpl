// Southern hake model

GLOBALS_SECTION
  #include <admodel.h>
  #include <time.h>
  #include <string.h>
  #undef depur
  #undef depuro
  #define depur(object) cout << #object "\n" << object << endl;
  #define depuro(object) cout << #object "\n" << object << endl; exit(1);


  #undef reporte
  #define reporte(object) report << #object "\n" << object << endl;

  #if defined(WIN32) && !defined(__linux__)
      const char* PLATFORM = "Windows";
  #else
      const char* PLATFORM = "Linux";
  #endif

  adstring BaseFileName;
  adstring ReportFileName;
  //adstring ResultsPath;

  adstring stripExtension(adstring fileName)
	{
		/* from Stock-Sintesis */
		const int length = fileName.size();
		for (int i=length; i>=0; --i)
		{
			if (fileName(i)=='.')
			{
				return fileName(1,i-1);
			}
		}
		return fileName;
	}


TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(7000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(7000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  arrmblsize=500000000;


DATA_SECTION
  !! cout<<"Modelo de Merluza del Sur corriendo bajo "<<PLATFORM<< endl;

  init_adstring DataFile;
  init_adstring ControlFile;
  //init_adstring ResultsFileName;

	!! BaseFileName = stripExtension(ControlFile);
	!! cout << "dat: " << " " << DataFile << endl;
	!! cout << "ctl: " << " " << ControlFile << endl;
	!! cout << "basefileName: " << " " << BaseFileName << endl;
	!! ReportFileName = BaseFileName + adstring(".rep");

  !! ad_comm::change_datafile_name(DataFile);
  init_int nyears;
  init_int nages;
  init_ivector vyears(1,nyears);
  int styr               //1977
  !! styr = vyears(1);
  int styr_pop           //1978
  !! styr_pop = vyears(2);
  int endyr
  !! endyr = vyears(nyears);
  init_ivector vages(1,nages);
  int stage
  !! stage = vages(1);
  int endage
  !! endage = vages(nages);
  init_vector surveyindex(styr,endyr);
  init_matrix cpueindex(styr,endyr,1,3);

  vector indxtrawl(styr,endyr);
  vector indxlongline(styr,endyr);
  vector indxartisanal(styr,endyr);
  !! indxtrawl = column(cpueindex,1);
  !! indxlongline = column(cpueindex,2);
  !! indxartisanal = column(cpueindex,3);

  init_matrix landing(styr,endyr,1,3);
  vector ytrawl(styr,endyr);
  vector ylongline(styr,endyr);
  vector yartisanal(styr,endyr);
  !! ytrawl = column(landing,1);
  !! ylongline = column(landing,2);
  !! yartisanal = column(landing,3);
  init_matrix catagetrwl(styr,endyr,stage,endage);
  init_matrix catagelongl(styr,endyr,stage,endage);
  init_matrix catageartisa(styr,endyr,stage,endage);
  init_matrix natagesurvey(styr,endyr,stage,endage);
  init_matrix watagefleet_pop(styr,endyr,stage,endage);
  init_matrix watagefleet_trawl(styr,endyr,stage,endage);
  init_matrix watagefleet_long(styr,endyr,stage,endage);
  init_matrix watagefleet_arti(styr,endyr,stage,endage);
  init_matrix watagefleet_surv(styr,endyr,stage,endage);
  matrix Wm_pop(styr,endyr,stage,endage);
  matrix Wm_trawl(styr,endyr,stage,endage);
  matrix Wm_long(styr,endyr,stage,endage);
  matrix Wm_arti(styr,endyr,stage,endage);
  matrix Wm_surv(styr,endyr,stage,endage);
  !! Wm_pop = (watagefleet_pop)/1e3;
  !! Wm_trawl = (watagefleet_trawl)/1e6;
  !! Wm_long = (watagefleet_long)/1e6;
  !! Wm_arti = (watagefleet_arti)/1e6;
  !! Wm_surv = (watagefleet_surv)/1e6;
  init_number M
  init_matrix msex(styr,endyr,stage,endage);
  init_number offset


  !! ad_comm::change_datafile_name(ControlFile);
  init_int phs_init
  init_int phs_R
  init_int phs_q
  init_int phs_Sel
  init_int phs_F
  init_number h
  init_int Ptrawl_1 //periodo con desembarque arrastre inicio
  init_int Ptrawl_2 //periodo con desembarque arrastre final
  int PeT
  !! PeT = Ptrawl_2 - Ptrawl_1 + 1;
  init_int Ppal_1 //periodo con desembarque palangre inicio
  init_int Ppal_2 //periodo con desembarque palangre final
  int PeP
  !! PeP = Ppal_2 - Ppal_1 + 1;
  init_int Pesp_1 //periodo con desembarque espinel inicio
  init_int Pesp_2 //periodo con desembarque espinel final
  int PeE
  !! PeE = Pesp_2 - Pesp_1 + 1;
  init_vector chQarr(1,3)
  init_vector chQpal(1,2)
  init_int inxtrawl
  init_vector strawl(1,inxtrawl)
  init_int inxpal
  init_vector spal(1,inxpal)
  init_int inxesp
  init_vector sesp(1,inxesp)
  init_int inxsurv
  init_vector ssurv(1,inxsurv)
  init_vector nss(1,4)
  init_int Pdarr_1
  init_int Pdarr_2
  init_int Pdpal_1
  init_int Pdpal_2
  init_matrix cv_matrix(styr,endyr,1,5)
  init_vector peso_index(1,4)

  vector cv_arr(styr,endyr)
  vector cv_pal(styr,endyr)
  vector cv_cru(styr,endyr)
  vector cv_art(styr,endyr)
  !! cv_arr = peso_index(1)*column(cv_matrix,2);
  !! cv_pal = peso_index(2)*column(cv_matrix,3);
  !! cv_cru = peso_index(3)*column(cv_matrix,4);
  !! cv_art = peso_index(4)*column(cv_matrix,5);

  init_vector cv_s(1,4)
  init_number cv_p
  init_vector rango_sa(1,4)
  init_vector rango_sl(1,4)
  init_number cv_sel_a
  init_number cv_sel_c

  init_int yr_sim
  init_int nFt
  init_vector mf(1,nFt)
  init_number recOp
  init_number rCap
  init_number Fpbr
  init_number offsetCt
  //!! depuro(offset)

INITIALIZATION_SECTION
  log_selA 2.7
  log_selB 1.5
  log_selC 200
  log_qpal -10
  log_qcru 0
  log_qarr -5
  log_Ro 18.5
  mu 0
  log_Farr -1.38
  log_Fesp -2.52
  log_Fpal -1.96

PARAMETER_SECTION
  init_bounded_vector log_selA(1,4,1.5,3.09,phs_Sel)
  init_bounded_vector log_selB(1,4,0.5,2.50,phs_Sel)
  init_bounded_vector log_selC(1,4,1,5.30,phs_Sel)

  init_vector log_qpal(1,3,phs_q)//(phs_q)
  init_number log_qcru(phs_q)
  init_number log_qesp(phs_q)//(1,2,phs_q)
  init_vector log_qarr(1,4,phs_q)//(phs_q)

  init_bounded_number log_Ro(17,20,phs_init)//10,20
  init_bounded_vector mu(styr_pop,endyr,-1,1,phs_R)
  init_bounded_vector log_Farr(Ptrawl_1,Ptrawl_2,-9,0.6,phs_F)
  init_bounded_vector log_Fesp(Pesp_1,Pesp_2,-9,0.6,phs_F)
  init_bounded_vector log_Fpal(Ppal_1,Ppal_2,-9,0.6,phs_F)


  likeprof_number Ro_pl
  sdreport_number Ro
  sdreport_number spr
  sdreport_number So
  sdreport_vector SB(styr,endyr)
  sdreport_vector R(styr,endyr)
  sdreport_vector BT(styr,endyr)
  sdreport_vector B6(styr,endyr)
  sdreport_vector S_pal(1,nages)
  sdreport_vector S_arr(1,nages)
  sdreport_vector S_esp(1,nages)
  sdreport_vector S_cru(1,nages)
  sdreport_vector Farr(Ptrawl_1,Ptrawl_2)
  sdreport_vector Fpal(Ppal_1,Ppal_2)
  sdreport_vector Fesp(Pesp_1,Pesp_2)
  sdreport_vector Ftot(styr,endyr)
  sdreport_vector muArr(styr,endyr)
  sdreport_vector muPal(styr,endyr)
  sdreport_vector muEsp(styr,endyr)
  sdreport_matrix BDp(endyr,endyr+yr_sim,1,nFt)
  sdreport_matrix Yproy(endyr,endyr+yr_sim,1,nFt)
  sdreport_matrix Fproy(endyr,endyr+yr_sim,1,nFt)
  sdreport_vector SBdeple(styr,endyr)
  sdreport_vector sprdeple(styr,endyr)
  sdreport_vector Bdepl(styr,endyr)

  number cuenta1
  number cuenta2
  number cuenta3
  number cuenta4

  // -----------------------------------------------
  number alpha
  number beta
  matrix No(styr,endyr,stage,endage)
  matrix NS(styr,endyr,stage,endage)
  vector NSo(stage,endage)
  vector uno_ages(1,nages)
  vector uno_years(styr,endyr)
  vector uno_years_arr(Ptrawl_1,Ptrawl_2)
  vector uno_years_pal(Ppal_1,Ppal_2)
  vector uno_years_esp(Pesp_1,Pesp_2)
  matrix Fcr_arr(Ptrawl_1,Ptrawl_2,stage,endage)
  matrix Fcr_pal(Ppal_1,Ppal_2,stage,endage)
  matrix Fcr_esp(Pesp_1,Pesp_2,stage,endage)
  matrix Fcr_total(styr,endyr,stage,endage)
  matrix Z(styr,endyr,stage,endage)
  matrix Surv(styr,endyr,stage,endage)
  matrix ND(styr,endyr,stage,endage)
  vector BDcru(styr,endyr)
  matrix Zarr(styr,endyr,stage,endage)
  matrix Zpal(styr,endyr,stage,endage)
  matrix Zesp(styr,endyr,stage,endage)
  vector BMVarr(styr,endyr)
  vector BMVpal(styr,endyr)
  vector BMVesp(styr,endyr)
  vector CPUEarr(styr,endyr)
  vector CPUEpal(styr,endyr)
  vector CPUEesp(styr,endyr)
  vector estBDcru(styr,endyr)
  matrix cageArr(styr,endyr,stage,endage)
  matrix cagePal(styr,endyr,stage,endage)
  matrix cageEsp(styr,endyr,stage,endage)
  matrix NSurvey(styr,endyr,stage,endage)
  vector YestArr(styr,endyr)
  vector YestPal(styr,endyr)
  vector YestEsp(styr,endyr)
  matrix pobsarr(1,inxtrawl,stage,endage)
  matrix pestarr(1,inxtrawl,stage,endage)
  matrix pobspal(1,inxpal,stage,endage)
  matrix pestpal(1,inxpal,stage,endage)
  matrix pobsesp(1,inxesp,stage,endage)
  matrix pestesp(1,inxesp,stage,endage)
  matrix pobssurv(1,inxsurv,stage,endage)
  matrix pestsurv(1,inxsurv,stage,endage)
  vector logL(1,12)
  vector penL(1,12)
  number a
  number sl
  number sr
  vector p(1,2)
  vector d_arr(1,2)
  vector d_pal(1,2)
  vector d_esp(1,2)
  vector d_cru(1,2)

  objective_function_value objF

  // proyecciones
  matrix Np(endyr,endyr+yr_sim,stage,endage)
  matrix Sp(endyr,endyr+yr_sim,stage,endage)
  matrix Fp(endyr,endyr+yr_sim,stage,endage)
  matrix Zp(endyr,endyr+yr_sim,stage,endage)
  matrix NSp(endyr,endyr+yr_sim,stage,endage)
  matrix Sep(endyr,endyr+yr_sim,stage,endage)
  matrix caep(endyr,endyr+yr_sim,stage,endage)

  vector Ftp(endyr,endyr+yr_sim)
  vector SBp(endyr,endyr+yr_sim)
  vector Rp(endyr,endyr+yr_sim)
  vector Yp(endyr,endyr+yr_sim)
  vector wp(stage,endage)
  vector msp(stage,endage)


PRELIMINARY_CALCS_SECTION
  uno_ages = 1;
  uno_years = 1;
  uno_years_arr = 1;
  uno_years_pal = 1;
  uno_years_esp = 1;
  Ro_pl.set_stepnumber(50);
  Ro_pl.set_stepsize(0.2);

RUNTIME_SECTION
  //convergence_criteria 1.e-1,1.e-01,1.e-03,1e-5,1e-5
  //maximum_function_evaluations 100,100,200,300,2500
  maximum_function_evaluations 10000, 30000, 100000, 500000
  convergence_criteria 1e-7,1e-8,1e-8, 1e-8

PROCEDURE_SECTION
  selectivity_exploitation_rate();
  initial_age_structure();
  selectivity_penalties();
  dynamics_abundance_per_fleet();
  biomass_and_mortality();
  estimates_cpue_fleet();
  catch_at_age();
  biomass_state();
  evaluate_objective_function();

  if(last_phase())
    {
    sim_Fcte();
    }

  if(mceval_phase())
    {
    ofstream out("model_msur_ctp2017.var.mcmc",ios::app);
    out << Ro << " " << objF << " " << Ftot << " " << SB << endl;
    out.close();
    }


FUNCTION selectivity_exploitation_rate
  int t,i;

     // Palangre
     a  = mfexp(log_selA(1));
     sl = mfexp(log_selB(1));
     sr = mfexp(log_selC(1));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_pal(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_pal(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_pal = S_pal / (max(S_pal)+1e-6);


     // Arrastre
     a  = mfexp(log_selA(2));
     sl = mfexp(log_selB(2));
     sr = mfexp(log_selC(2));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_arr(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_arr(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_arr = S_arr / (max(S_arr)+1e-6);


     // Espinel
     a  = mfexp(log_selA(3));
     sl = mfexp(log_selB(3));
     sr = mfexp(log_selC(3));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_esp(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_esp(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_esp = S_esp / (max(S_esp)+1e-6);


     // Cruceros
     a  = mfexp(log_selA(4));
     sl = mfexp(log_selB(4));
     sr = mfexp(log_selC(4));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_cru(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_cru(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_cru = S_cru / (max(S_cru)+1e-6);


  Fpal = mfexp(log_Fpal); //-1.96
  Fesp = mfexp(log_Fesp);
  Farr = mfexp(log_Farr);

  Fcr_pal = elem_prod(outer_prod(uno_years_pal,S_pal),outer_prod(Fpal,uno_ages));
  Fcr_esp = elem_prod(outer_prod(uno_years_esp,S_esp),outer_prod(Fesp,uno_ages));
  Fcr_arr = elem_prod(outer_prod(uno_years_arr,S_arr),outer_prod(Farr,uno_ages));

  for(t=styr; t<=endyr; t++)
    {
    if(t < Pesp_1)      // si t < 1981
      { Fcr_total(t) = Fcr_arr(t);}
    else if (t >= Pesp_1 & t < Ppal_1) // si t esta entre 1981 y 1987
      { Fcr_total(t) = Fcr_arr(t) + Fcr_esp(t);}
    else                                           // todo lo demas
      {Fcr_total(t) = Fcr_arr(t) + Fcr_esp(t) + Fcr_pal(t);}
    }

  Z = Fcr_total + M;  //Fcr_total matriz que va desde 1 nyears 1 a edades, mortalidad por pesca total a la edad
  Surv = mfexp(-1.0 * Z);

      for(t=styr; t<=endyr; t++)
    {
    if(t < Pesp_1)
      { Ftot(t) = Farr(t);}
    else if (t >= Pesp_1 & t < Ppal_1)
      { Ftot(t) = Farr(t) + Fesp(t);}
    else
      { Ftot(t) = Farr(t) + Fesp(t) + Fpal(t);} //moratlidad por pesca total anual, vector 1, nyears
    }

FUNCTION initial_age_structure
  int a;

  Ro = mfexp(log_Ro);
  Ro_pl = Ro;
  No(styr,stage) = Ro;
  for(a=stage+1; a<=endage; a++)
    {
    No(styr,a) = No(styr,a-1)*mfexp(-M);
       if (a==endage)
       {
       No(styr,a) += No(styr,a)/(1.0-mfexp(-M));
        }
    }
  NSo = elem_prod(elem_prod(extract_row(No,styr),msex(styr)),Wm_pop(styr))*mfexp(-M*9/12);
  So = sum(NSo);
  alpha = (So/Ro)*(1.0-h)/(4.0*h);
  beta = (5.0*h-1.0)/(4.0*h*Ro);
  spr = (So/Ro)*1e6;


FUNCTION selectivity_penalties
  int t;

     // Palangre
     a  = mfexp(log_selA(1));
     sl = mfexp(log_selB(1));
     sr = mfexp(log_selC(1));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_pal(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
           {d_pal(t) = pow(2,-1*square((p(t)-a)/sr));}
     }

     // Arrastre
     a  = mfexp(log_selA(2));
     sl = mfexp(log_selB(2));
     sr = mfexp(log_selC(2));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_arr(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
           {d_arr(t) = pow(2,-1*square((p(t)-a)/sr));}
     }

     // Espinel
     a  = mfexp(log_selA(3));
     sl = mfexp(log_selB(3));
     sr = mfexp(log_selC(3));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_esp(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
           {d_esp(t) = pow(2,-1*square((p(t)-a)/sr));}
     }

     // Cruceros
     a  = mfexp(log_selA(4));
     sl = mfexp(log_selB(4));
     sr = mfexp(log_selC(4));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_cru(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
         {d_cru(t) = pow(2,-1*square((p(t)-a)/sr));}
     }


FUNCTION dynamics_abundance_per_fleet
  int t;

  R(styr) = So/(alpha+(beta*So));
  NS.rowfill(styr,(elem_prod(elem_prod(elem_prod(extract_row(No,styr),msex(styr)),Wm_pop(styr)),mfexp(-1.0*Z(styr)*9/12))));
  SB(styr) = sum(extract_row(NS,styr));
  SBdeple(styr) = SB(styr)/So;
  sprdeple(styr) = (SB(styr)/R(styr))*1e6;

  for(t=styr_pop; t<=endyr; t++) {
    R(t) = SB(t-1)/(alpha+(beta*SB(t-1)))*mfexp(mu(t));
    No(t,1) = R(t);
    No(t)(stage+1,endage) =  ++elem_prod(No(t-1)(stage, endage - 1),Surv(t-1)(stage, endage - 1));
    No(t,endage) += No(t,endage)/(1-Surv(t-1,endage));
    NS(t) = elem_prod(elem_prod(elem_prod(No(t),msex(t)),Wm_pop(t)),mfexp(-1.0*Z(t)*9/12));
    SB(t) =  sum(extract_row(NS,t));
    SBdeple(t) = SB(t)/So;
    sprdeple(t) = (SB(t)/R(t))*1e6;
  }


FUNCTION biomass_and_mortality
  int t;
  ND = elem_prod(No,msex)*mfexp(-M*9/12);
  BDcru = rowsum(elem_prod(elem_prod(ND,outer_prod(uno_years,S_cru)),Wm_surv));

  Zarr = Fcr_arr + M;
  BMVarr = rowsum(elem_prod(elem_div(1-mfexp(-1.0*Zarr),Zarr),elem_prod(elem_prod(No,outer_prod(uno_years,S_arr)),Wm_trawl)));
  muArr = elem_div(ytrawl+1e-6,BMVarr);

  Zpal = M; for(t=Ppal_1; t<=Ppal_2; t++) {Zpal(t) += Fcr_pal(t);}
  BMVpal = rowsum(elem_prod(elem_div(1-mfexp(-1.0*Zpal),Zpal),elem_prod(elem_prod(No,outer_prod(uno_years,S_pal)),Wm_long)));
  muPal = elem_div(ylongline+1e-6,BMVpal);

  Zesp = M; for(t=Pesp_1; t<=Pesp_2; t++) {Zesp(t) += Fcr_esp(t);}
  BMVesp = rowsum(elem_prod(elem_div(1-exp(-1.0*Zesp),Zesp),elem_prod(elem_prod(No,outer_prod(uno_years,S_esp)),Wm_arti)));
  muEsp = elem_div(yartisanal+1e-6,BMVesp);


FUNCTION estimates_cpue_fleet
  CPUEarr(styr,chQarr(1))        = mfexp(log_qarr(1))*BMVarr(styr,chQarr(1));
  CPUEarr(chQarr(1)+1,chQarr(2)) = mfexp(log_qarr(2))*BMVarr(chQarr(1)+1,chQarr(2));
  CPUEarr(chQarr(2)+1,chQarr(3)) = mfexp(log_qarr(3))*BMVarr(chQarr(2)+1,chQarr(3));
  CPUEarr(chQarr(3)+1,endyr)     = mfexp(log_qarr(4))*BMVarr(chQarr(3)+1,endyr);

  CPUEpal(styr,chQpal(1))        = mfexp(log_qpal(1))*BMVpal(styr,chQpal(1));
  CPUEpal(chQpal(1)+1,chQpal(2)) = mfexp(log_qpal(2))*BMVpal(chQpal(1)+1,chQpal(2));
  CPUEpal(chQpal(2)+1,endyr)     = mfexp(log_qpal(3))*BMVpal(chQpal(2)+1,endyr);

  CPUEesp  = mfexp(log_qesp)*BMVesp;
  estBDcru = mfexp(log_qcru)*BDcru;

FUNCTION catch_at_age
  double tiny = 1e-6;
  int t;

  cageArr = elem_prod(elem_prod(No,Zarr - M),elem_div(1-mfexp(-1.0*Z),Z));
          YestArr = rowsum(elem_prod(Wm_trawl,cageArr));

  cagePal = elem_prod(elem_prod(No,Zpal - M),elem_div(1-mfexp(-1.0*Z),Z));
          YestPal = rowsum(elem_prod(Wm_long,cagePal));

  cageEsp = elem_prod(elem_prod(No,Zesp - M),elem_div(1-mfexp(-1.0*Z),Z));
          YestEsp = rowsum(elem_prod(Wm_arti,cageEsp));

  NSurvey = elem_prod(elem_prod(No,mfexp(-1.0*Z*9/12)),outer_prod(uno_years,S_cru));

  for(t=1; t<=inxtrawl; t++)
    {
    pobsarr(t) = catagetrwl(strawl(t))/sum(catagetrwl(strawl(t)+tiny));
    pestarr(t) = cageArr(strawl(t))/sum(cageArr(strawl(t)));
    }

  for(t=1; t<=inxpal; t++)
    {
    pobspal(t) = catagelongl(spal(t))/sum(catagelongl(spal(t)+tiny));
    pestpal(t) = cagePal(spal(t))/sum(cagePal(spal(t)));
    }

  for(t=1; t<=inxesp; t++)
    {
    pobsesp(t) = catageartisa(sesp(t))/sum(catageartisa(sesp(t)+tiny));
    pestesp(t) = cageEsp(sesp(t))/sum(cageEsp(sesp(t)));
    }

  for(t=1; t<=inxsurv; t++)
    {
    pobssurv(t) = natagesurvey(ssurv(t))/sum(natagesurvey(ssurv(t)+tiny));
    pestsurv(t) = NSurvey(ssurv(t))/sum(NSurvey(ssurv(t)));
    }


FUNCTION biomass_state
  int t;

  BT = rowsum(elem_prod(No,Wm_pop));
  for (t=styr; t<=endyr; t++)
   {
   B6(t) = sum(elem_prod(No(t)(stage+5, endage),Wm_pop(t)(stage+5, endage)));
   Bdepl(t)=SB(t)/So;
  }


FUNCTION evaluate_objective_function
  int t;

  logL.initialize();

  logL(1) = -1.*nss(1)*sum(elem_prod(pobsarr,log(pestarr)));
  logL(2) = -1.*nss(2)*sum(elem_prod(pobspal,log(pestpal)));
  logL(3) = -1.*nss(3)*sum(elem_prod(pobsesp,log(pestesp)));
  logL(4) = -1.*nss(4)*sum(elem_prod(pobssurv,log(pestsurv)));

   for (t=styr; t<=endyr; t++)
    {
    if (indxtrawl(t)>0)
       {logL(5) += 0.5 * square(log(indxtrawl(t))-log(CPUEarr(t)))/square(cv_arr(t));}

    if (indxlongline(t)>0)
       {logL(6) += 0.5 * square(log(indxlongline(t))-log(CPUEpal(t)))/square(cv_pal(t));}

    if (surveyindex(t)>0)
      {logL(7) += 0.5 * square(log(surveyindex(t))-log(estBDcru(t)))/square(cv_cru(t));}

    if (indxartisanal(t)>0)
       {logL(8) += 0.5 * square(log(indxartisanal(t))-log(CPUEesp(t)))/square(cv_art(t));}

    if (ytrawl(t)>0)
       {logL(9) += 0.5 * square(log(ytrawl(t))-log(YestArr(t)))/square(cv_s(1));}

    if (ylongline(t)>0)
       {logL(10) += 0.5 * square(log(ylongline(t))-log(YestPal(t)))/square(cv_s(2));}

    if (yartisanal(t)>0)
       {logL(11) += 0.5 * square(log(yartisanal(t))-log(YestEsp(t)))/square(cv_s(3));}

    }

  logL(12) = 0.5 * norm2(mu)/square(cv_s(4));//  + size_count(mu)*log(cvs(7));

  // palangre arrastre espinel crucero
  penL(1) = norm2(d_pal - 0.5)/(2*square(cv_p));
  penL(2) = norm2(d_arr - 0.5)/(2*square(cv_p));
  penL(3) = norm2(d_esp - 0.5)/(2*square(cv_p));
  penL(4) = norm2(d_cru - 0.5)/(2*square(cv_p));

  penL(5) = 0.5 * (square(mfexp(log_selA(1))-rango_sa(1))/(2*square(cv_sel_a)));
  penL(6) = 0.5 * (square(mfexp(log_selA(2))-rango_sa(2))/(2*square(cv_sel_a)));
  penL(7) = 0.5 * (square(mfexp(log_selA(3))-rango_sa(3))/(2*square(cv_sel_a)));
  penL(8) = 0.5 * (square(mfexp(log_selA(4))-rango_sa(4))/(2*square(cv_sel_a)));

  penL(9) = 0.5 * (square(mfexp(log_selC(1))-rango_sl(1))/(2*square(10.59*cv_sel_c)));
  penL(10) = 0.5 * (square(mfexp(log_selC(2))-rango_sl(2))/(2*square(10.59*cv_sel_c)));
  penL(11) = 0.5 * (square(mfexp(log_selC(3))-rango_sl(3))/(2*square(cv_sel_c)));
  penL(12) = 0.5 * (square(mfexp(log_selC(4))-rango_sl(4))/(2*square(cv_sel_c)));

  objF = sum(logL) + penL(9) + penL(10) + penL(11) + penL(12);


FUNCTION sim_Fcte
  int j,t;

  for (j=1; j<=nFt; j++)
	{
		Np(endyr) = No(endyr);
    Sp(endyr) = Surv(endyr);
    Fp(endyr) = Fcr_total(endyr);
    Zp(endyr) = Z(endyr);
    NSp(endyr) = NS(endyr);
    Sep(endyr) = Fcr_total(endyr)/max(Fcr_total(endyr));
    caep(endyr) = cageArr(endyr) + cagePal(endyr) + cageEsp(endyr);

    Ftp(endyr) = Ftot(endyr);
    SBp(endyr) = SB(endyr);
    Rp(endyr)  = R(endyr);
    Yp(endyr)  = YestArr(endyr) + YestPal(endyr) + YestEsp(endyr);

		msp = msex(endyr);
		wp  = Wm_pop(endyr);

		for (t=endyr; t<=endyr+yr_sim-1; t++)
		{
      if (recOp == 1){
        Rp(t+1) = mean(R(endyr-5,endyr));
      } else if (recOp == 2) {
        Rp(t+1) = mean(R(endyr-26,endyr-12));
      } else {
        Rp(t+1) = mfexp(log_Ro-0.5*square(cv_s(4)));
      }

      Np(t+1,1) = Rp(t+1);
      Np(t+1)(stage+1,endage) =  ++elem_prod(Np(t)(stage, endage - 1),Sp(t)(stage, endage - 1));
      Np(t+1,endage) += Np(t+1,endage)/(1-Sp(t,endage));

      if (t==endyr){
        Fp(t)    = Ftp(t)*Sep(t);
        Ftp(t+1) = rCap*Ftp(t);
        Sep(t+1) = Sep(t);
        Fp(t+1)  = Ftp(t+1)*Sep(t+1);
      } else {
        Ftp(t+1) = mf(j)*Fpbr;
        Sep(t+1) = Sep(endyr);
        Fp(t+1) = Ftp(t+1)*Sep(t+1);
      }

      Zp(t+1) = Fp(t+1) + M;
      Sp(t+1) = mfexp(-1.0 * Zp(t+1));

      NSp(t+1) = elem_prod(elem_prod(elem_prod(Np(t+1),msp),wp),mfexp(-1.0*Zp(t+1)*9/12));
      SBp(t+1) = sum(extract_row(NSp,t+1));

      caep(t+1) = elem_prod(elem_div(Fp(t+1),Zp(t+1)),elem_prod(1.0-Sp(t+1),Np(t+1)));
			Yp(t+1) 	= sum(elem_prod(caep(t+1),wp));
      Yproy(endyr,j) = Yp(endyr);
			Yproy(t+1,j) 	 = Yp(t+1);
      Fproy(endyr,j) = Ftp(endyr);
			Fproy(t+1,j)	 = Ftp(t+1);
      BDp(endyr,j)   = SBp(endyr);
      BDp(t+1,j)   = SBp(t+1);
		}
	}

    if(mceval_phase())
    {
    ofstream pry("model_msur_ctp2017.pry.mcmc",ios::app);
    for (int i=endyr+1; i<=endyr+yr_sim; i++)
    {
      pry << "Captura " << i << Yproy(i) << endl;
    }
    pry.close();
    }


FUNCTION Write_proj
  int k, i;
  adstring report_name;
    {
      report_name = "proyecciones.jcq";
      ofstream R_report(report_name);

      R_report << "SBdeple" << endl; for (int i=styr; i<=endyr; i++){
        double lb=value(SBdeple(i)/exp(2.*sqrt(log(1+square(SBdeple.sd(i))/square(SBdeple(i))))));
        double ub=value(SBdeple(i)*exp(2.*sqrt(log(1+square(SBdeple.sd(i))/square(SBdeple(i))))));
        R_report << i << " " << SBdeple(i) << " " << SBdeple.sd(i) << " " << lb << " " << ub << endl;
      }

      R_report << "sprdeple" << endl; for (int i=styr; i<=endyr; i++){
        double lb=value(sprdeple(i)/exp(2.*sqrt(log(1+square(sprdeple.sd(i))/square(sprdeple(i))))));
        double ub=value(sprdeple(i)*exp(2.*sqrt(log(1+square(sprdeple.sd(i))/square(sprdeple(i))))));
        R_report << i << " " << sprdeple(i) << " " << sprdeple.sd(i) << " " << lb << " " << ub << endl;
      }

      R_report << "Ftot" << endl; for (int i=styr; i<=endyr; i++){
        double lb=value(Ftot(i)/exp(2.*sqrt(log(1+square(Ftot.sd(i))/square(Ftot(i))))));
        double ub=value(Ftot(i)*exp(2.*sqrt(log(1+square(Ftot.sd(i))/square(Ftot(i))))));
        R_report << i << " " << Ftot(i) << " " << Ftot.sd(i) << " " << lb << " " << ub << endl;
      }

      R_report << "spr" << endl;
      R_report << spr << " " << spr.sd << endl;

      R_report << "Ro" << endl;
        R_report << Ro << " " << Ro.sd << endl;

      R_report << "So" << endl;
        R_report << So << " " << So.sd << endl;


      R_report<<"Catch_fut"<< endl;
      for (k=1; k<=nFt; k++){
        for (i=endyr; i<=endyr+yr_sim; i++)
        {
          if (k==1){
            R_report << k << " " << i << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
          } else {
            double lb=value(Yproy(i,k)/exp(2.*sqrt(log(1+square(Yproy.sd(i,k))/square(Yproy(i,k))))));
            double ub=value(Yproy(i,k)*exp(2.*sqrt(log(1+square(Yproy.sd(i,k))/square(Yproy(i,k))))));
            R_report << k << " " << i << " " << Yproy(i,k) << " " << Yproy.sd(i,k) << " " << lb << " " << ub << endl;
          }
        }
      }

      // R_report<<"Ft_fut"<< endl;
      // for (k=1; k<=nFt; k++){
      //   for (i=endyr; i<=endyr+2; i++)
      //   {
      //     if (k==1){
      //       R_report << k << " " << i << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
      //     } else {
      //       double lb=value(Fproy(i,k)/exp(2.*sqrt(log(1+square(Fproy.sd(i,k))/square(Fproy(i,k))))));
      //       double ub=value(Fproy(i,k)*exp(2.*sqrt(log(1+square(Fproy.sd(i,k))/square(Fproy(i,k))))));
      //       R_report << k << " " << i << " " << Fproy(i,k) << " " << Fproy.sd(i,k) << " " << lb << " " << ub << endl;
      //     }
      //   }
      // }

      R_report<<"BD_fut"<< endl;
      for (k=1; k<=nFt; k++){
        for (i=endyr; i<=endyr+yr_sim; i++)
        {
          if (k==1){
            R_report << k << " " << i << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
          } else {
            double lb=value(BDp(i,k)/exp(2.*sqrt(log(1+square(BDp.sd(i,k))/square(BDp(i,k))))));
            double ub=value(BDp(i,k)*exp(2.*sqrt(log(1+square(BDp.sd(i,k))/square(BDp(i,k))))));
            R_report << k << " " << i << " " << BDp(i,k) << " " << BDp.sd(i,k) << " " << lb << " " << ub << endl;
          }
        }
      }

      R_report.close();
    }


REPORT_SECTION
  reporte(vages);
  reporte(vyears);
  reporte(PeT);
  reporte(PeP);
  reporte(PeE);
  reporte(strawl);
  reporte(spal);
  reporte(sesp);
  reporte(ssurv);
  reporte(logL);
  reporte(penL);
  reporte(S_arr);
  reporte(S_pal);
  reporte(S_esp);
  reporte(S_cru);
  reporte(surveyindex);
  reporte(indxtrawl);
  reporte(indxlongline);
  reporte(indxartisanal);
  reporte(CPUEarr);
  reporte(CPUEpal);
  reporte(CPUEesp);
  reporte(estBDcru);
  reporte(ytrawl);
  reporte(ylongline);
  reporte(yartisanal);
  reporte(YestArr);
  reporte(YestPal);
  reporte(YestEsp);
  reporte(pobsarr);
  reporte(pestarr);
  reporte(pobspal);
  reporte(pestpal);
  reporte(pobsesp);
  reporte(pestesp);
  reporte(pobssurv);
  reporte(pestsurv);
  reporte(mu);
  reporte(R);
  reporte(SB);
  reporte(Bdepl);
  reporte(BDcru);
  reporte(BMVpal);
  reporte(BMVesp);
  reporte(BMVarr);
  reporte(BT);
  reporte(B6);
  reporte(Farr);
  reporte(Fpal);
  reporte(Fesp);
  reporte(Fcr_total);
  reporte(Yproy);
  reporte(BDp);
  reporte(mf);
  reporte(Ftot);
  reporte(alpha);
  reporte(beta);
  reporte(watagefleet_pop);
  reporte(watagefleet_trawl);
  reporte(watagefleet_long);
  reporte(watagefleet_arti);
  reporte(watagefleet_surv);
  reporte(No);
  reporte(So);
  reporte(Fproy);
  reporte(Fp);
  reporte(msex);
  reporte(Surv);
  reporte(surveyindex);
  reporte(cpueindex);
  reporte(catagetrwl);
  reporte(catagelongl);
  reporte(catageartisa);
  reporte(natagesurvey);
  reporte(Wm_pop);
  reporte(Wm_trawl);
  reporte(Wm_long);
  reporte(Wm_arti);
  reporte(Wm_surv);
  reporte(yr_sim);
  reporte(chQarr);
  reporte(chQpal);
  reporte(SBdeple);
  reporte(offset);
  reporte(offsetCt);
  reporte(logL);
  reporte(Np);
  reporte(Sp);
  reporte(Fp);
  reporte(Fpbr);


FINAL_SECTION
  Write_proj();
