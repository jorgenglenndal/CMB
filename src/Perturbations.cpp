#include"Perturbations.h"



//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  //relevant for next milestone

  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");
  
  //===================================================================
  // DONE: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  arma::vec k_array_arma = arma::logspace(log10(k_min),log10(k_max),n_k);
  
  Vector k_array;
  for (int i=0;i<k_array_arma.size();i++){
      k_array.push_back(k_array_arma[i]);
  }
  

  
 
  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }
    
    // Current value of k
    double k = k_array[ik];
    ////double k = 10./Constants.Mpc;
    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);
    
    //checking index to see how large the sub vector of x_array should be, i.e. in tight coupling
    int tight_index;
    for (int i =0;i < x_array.size();i++){
      if (x_array[i] > x_end_tight){
        tight_index = i;
        break;  
      }
      else{}
    }
    Vector tight_coupling_x_array;
    tight_coupling_x_array.clear();
    for (int i=0;i<tight_index;i++){
        tight_coupling_x_array.push_back(x_array[i]);
    }
    Vector from_and_after_tc_x_array;
    from_and_after_tc_x_array.clear();
    for (int i=0;i < n_x-tight_index+1;i++){
      from_and_after_tc_x_array.push_back(x_array[tight_index-1+i]);
    }
    //there is overlap at x_end_tight. The last value in tc is therefore removed later... 

    //===================================================================
    // DONE: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================




    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);
    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx_tight_coupling){
      return rhs_tight_coupling_ode(x, k, y, dydx_tight_coupling);
    };
    // Integrate from x_start -> x_end_tight
    ODESolver tight_coupling_ODE;
    std::vector<Vector> tight_coupling_data;
    tight_coupling_data.clear();
    //Vector tight_coupling_data;
    if (tight_coupling_x_array.size() >= 2){
        tight_coupling_ODE.solve(dydx_tight_coupling,tight_coupling_x_array,y_tight_coupling_ini);
        tight_coupling_data = tight_coupling_ODE.get_data();
    }
    else{
      tight_coupling_data.push_back(y_tight_coupling_ini);
    }
    
    
    //we only have 7 parameters in tight coupling
    Vector y_tight_coupling_end(7);
    for (int i=0;i < 7;i++){
      y_tight_coupling_end[i] = tight_coupling_data.back()[i];
    }


    

    //====i===============================================================
    // DONE: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling_end, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx_full){
      return rhs_full_ode(x, k, y, dydx_full);
    };
    ODESolver full_ODE;
    std::vector<Vector> after_tc_data;
    after_tc_data.clear();
    if (from_and_after_tc_x_array.size() >= 2){
      full_ODE.solve(dydx_full,from_and_after_tc_x_array,y_full_ini);
      after_tc_data = full_ODE.get_data();
    }
    else{
      if (from_and_after_tc_x_array.size() == 1){
        //after_tc_data[0] = tight_coupling_data.back();
        after_tc_data.push_back(tight_coupling_data.back());
      }
      else{
        //nothing
      }
    }
    //removing the overlap at x_end_tight
    tight_coupling_data.pop_back();  
    
    //storing the data in correct format for spline
    //tight_coupling_data.size() + after_tc_data.size() = x_array.size() = const. 
    //all solved quantities
    for (int i = 0; i < tight_coupling_data.size() + after_tc_data.size(); i++){
        if (i<tight_coupling_data.size()){
          delta_cdm_vector.push_back(tight_coupling_data[i][0]);
          v_cdm_vector.push_back(tight_coupling_data[i][1]);
          delta_b_vector.push_back(tight_coupling_data[i][2]);
          v_b_vector.push_back(tight_coupling_data[i][3]);
          Phi_vector.push_back(tight_coupling_data[i][4]);
          
          theta0_vector.push_back(tight_coupling_data[i][5]);
          theta1_vector.push_back(tight_coupling_data[i][6]);
          
          //higher order thetas are zero during tight coupling
          theta2_vector.push_back(0.);
          theta3_vector.push_back(0.);
          theta4_vector.push_back(0.);
          theta5_vector.push_back(0.);
          theta6_vector.push_back(0.);
          theta7_vector.push_back(0.);
        }
        else{
          delta_cdm_vector.push_back(after_tc_data[i-tight_coupling_data.size()][0]);
          v_cdm_vector.push_back(after_tc_data[i-tight_coupling_data.size()][1]);
          delta_b_vector.push_back(after_tc_data[i-tight_coupling_data.size()][2]);
          v_b_vector.push_back(after_tc_data[i-tight_coupling_data.size()][3]);
          Phi_vector.push_back(after_tc_data[i-tight_coupling_data.size()][4]);

          theta0_vector.push_back(after_tc_data[i-tight_coupling_data.size()][5]);
          theta1_vector.push_back(after_tc_data[i-tight_coupling_data.size()][6]);

          theta2_vector.push_back(after_tc_data[i-tight_coupling_data.size()][7]);
          theta3_vector.push_back(after_tc_data[i-tight_coupling_data.size()][8]);
          theta4_vector.push_back(after_tc_data[i-tight_coupling_data.size()][9]);
          theta5_vector.push_back(after_tc_data[i-tight_coupling_data.size()][10]);
          theta6_vector.push_back(after_tc_data[i-tight_coupling_data.size()][11]);
          theta7_vector.push_back(after_tc_data[i-tight_coupling_data.size()][12]);
        }
    }

    //storing Psi
    
    //using the initial condition
    Psi_vector.push_back(-2./3.);
    const double OmegaR0 = cosmo->get_OmegaR(0.);
    const double Omegab0 = cosmo->get_OmegaB(0.);
    const double OmegaCDM0 = cosmo->get_OmegaCDM(0.);
    const double H0 = cosmo->get_H0();

    for (int i = 1; i < x_array.size(); i++){
      double a = exp(x_array[i]);
      Psi_vector.push_back(-Phi_vector[i+n_x * ik]-12.*H0*H0/(Constants.c*Constants.c*k*k*a*a)*OmegaR0*theta2_vector[i + n_x*ik]);
      ////Psi_vector.push_back(-Phi_vector[i+n_x * 0]-12.*H0*H0/(Constants.c*Constants.c*k*k*a*a)*OmegaR0*theta2_vector[i + n_x*0]);
    }


  }
  


  //only relevant when calculatong for single k value
  ////std::string filename; 
  ////std::ofstream ofile;
  ////filename = "safeguard___perturbations_k10.txt";        
  ////int width = 40;
  ////int prec = 12; 
  ////ofile.open(filename);
  ////for (int i = 0; i < x_array.size(); i++){
  ////    ofile << std::setw(width)<<std::setprecision(prec) << x_array[i] <<
  ////    std::setw(width) << std::setprecision(prec) << theta0_vector[i]  <<
  ////    std::setw(width)<<std::setprecision(prec) <<   theta1_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec)  << theta2_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec) << Phi_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec) << Psi_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec) << delta_cdm_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec) << delta_b_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec) << v_cdm_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec) << v_b_vector[i] <<
  ////    std::setw(width)<<std::setprecision(prec) << cosmo->eta_of_x(x_array[i]) <<
  ////    std::endl;      
  ////}
  ////ofile.close();



  //=============================================================================
  // DONE: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================

  delta_cdm_spline.create(x_array,k_array,delta_cdm_vector,"delta_cdm_spline");
  delta_b_spline.create(x_array,k_array,delta_b_vector,"delta_b_spline");
  v_cdm_spline.create(x_array,k_array,v_cdm_vector,"v_cdm_spline");
  v_b_spline.create(x_array,k_array,v_b_vector,"v_b_spline");
  Phi_spline.create(x_array,k_array,Phi_vector,"Phi_spline");
  Psi_spline.create(x_array,k_array,Psi_vector,"Psi_spline");

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  Theta_spline[0].create(x_array,k_array,theta0_vector,"theta0_spline");
  Theta_spline[1].create(x_array,k_array,theta1_vector,"theta1_spline");
  Theta_spline[2].create(x_array,k_array,theta2_vector,"theta2_spline");
  Theta_spline[3].create(x_array,k_array,theta3_vector,"theta3_spline");
  Theta_spline[4].create(x_array,k_array,theta4_vector,"theta4_spline");
  Theta_spline[5].create(x_array,k_array,theta5_vector,"theta5_spline");
  Theta_spline[6].create(x_array,k_array,theta6_vector,"theta6_spline");
  Theta_spline[7].create(x_array,k_array,theta7_vector,"theta7_spline");
  Utils::EndTiming("integrateperturbation");

}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  //double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  //double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //manually setting the indices for theta0 and theta1
  double &theta0       =  y_tc[Constants.ind_start_theta_tc];
  double &theta1       =  y_tc[Constants.ind_start_theta_tc+1];     
  //=============================================================================
  // DONE: Set the initial conditions in the tight coupling regime
  //=============================================================================

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  //value of Psi
  
  const double Psi = -2./3.;
  //Psi_vector.clear();
  //Psi_vector.push_back(Psi);
  


  //Psi_vector.push_back(-2./3.);
  Phi = -Psi;
  delta_cdm = -3./2.*Psi;
  delta_b = -3./2.*Psi;
  v_cdm = -Constants.c*k/(2.*cosmo->Hp_of_x(x))*Psi;
  v_b = -Constants.c*k/(2.*cosmo->Hp_of_x(x))*Psi;
  
  // SET: Photon temperature perturbations (Theta_ell)
  theta0 = -1./2.*Psi;
  theta1 = Constants.c*k/(6.*cosmo->Hp_of_x(x))*Psi;
  
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double &theta0_tc        = y_tc[Constants.ind_start_theta_tc];
  const double &theta1_tc        = y_tc[Constants.ind_start_theta_tc+1];
  //const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double &theta0           = y[Constants.ind_start_theta_tc];
  double &theta1           = y[Constants.ind_start_theta_tc + 1];
  double &theta2           = y[Constants.ind_start_theta_tc + 2];
  double &theta3           = y[Constants.ind_start_theta_tc + 3];
  double &theta4           = y[Constants.ind_start_theta_tc + 4];
  double &theta5           = y[Constants.ind_start_theta_tc + 5];
  double &theta6           = y[Constants.ind_start_theta_tc + 6];
  double &theta7           = y[Constants.ind_start_theta_tc + 7];
  //double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  //double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // DONE: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  // ...
  // ...
  double a = exp(x);
  const double OmegaR0 = cosmo->get_OmegaR(0.);
  const double Omegab0 = cosmo->get_OmegaB(0.);
  const double OmegaCDM0 = cosmo->get_OmegaCDM(0.);
  const double H0 = cosmo->get_H0();
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_ddx = rec->ddtauddx_of_x(x);
  double Hp = cosmo->Hp_of_x(x); 
  double dHp_dx = cosmo->dHpdx_of_x(x);
  double ck_Hp = Constants.c*k/Hp;
  
  



  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;
  Phi = Phi_tc;


  // SET: Photon temperature perturbations (Theta_ell)
  theta0 = theta0_tc;
  theta1 = theta1_tc;
  theta2 = -20./45.*ck_Hp/dtau_dx*theta1_tc;


  for (int i = 3; i < n_ell_theta;i++){
    y[Constants.ind_start_theta_tc+i] = -(i/(2.*i+1.))*ck_Hp/dtau_dx*y[Constants.ind_start_theta_tc+i-1];
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k){
    for (double time_index = 0; time_index < x_array.size(); time_index++){
        if (abs(rec->dtaudx_of_x(x_array[time_index])) > 10. && abs(rec->dtaudx_of_x(x_array[time_index])) > 10.*Constants.c*k/cosmo->Hp_of_x(x_array[time_index])){}
        
        else{
          if (x_array[time_index] <= -8.3){
            return x_array[time_index];
          }
          else{
            return -8.3;
          }
        }
    } 
    //=============================================================================
  // DONE: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  // ...
  std::cout<< "tight coupling time is wrong" << std::endl;
  return 0;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  //arma::vec k_array_arma = arma::logspace(log10(k_min),log10(k_max),n_k);
  //
  //Vector k_array;
  //for (int i=0;i<k_array_arma.size();i++){
  //    k_array.push_back(k_array_arma[i]);
  //}
  //k_array.push_back(k_max);
  Vector k_array = Utils::linspace(k_min,k_max,n_k);
  //double k_min = k_min*Mpc;
  //double k_max = k_max*Mpc;

  //for (int i=0;i<n_k;i++){
  //  k_array.push_back(k_min_+(k_max_-k_min_)*(i/(n_k-1))*(i/(n_k-1)));
  //  std::cout <<  k_array[i] << std::endl;
  //}
  //k_array.push_back(k_max);
  
  Vector x_array = Utils::linspace(-10.,x_end,round(1e3));
  
  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array;//(k_array.size() * x_array.size());
  ST_array.clear();
  //Vector SE_array(k_array.size() * x_array.size());
  double dx;
  //const double c             = Constants.c;
  const double c = Constants.c;
  //Vector x_array_source;
  // Compute source functions

  for(auto ik = 0; ik < k_array.size(); ik++){
    const double k = k_array[ik];
    double ck = c*k;
    //x_array_source.clear();
    //x_array_source.push_back(-10.);
    for(auto ix = 0; ix < x_array.size(); ix++){
      const double x = x_array[ix];
      const double Hp            = cosmo->Hp_of_x(x);
      const double tau           = rec->tau_of_x(x);
      const double g_tilde       = rec->g_tilde_of_x(x);
      
      const double dHpdx         = cosmo->dHpdx_of_x(x);
      const double dg_dxtilde    = rec->dgdx_tilde_of_x(x);
      const double ddg_ddxtilde  = rec->ddgddx_tilde_of_x(x);
      const double ddHpddx       = cosmo->ddHpddx_of_x(x);
      const double theta2        = get_Theta(x,k,2);
      const double dtheta2dx     = get_dTheta2dx(x,k);
      const double ddtheta2ddx   = get_ddTheta2ddx(x,k);

      //dx = 2.*M_PI*cosmo->Hp_of_x(x)/(10.*ck);
      //if (x+dx >= 0.){
      //    x_array_source.push_back(0.);
      //    ST_array.push_back(g_tilde*(get_Theta(x,k,0)+get_Psi(x,k)+1./4.*get_Theta(x,k,2)) + exp(-tau)*(get_dPsidx(x,k)-get_dPhidx(x,k))
      //     -1./(c*k)*(dHpdx*g_tilde*get_v_b(x,k)+Hp*dg_dxtilde*get_v_b(x,k)+Hp*g_tilde*get_dv_bdx(x,k))+ 3./(4.*c*c*k*k)
      //     *((dHpdx*dHpdx + Hp*ddHpddx)*g_tilde*theta2 + 3.*Hp*dHpdx*(dg_dxtilde*theta2 + g_tilde*dtheta2dx) + Hp*Hp*(ddg_ddxtilde*theta2 + 2.*dg_dxtilde*dtheta2dx + g_tilde*ddtheta2ddx)));
      //    break;
      //}
      //x_array_source.push_back(x+dx);





        // NB: This is the format the data needs to be stored 
        // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      //const int index = ix + n_x * ik;

        //=============================================================================
        // TODO: Compute the source functions
        //=============================================================================
        // Fetch all the things we need...

         
         //const double dx            = 2.*M_PI*Hp/(c*k)



      // ...
      // ...

      // Temperatur source using splines for the derivatives
      ST_array.push_back(g_tilde*(get_Theta(x,k,0)+get_Psi(x,k)+1./4.*get_Theta(x,k,2)) + exp(-tau)*(get_dPsidx(x,k)-get_dPhidx(x,k))
           -1./(c*k)*(dHpdx*g_tilde*get_v_b(x,k)+Hp*dg_dxtilde*get_v_b(x,k)+Hp*g_tilde*get_dv_bdx(x,k))+ 3./(4.*c*c*k*k)
           *((dHpdx*dHpdx + Hp*ddHpddx)*g_tilde*theta2 + 3.*Hp*dHpdx*(dg_dxtilde*theta2 + g_tilde*dtheta2dx) + Hp*Hp*(ddg_ddxtilde*theta2 + 2.*dg_dxtilde*dtheta2dx + g_tilde*ddtheta2ddx)));
           //(dHpdx*(dHpdx*g_tilde*theta2+Hp*dg_dxtilde*theta2+Hp*g_tilde*dtheta2dx) + Hp*(ddHpddx*g_tilde*theta2+dHpdx*dg_dxtilde*theta2+dHpdx*g_tilde*dtheta2dx
           //+dHpdx*dg_dxtilde*theta2+Hp*ddg_ddxtilde*theta2+Hp*dg_dxtilde*dtheta2dx+dHpdx*g_tilde*dtheta2dx+Hp*dg_dxtilde*dtheta2dx+Hp*g_tilde*ddtheta2ddx));     
      //std::cout << ST_array[ix+n_x*ik] << std::endl; 
      // Polarization source
      //if(Constants.polarization){
      //  //SE_array[index] = 0.0;
      //}
    }
  }
  //if (x_array_source[x_array_source.size()-1] < 0.){
  //      std::cout << "x did not reach 0 in perturbation souurce" << std::endl;
  //    }
  // Spline the source functions
  //std::cout << "565" << std::endl;
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  //std::cout << "584" << std::endl;
  //if(Constants.polarization){
  //  SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  //}

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double &theta0          =  y[Constants.ind_start_theta_tc];
  const double &theta1          =  y[Constants.ind_start_theta_tc+1];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double &dtheta0dx       =  dydx[Constants.ind_start_theta_tc];
  double &dtheta1dx       =  dydx[Constants.ind_start_theta_tc+1];

  //=============================================================================
  // DONE: fill in the expressions for all the derivatives
  //=============================================================================
  double a = exp(x);
  const double OmegaR0 = cosmo->get_OmegaR(0.);
  const double Omegab0 = cosmo->get_OmegaB(0.);
  const double OmegaCDM0 = cosmo->get_OmegaCDM(0.);
  const double H0 = cosmo->get_H0();
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_ddx = rec->ddtauddx_of_x(x);
  double Hp = cosmo->Hp_of_x(x); 
  double dHp_dx = cosmo->dHpdx_of_x(x);
  double ck_Hp = Constants.c*k/Hp;
  double R = 4.*OmegaR0/(3.*Omegab0*a);
  double theta2 = -20./45.*ck_Hp/dtau_dx*theta1;
  double Psi = -Phi-12.*H0*H0/(Constants.c*Constants.c*k*k*a*a)*OmegaR0*theta2;
  

  dPhidx = Psi-ck_Hp*ck_Hp*Phi/3.+ H0*H0/(2.*Hp*Hp)*(OmegaCDM0/a*delta_cdm+Omegab0/a*delta_b+4.*OmegaR0/(a*a)*theta0);
  dtheta0dx = -ck_Hp*theta1-dPhidx;
  double q = (-((1.-R)*dtau_dx+(1.+ R)*ddtau_ddx)*(3.*theta1+v_b)-ck_Hp*Psi+(1.-dHp_dx/Hp)*ck_Hp*(-theta0+2.*theta2)-ck_Hp*dtheta0dx)
             /((1. + R)*dtau_dx+dHp_dx/Hp-1.);
  
    
  // SET: Scalar quantities (Phi, delta, v, ...)
  //dv_bdx = 1./(1. + R)*(-v_b-ck_Hp*Psi+R*(q+ck_Hp*(-theta0+2.*theta2)-ck_Hp*Psi));
  dv_bdx = -v_b-ck_Hp*Psi+R*(1./(1.+R))*(q + v_b + ck_Hp*(-theta0+2.*theta2));
  
  
  ddelta_cdmdx = ck_Hp*v_cdm-3.*dPhidx;
  dv_cdmdx = -v_cdm-ck_Hp*Psi;
  ddelta_bdx = ck_Hp*v_b-3.*dPhidx;
  
  
  // SET: Photon multipoles (Theta_ell)
  dtheta1dx = 1./3.*(q-dv_bdx);
  
  // ...
  // ...

  //Psi_vector.push_back(Psi);

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  //std::cout<<x << std::endl;
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const int n_ell_tot_full      = Constants.n_ell_tot_full;    
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double &theta0          =  y[Constants.ind_start_theta_tc];
  const double &theta1          =  y[Constants.ind_start_theta_tc + 1];
  const double &theta2          =  y[Constants.ind_start_theta_tc + 2];
  const double &theta3          =  y[Constants.ind_start_theta_tc + 3];
  const double &theta4          =  y[Constants.ind_start_theta_tc + 4];
  const double &theta5          =  y[Constants.ind_start_theta_tc + 5];
  const double &theta6          =  y[Constants.ind_start_theta_tc + 6];
  const double &theta7          =  y[Constants.ind_start_theta_tc + 7];  
  
  //const double *Theta_p         = &y[Constants.ind_start_thetap];
  //const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double &dtheta0_dx      = dydx[Constants.ind_start_theta_tc];
  double &dtheta1_dx      = dydx[Constants.ind_start_theta_tc + 1];
  double &dtheta2_dx      = dydx[Constants.ind_start_theta_tc + 2];
  double &dtheta3_dx      = dydx[Constants.ind_start_theta_tc + 3];
  double &dtheta4_dx      = dydx[Constants.ind_start_theta_tc + 4];
  double &dtheta5_dx      = dydx[Constants.ind_start_theta_tc + 5];
  double &dtheta6_dx      = dydx[Constants.ind_start_theta_tc + 6];
  double &dtheta7_dx      = dydx[Constants.ind_start_theta_tc + 7];
  
  
  //double *dThetadx        = &dydx[Constants.ind_start_theta];
  //double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  //double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  // 
  // ...

  double a = exp(x);
  const double OmegaR0 = cosmo->get_OmegaR(0.);
  const double Omegab0 = cosmo->get_OmegaB(0.);
  const double OmegaCDM0 = cosmo->get_OmegaCDM(0.);
  const double H0 = cosmo->get_H0();
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_ddx = rec->ddtauddx_of_x(x);
  double Hp = cosmo->Hp_of_x(x); 
  double dHp_dx = cosmo->dHpdx_of_x(x);
  double eta_x = cosmo->eta_of_x(x);
  double ck_Hp = Constants.c*k/Hp;
  double R = 4.*OmegaR0/(3.*Omegab0*a);
  double Psi = -Phi-12.*H0*H0/(Constants.c*Constants.c*k*k*a*a)*OmegaR0*theta2;
  //Psi_vector.push_back(Psi);

  //=============================================================================
  // DONE: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx = Psi-ck_Hp*ck_Hp*Phi/3.+ H0*H0/(2.*Hp*Hp)*(OmegaCDM0/a*delta_cdm+Omegab0/a*delta_b+4.*OmegaR0/(a*a)*theta0);
  ddelta_cdmdx = ck_Hp*v_cdm-3.*dPhidx;
  dv_cdmdx = -v_cdm-ck_Hp*Psi;
  ddelta_bdx = ck_Hp*v_b-3.*dPhidx;
  dv_bdx = -v_b-ck_Hp*Psi+dtau_dx*R*(3.*theta1+v_b);
  
  // SET: Photon multipoles (Theta_ell)
  dtheta0_dx = -ck_Hp*theta1-dPhidx;
  dtheta1_dx = ck_Hp/3.*theta0-2.*ck_Hp/3.*theta2+ck_Hp/3.*Psi+ dtau_dx*(theta1+1./3.*v_b);
  double delta;
  for (int i = 2;i < n_ell_theta-1 ;i++){
    if (i == 2){
      delta = 1.;
    }
    else{
      delta = 0.;
    }
    dydx[Constants.ind_start_theta+i] = i*ck_Hp/(2.*i+1.)*y[Constants.ind_start_theta+i-1]-(i+1.)*ck_Hp/(2.*i+1.)*y[Constants.ind_start_theta+i+1]
                                        +dtau_dx*(y[Constants.ind_start_theta+i]-1./10.*theta2*delta); //PI = theta2 for master students
  }

  dtheta7_dx = ck_Hp*theta6 - Constants.c*(7.+1.)/(Hp*eta_x)*theta7 + dtau_dx*theta7;
 

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================



double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_dv_bdx(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
//double Perturbations::get_Source_E(const double x, const double k) const{
//  return SE_spline(x,k);
//}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}

double Perturbations::get_dTheta2dx(const double x, const double k) const{
  return Theta_spline[2].deriv_x(x,k);
}

double Perturbations::get_ddTheta2ddx(const double x, const double k) const{
  return Theta_spline[2].deriv_xx(x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================


void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    fp << cosmo->eta_of_x(x) << " ";

    //fp << get_Pi(x,k)        << " ";
    //fp << get_Source_T(x,k)  << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}