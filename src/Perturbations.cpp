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

  //compute_source_functions();
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
  arma::vec k_array = arma::logspace(log10(k_min),log10(k_max),n_k);


  
  //Vector k_array = {0.001/Constants.Mpc,0.01/Constants.Mpc,0.1/Constants.Mpc};
  //double n_k = 3;
  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);
    
    //checking index to see how large the sub vector of x_array should be
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
    tight_coupling_ODE.solve(dydx_tight_coupling,tight_coupling_x_array,y_tight_coupling_ini);
    auto tight_coupling_data = tight_coupling_ODE.get_data();
    
    //we only have 7 parameters in tight coupling
    Vector y_tight_coupling_end(7);
    for (int i=0;i < 7;i++){
      y_tight_coupling_end[i] = tight_coupling_data[tight_coupling_data[i].size()-1][i];
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
    full_ODE.solve(dydx_full,from_and_after_tc_x_array,y_full_ini);
    auto after_tc_data = full_ODE.get_data();

    //removing the overlap at x_end_tight
    for (int i=0;i < tight_coupling_data.size();i++){
     tight_coupling_data[i].pop_back();
    }


    //storing the data
    //tight_coupling_data[0].size()+after_tc_data[0].size() = x_array.size() = const.
    //x_end_tight only changes with k

    //delta_cdm
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //    delta_cdm_vector.push_back(tight_coupling_data[i][0]);
    //    }
    //    else{
    //      delta_cdm_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][0]);
    //    }
    //}
    
    //v_cdm
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      v_cdm_vector.push_back(tight_coupling_data[i][1]);
    //    }
    //    else{
    //      v_cdm_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][1]);
    //    }
    //}

    //delta_b
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      delta_b_vector.push_back(tight_coupling_data[i][2]);
    //    }
    //    else{
    //      delta_b_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][2]);
    //    }
    //}

    //v_b
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      v_b_vector.push_back(tight_coupling_data[i][3]);
    //    }
    //    else{
    //      v_b_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][3]);
    //    }
    //}

    //Phi
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      Phi_vector.push_back(tight_coupling_data[i][4]);
    //    }
    //    else{
    //      Phi_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][4]);
    //    }
    //}

    //theta0
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta0_vector.push_back(tight_coupling_data[i][5]);
    //    }
    //    else{
    //      theta0_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][5]);
    //    }
    //}

     //theta1
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta1_vector.push_back(tight_coupling_data[i][6]);
    //    }
    //    else{
    //      theta1_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][6]);
    //    }
    //}

    //all other thetas are zero during tight coupling
    
    //theta2
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta2_vector.push_back(0.);
    //    }
    //    else{
    //      theta2_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][7]);
    //    }
    //}


    //theta3
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta3_vector.push_back(0.);
    //    }
    //    else{
    //      theta3_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][8]);
    //    }
    //}
//
    ////theta4
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta4_vector.push_back(0.);
    //    }
    //    else{
    //      theta4_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][9]);
    //    }
    //}
//
    ////theta5
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta5_vector.push_back(0.);
    //    }
    //    else{
    //      theta5_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][10]);
    //    }
    //}
    //      
//
    ////theta6
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta6_vector.push_back(0.);
    //    }
    //    else{
    //      theta6_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][11]);
    //    }
    //}
//
    ////theta7
    //for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
    //    if (i<tight_coupling_data[0].size()){
    //      theta7_vector.push_back(0.);
    //    }
    //    else{
    //      theta7_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][12]);
    //    }
    //}

    //all quantities
    for (int i=0;i<tight_coupling_data[0].size()+after_tc_data[0].size();i++){
        if (i<tight_coupling_data[0].size()){
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
          delta_cdm_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][0]);
          v_cdm_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][1]);
          delta_b_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][2]);
          v_b_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][3]);
          Phi_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][4]);

          theta0_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][5]);
          theta1_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][6]);

          theta2_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][7]);
          theta3_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][8]);
          theta4_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][9]);
          theta5_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][10]);
          theta6_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][11]);
          theta7_vector.push_back(after_tc_data[i-tight_coupling_data[0].size()][12]);
          
        }
    }



    
    
    
    
    //ready to be splined:
    
    
    //Phi_vector is
    //for(int ix = 0; ix < nx; ix++)
    //for(int iy = 0; iy < ny; iy++)
    //  z_array[ix + nx * iy] = function(x_array[ix], y_array[iy]);

    //Vector delta_cdm;
    //for 


    //Vector v_cdm;
    //Vector delta_b;
    //Vector v_b;
    //Vector Phi;
    //Vector 

  



   //    tight_coupling_data[0].push_back(after_tc_data[0][0]);

    //combining the results into one solution array
    //  for(int i = 0; i < all_data.size(); i++){
    //  auto y_exact = analytical_solution(y_ic, x_array[i]);
    //  std::cout << " xi     = " << std::setw(10) << x_array[i];
    //  std::cout << " y0(xi) = " << std::setw(10) << all_data[i][0]  << " y0_exact(xi) = " << std::setw(10) << y_exact[0];
    //  std::cout << " y1(xi) = " << std::setw(10) << all_data[i][1]  << " y1_exact(xi) = " << std::setw(10) << y_exact[1] << "\n";
    //}


    //double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
    // double &delta_b      =  y_tc[Constants.ind_deltab_tc];
    // double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
    // double &v_b          =  y_tc[Constants.ind_vb_tc];
    // double &Phi          =  y_tc[Constants.ind_Phi_tc];
    



    //removing the solution at x_end_tight_coupling since there is overlap with the full solution at this point
    //theta_l with l >= 2 is zero
  

       //for (int j=0; j < after_tc_data[j].size();j++){
       //   tight_coupling_data[i].push_back(after_tc_data[i][j])
       //   }
    

    //for (int j=0; j < after_tc_data[j].size();j++){
    //    after_tc_data[j].push_back()
    //
    //
    //




    //Vector combined_delta_b;
    //Vector combined_v_cdm;
    //Vector combined_v_b;
    //Vector combined_Phi;

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
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
  if (k > 0.15/Constants.Mpc){
    for (double time_index = 0; time_index < x_array.size(); time_index++){
        if (rec->dtaudx_of_x(x_array[time_index]) > 10. && rec->dtaudx_of_x(x_array[time_index]) > 10.*Constants.c*k/cosmo->Hp_of_x(x_array[time_index])){}
        
        else{
          if (x_array[time_index] <= -8.3){
            return x_array[time_index];
          }
          else{
            return -8.3;
          }
        }
    } 
  }
  
  else{
    return -8.3;
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
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

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
  double q = (-((1.-R)*dtau_dx+(1.+ R)*ddtau_ddx)*(3.*theta1+v_b)-ck_Hp*Psi+(1.-dHp_dx/Hp)*ck_Hp*(-theta0+2.*theta2)-ck_Hp*dtheta0dx)
             /((1. + R)*dtau_dx+dHp_dx/Hp-1.);
  
    
  // SET: Scalar quantities (Phi, delta, v, ...)
  dv_bdx = 1./(1. + R)*(-v_b-ck_Hp*Psi+R*(q+ck_Hp*(-theta0+2.*theta2)-ck_Hp*Psi));
  
  
  dPhidx = Psi-ck_Hp*ck_Hp*Phi/3.+ H0*H0/(2.*Hp*Hp)*(OmegaCDM0/a*delta_cdm+Omegab0/a*delta_b+4.*OmegaR0/(a*a)*theta0);
  ddelta_cdmdx = ck_Hp*v_cdm-3.*dPhidx;
  dv_cdmdx = -v_cdm-ck_Hp*Psi;
  ddelta_bdx = ck_Hp*v_b-3.*dPhidx;
  
  
  // SET: Photon multipoles (Theta_ell)
  dtheta1dx = 1./3.*(q-dv_bdx);
  dtheta0dx = -ck_Hp*theta1-dPhidx;
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
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
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
    fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}