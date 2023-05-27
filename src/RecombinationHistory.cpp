#include"RecombinationHistory.h"
//#include<iostream>

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
 
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
  solve_number_density_electrons_with_separate_solutions();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();

  //sound horizon
  sound_horizon();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // DONE: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start,x_today,npts_rec_arrays); 
  
  //linspaces from 0 to 1 because I just want a vector with npts_rec_arrays elements
  Vector Xe_arr  = Utils::linspace(0.,1.,npts_rec_arrays);
  Vector ne_arr  = Utils::linspace(0.,1.,npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // DONE: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      //std::cout << Xe_current << std::endl;
      saha_regime = false;
      

    if(saha_regime){
      //std::cout << "saha" << std::endl; 
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
      //=============================================================================
      // DONE: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...

    } else {

      //==============================================================
      // DONE: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // DONE: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
      Vector Xe_init{Xe_current};
      const double OmegaB0      = cosmo->get_OmegaB(0.0);
      const double G           = Constants.G;

      const double rho_c0 = 3.*pow(cosmo->get_H0(),2.)/(8.*M_PI*G);
      const double m_H         = Constants.m_H;

      //solving the ODE on the remaining x_array called sub_x
      Vector sub_x = {x_array.begin() + i,x_array.end()};   
      peebles_Xe_ode.solve(dXedx, sub_x, Xe_init);
      auto Xe_peebles = peebles_Xe_ode.get_data_by_component(0);
      // the for loop is used to put the ODE solution into the correct place in the pre existing solution array
      for (int j=i; j < npts_rec_arrays;j++){
        Xe_arr[j] = Xe_peebles[j-i];        
        double nb = OmegaB0*rho_c0/(m_H*exp(3.*x_array[j]));
        ne_arr[j] = Xe_peebles[j-i]*nb;

      }
      break;
    }
  }
  //taking the logarithm of ne since it varies a lot (for the spline). Xe does not vary to the same extent...
  Vector log_ne = Utils::linspace(0.,1,npts_rec_arrays); 

  for (int i=0;i < npts_rec_arrays;i++){
     log_ne[i] = log(ne_arr[i]);
  }
  
  Xe_of_x_spline.create(x_array,Xe_arr,"Xe");
  log_ne_of_x_spline.create(x_array,log_ne,"ne");
  //std::cout<< Xe_of_x_spline(-6.9859707941588) << std::endl;


  //=============================================================================
  // DONE: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...

  Utils::EndTiming("Xe");
}
void RecombinationHistory::solve_number_density_electrons_with_separate_solutions(){
  Utils::StartTiming("Xe separate solutions");
  
  //=============================================================================
  // DONE: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start,x_today,npts_rec_arrays); 
  
  //linspaces from 0 to 1 because I just want a vector with npts_rec_arrays elements
  Vector Xe_arr_saha  = Utils::linspace(0.,1.,npts_rec_arrays);
  Vector ne_arr_saha  = Utils::linspace(0.,1.,npts_rec_arrays);

  //Vector Xe_arr_peebels  = Utils::linspace(0.,1.,npts_rec_arrays);
  Vector ne_arr_peebles  = Utils::linspace(0.,1.,npts_rec_arrays);

  // Calculate recombination history
  //bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){
    //std::cout << i << std::endl;

    //==============================================================
    // DONE: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;
    
    Xe_arr_saha[i] = Xe_current;
    ne_arr_saha[i] = ne_current;
  }

    // Are we still in the Saha regime?
    //if(Xe_current < Xe_saha_limit)
    //  //std::cout << Xe_current << std::endl;
    //  saha_regime = false;
      

    //if(saha_regime){
    //  //std::cout << "saha" << std::endl; 
    //  Xe_arr[i] = Xe_current;
    //  ne_arr[i] = ne_current;
    //  //=============================================================================
    //  // DONE: Store the result we got from the Saha equation
    //  //=============================================================================
    //  //...
    //  //...

    

      //==============================================================
      // DONE: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...
      //std::cout << "202" << std::endl;
      // The Peebles ODE equation
//      ODESolver peebles_Xe_ode;
//      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
//        return rhs_peebles_ode(x, Xe, dXedx);
//      };
//      //std::cout << "208" << std::endl;
//      //=============================================================================
//      // DONE: Set up IC, solve the ODE and fetch the result 
//      //=============================================================================
//      //...
//      //...
//      //Vector Xe_init{Xe_current};
//      const double OmegaB0      = cosmo->get_OmegaB(0.0);
//      const double G           = Constants.G;
//
//      const double rho_c0 = 3.*pow(cosmo->get_H0(),2.)/(8.*M_PI*G);
//      const double m_H         = Constants.m_H;
//
//      Vector init;
//      init.push_back(Xe_arr_saha[0]);  std::cout << init[0] << std::endl;
//      peebles_Xe_ode.solve(dXedx, x_array, init); //std::cout << "223" << std::endl;
//      auto Xe_peebles = peebles_Xe_ode.get_data_by_component(0);
//      for (int i=0; i < npts_rec_arrays;i++){
//        double nb = OmegaB0*rho_c0/(m_H*exp(3.*x_array[i]));
//        ne_arr_peebles[i] = Xe_peebles[i]*nb;
//
//      }
    
  
  //taking the logarithm of ne since it varies a lot (for the spline). Xe does not vary to the same extent...
  Vector log_ne_saha = Utils::linspace(0.,1,npts_rec_arrays); 
  //Vector log_ne_peebles = Utils::linspace(0.,1,npts_rec_arrays); 

  for (int i=0;i < npts_rec_arrays;i++){
     log_ne_saha[i] = log(ne_arr_saha[i]);
     //log_ne_peebles[i] = log(ne_arr_peebles[i]);
  }
  
  Xe_of_x_saha_spline.create(x_array,Xe_arr_saha,"Xe_saha");
  //std::cout<< "Xe_spline here" << std::endl;
  log_ne_of_x_saha_spline.create(x_array,log_ne_saha,"ne_saha");

  //Xe_of_x_peebles_spline.create(x_array,Xe_peebles,"Xe_peebles");
  //log_ne_of_x_peebles_spline.create(x_array,log_ne_peebles,"ne_peebles");


  //=============================================================================
  // DONE: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...

  Utils::EndTiming("Xe separate solutions");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB0      = cosmo->get_OmegaB(0.0);
  
  
  double TCMB      = cosmo->get_TCMB(x);
  const double rho_c0 = 3.*pow(cosmo->get_H0(),2.)/(8*M_PI*G);

  // Electron fraction and number density
  double Xe;// = 0.0;
  double ne;// = 0.0;
  
  //=============================================================================
  // DONE: Compute Xe and ne from the Saha equation
  //=============================================================================
  //...
  //...

  double nb = OmegaB0*rho_c0/(m_H*a*a*a);
  double b = 1./nb*pow(m_e*TCMB*k_b/(2.*M_PI*hbar*hbar),3./2.)*exp(-epsilon_0/(TCMB*k_b));
  
  // rewriting the square root to avoid Huge - Huge.
  Xe = 2./(sqrt(1.+4./b)+1.);
  ne = Xe*nb; //nH = nb

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){
  //std::cout << x << std::endl;
  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double rho_c0 = 3.*pow(cosmo->get_H0(),2.)/(8*M_PI*G);
  const double OmegaB0      = cosmo->get_OmegaB(0.0);
  double nb = OmegaB0*rho_c0/(m_H*a*a*a);
  double H = cosmo->H_of_x(x);
  double TCMB = cosmo->get_TCMB(x);
  

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...

  //=============================================================================
  // DONE: Write the expression for dXedx
  //=============================================================================
  //...
  //...

  
  double phi_2 = 0.448*log(epsilon_0/(TCMB*k_b)); //fixed
  
  double alpha_2 = c*8./sqrt(3.*M_PI) * sigma_T * sqrt(epsilon_0/(TCMB*k_b))*phi_2; //fixed?
  
  double beta = alpha_2*pow(m_e*TCMB*k_b/(2.*M_PI*hbar*hbar),3./2.)*exp(-epsilon_0/(TCMB*k_b)); //fixed
  
  double beta_2 = alpha_2*pow(m_e*TCMB*k_b/(2.*M_PI*hbar*hbar),3./2.)*exp(-1./4.*epsilon_0/(TCMB*k_b)); //fixed
  
  double n_1s = (1. - X_e)*nb; //fixed
  
  double lambda_a = H * pow(3.*epsilon_0/(hbar*c), 3.)/(pow(8.*M_PI,2.)*n_1s);//fixed
  
  double Cr = (lambda_2s1s + lambda_a)/(lambda_2s1s + lambda_a + beta_2); //fixed
  
  double rhs = Cr/H * (beta*(1.- X_e) - nb*alpha_2*X_e*X_e); //fixed

  
  dXedx[0] = rhs;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1e3;
  Vector x_array = Utils::linspace(x_start, x_today, npts);
  Vector x_array_reverse = Utils::linspace(x_today, x_start, npts);
 
  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    const double H = cosmo->H_of_x(x);
    const double ne = ne_of_x(x);
    //=============================================================================
    // DONE: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = -c*ne*sigma_T/H;


    return GSL_SUCCESS;
  };

  //=============================================================================
  // DONE: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...
  ODESolver tau_of_x_ode;
  Vector tau_init = {0.0};
  tau_of_x_ode.solve(dtaudx,x_array_reverse,tau_init);
  auto tau_of_x_solution = tau_of_x_ode.get_data_by_component(0);
  
  Vector tau_of_x_solution_increasing = Utils::linspace(0.,1.,npts);

  
  for (int i =0;i<npts;i++){
    tau_of_x_solution_increasing[i] = tau_of_x_solution[npts-1-i];
  }
  tau_of_x_spline.create(x_array,tau_of_x_solution_increasing,"tau");

  //=============================================================================
  // DONE: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...
  Vector g_tilde_array = Utils::linspace(0.,1.,npts);
  for (int i=0;i<npts;i++){

    g_tilde_array[i] = g_tilde(x_array[i]);
  }
  g_tilde_of_x_spline.create(x_array,g_tilde_array,"g");
  
  

  Utils::EndTiming("opticaldepth");
}

double RecombinationHistory::g_tilde(double x) const{

    return -dtaudx_of_x(x)*exp(-tau_of_x(x));
  }


//sound_horizon
void RecombinationHistory::sound_horizon(){
   //may be overkill with x_very_very_early = -100
    Vector x_array = Utils::linspace(x_start,x_today,npts_rec_arrays);
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    const double OmegaR0 = cosmo->get_OmegaR(0.0);
    const double OmegaB0 = cosmo->get_OmegaB(0.0);
    ODEFunction dsdx = [&](double x, const double *s, double *dsdx){  

      const double a = exp(x);
      const double Hp = cosmo->Hp_of_x(x);
      double R = 4.*OmegaR0/(3.*OmegaB0*a);
      double cs = c*sqrt(R/(3.*(1.+R)));

    // Set the derivative for photon optical depth
      dsdx[0] = cs/Hp;
      return GSL_SUCCESS;
    };
  double R_init = 4.*OmegaR0/(3.*OmegaB0*exp(x_start));
  double cs_init = c*sqrt(R_init/(3.*(1.+R_init)));
  double Hp_init = cosmo->Hp_of_x(x_start);

  ODESolver sound_horizon_ode;
  Vector sound_horizon_init = {cs_init/Hp_init};
  sound_horizon_ode.solve(dsdx,x_array,sound_horizon_init);
  auto sound_horizon_solution = sound_horizon_ode.get_data_by_component(0);
  sound_horizon_spline.create(x_array,sound_horizon_solution,"SH");

}

//====================================================
// Get methods
//====================================================


double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}
double RecombinationHistory::get_sound_horizon(double x) const{
  return sound_horizon_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // DONE: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...
  
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return Xe_of_x_spline(x);
}
double RecombinationHistory::Xe_of_x_peebles(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return Xe_of_x_peebles_spline(x);
}
double RecombinationHistory::Xe_of_x_saha(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return Xe_of_x_saha_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info(){
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:                                            " << Yp                                                    << "\n";
  std::cout << "Sound Horizon [Gpc]:                           " << get_sound_horizon(x_decoupling)/(Constants.Mpc*1e3)               << "\n";
  std::cout << "Sound Horizon over conformal time:             " << get_sound_horizon(x_decoupling)/cosmo->eta_of_x(x_decoupling)     << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_today;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp <<std::setprecision(15) <<Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp <<std::setprecision(15) <<Xe_of_x_saha(x)           << " ";
    //fp <<std::setprecision(15) <<Xe_of_x_peebles(x)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

