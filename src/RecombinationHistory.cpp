#include"RecombinationHistory.h"
//#include<iostream>

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) : //, double rho_c0, double OmegaB0) :
  cosmo(cosmo),
  Yp(Yp)//,
  //rho_c0(rho_c0),
  //OmegaB0(OmegaB0)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // DONE: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start,x_end,npts_rec_arrays); 
  
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

      const double rho_c0 = 3.*pow(cosmo->get_H0(),2.)/(8*M_PI*G);
      const double m_H         = Constants.m_H;


      //peebles_Xe_ode.solve(dXedx, x_array, Xe_init);
      Vector sub_x = {x_array.begin() + i,x_array.end()};   
      peebles_Xe_ode.solve(dXedx, sub_x, Xe_init);
      auto Xe_peebles = peebles_Xe_ode.get_data_by_component(0);
      for (int j=i; j < npts_rec_arrays;j++){
        Xe_arr[j] = Xe_peebles[j-i];
        
        double nb = OmegaB0*rho_c0/(m_H*exp(3.*x_array[j]));

        ne_arr[j] = Xe_peebles[j-i]*nb;

      }


       
      break;
    
    }
  }
  //auto t_array = ode_for_t_of_x.get_data_by_component(0);
  //Vector log_Xe = Utils::linspace(0.,1,npts_rec_arrays); 
  Vector log_ne = Utils::linspace(0.,1,npts_rec_arrays); 

  for (int i=0;i<npts_rec_arrays;i++){
  log_ne[i] = log(ne_arr[i]);
  //log_Xe[i] =log(Xe_arr[i]);
  //std::cout << std::setprecision(15) << ne_arr[i] << std::endl;
  }
  
  Xe_of_x_spline.create(x_array,Xe_arr,"Xe");
  log_ne_of_x_spline.create(x_array,log_ne,"ne");
  //Xe_of_x_spline.create(log_Xe_of_x_spline))  


  //=============================================================================
  // DONE: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...

  Utils::EndTiming("Xe");
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
  
  //Xe = (-b+pow(b*b+4*b,1./2.))/2.;
  Xe = 2./(sqrt(1.+4./b)+1.);
  ne = Xe*nb; //nH = nb

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

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
  const int npts = 4000;
  Vector x_array = Utils::linspace(x_start, x_today, npts);
  Vector x_array_reverse = Utils::linspace(x_today, x_start, npts);
  Vector z_array = Utils::linspace(0.,1.,npts);
 
  for (int i =0;i<npts;i++){
  //std::cout << x_array_reverse[i] << std::endl;
  z_array[i] = exp(-x_array_reverse[i])-1.;
  std::cout << z_array[i] << std::endl;
  }

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudz = [&](double z, const double *tau, double *dtaudz){
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    const double x = -log(1+z);
    const double H = cosmo->H_of_x(x);
    const double ne = ne_of_x(x);
    //const double z = exp(-x)-1.; 
    //=============================================================================
    // DONE: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth

    dtaudz[0] = c*ne*sigma_T/(H*(1.+z));

    return GSL_SUCCESS;
  };

  //=============================================================================
  // DONE: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...
  ODESolver tau_of_z_ode;
  Vector tau_init = {0.0};
  tau_of_z_ode.solve(dtaudz,z_array,tau_init);
  auto tau_of_z_solution = tau_of_z_ode.get_data_by_component(0);
  tau_of_z_spline.create(z_array,tau_of_z_solution,"tau_z");
  Vector tau_of_x_solution = Utils::linspace(0.,1.,npts);
  
  for (int i =0;i<npts;i++){
    tau_of_x_solution[i] = tau_of_z_spline(-log(1.+z_array[i]));
  }
  //tau_of_x_spline.create(x_array_reverse,tau_of_x_solution,"tau");


   



  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // DONE: Implement
  //=============================================================================
  //...
  //...

  return Xe_of_x_spline(x);
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
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp <<std::setprecision(15) <<Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    //fp << dtaudx_of_x(x)       << " ";
    //fp << ddtauddx_of_x(x)     << " ";
    //fp << g_tilde_of_x(x)      << " ";
    //fp << dgdx_tilde_of_x(x)   << " ";
    //fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

