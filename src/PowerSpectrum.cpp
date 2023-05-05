#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // DONE: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector k_array = Utils::linspace(k_min,k_max,n_k);
  Vector log_k_array = log(k_array);


  //=========================================================================
  // DONE: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // DONE: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // #TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  Vector z_array = Utils::linspace(0,40000,40001);
  Vector bessel;
  for (int l=0;l < ells.size();l++){
    int ell = ells[l];
    bessel.clear();
    for (int z_index=0;z_index < z_array.size();z_index++){
      int z = z_array[z_index];
      bessel.push_back(Utils::j_ell(ell, z));
    }
    j_ell_splines[l].create(z_array,bessel);
  }
  //std::cout << "success" << std::endl;
  //j_ell(z) = Utils::j_ell(ell, z)
    
  //=============================================================================
  // DONE: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

//  for(size_t i = 0; i < ells.size(); i++){
//    const int ell = ells[i];
//
//    // ...
//    // ...
//    // ...
//    // ...
//
//    // Make the j_ell_splines[i] spline
//  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");
  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  const double eta_0 = cosmo->eta_of_x(0);
  //const double dx = 2.*M_PI/6.;
  //const double n_x = 10./dx 
  //Vector x_array = Utils::linspace(-10,0,100)
 
  //k*(eta_0-eta)


  double eta_p0;
  double eta_p1;
  double dx;
  for (int l = 0;l<ells.size();l++){
    double ell = ells[l];
    for(size_t ik = 0; ik < k_array.size(); ik++){
      double k = k_array[ik];
    
      //=============================================================================
      // TODO: Implement to solve for the general line of sight integral 
      // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
      // given value of k
      //=============================================================================
      // ...
      // ...
      
      //integrate from -10 to 0...
      double area = 0.0;
      //trapezoid rule
      for (int x_index = 0;x_index<x_array.size()-1;x_index++){
        dx = x_array[x_index+1] - x_array[x_index];
        eta_p0 = cosmo->eta_of_x(x_array[x_index]);
        eta_p1 = cosmo->eta_of_x(x_array[x_index+1]);
        area += (source_function(x_array[x_index+1],k)*j_ell_splines[l](k*(eta_0-eta_p1)))+(source_function(x_array[x_index],k)*j_ell_splines[l](k*(eta_0-eta_p0)))/2.*dx;

        //A = (a+b)/2*dx

      }
      
      //F_ell(k) = Int dx jell(k(eta  -eta0)) * S(x,k)
      // Store the result for Source_ell(k) in results[ell][ik]
      result[l][ik] = area;
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // DONE: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  // ...
  // ...
  // ...
  // ...
  for (int i=0;i<nells;i++){
    thetaT_ell_of_k_spline[i].create(k_array,thetaT_ell_of_k[i]); 
  }
  //============================================================================
  // #TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){

    // ...
    // ...
    // ...
    // ...

  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
    Vector k_array = exp(log_k_array);
  const int nells      = ells.size();
  Vector result;
  double area = 0.;
  double du;
  for (int l =0;l<nells;l++){
    for (int index = 0;index<log_k_array.size()-1;index++){
      du = log_k_array[index+1]-log_k_array[index];
      area += (primordial_power_spectrum(k_array[index+1])*f_ell_spline[l](k_array[index+1])*f_ell_spline[l](k_array[index+1])
            + primordial_power_spectrum(k_array[index])*f_ell_spline[l](k_array[index])*f_ell_spline[l](k_array[index]))/2.*du;
    }
  result.push_back(4.*M_PI*area);
  }
  //============================================================================
  // DONE: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  // ...
  // ...
  // ...
  // ...

  

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_j_ell(const double z,const int ell) const{
  return j_ell_splines[ell](z);
}
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

