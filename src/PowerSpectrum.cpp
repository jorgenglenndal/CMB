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
  //arma::vec k_array_arma = arma::logspace(log10(k_min),log10(k_max),n_k);
  
  //Vector k_array;
  //for (int i=0;i<n_k;i++){
  //  k_array.push_back(k_min+(k_max-k_min)*(i/100.)*(i/100.));
  //}
  //Vector k_array = Utils::linspace(k_min,k_max,n_k);
  //for (int i=0;i<k_array_arma.size();i++){
  //    k_array.push_back(k_array_arma[i]);
  //}
  ////Vector k_array = Utils::linspace(k_min,k_max,n_k);
  //Vector log_k_array = log(k_array);


  //=========================================================================
  // DONE: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  ////generate_bessel_function_splines();

  //=========================================================================
  // DONE: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  Vector k_array_new = Utils::linspace(k_min,k_max,n_k);
  Vector k_array_zone_0 = Utils::linspace(k_min_zone_0,k_max_zone_0,n_k_zone_0);
  Vector k_array_zone_1 = Utils::linspace(k_min_zone_1,k_max_zone_1,n_k_zone_1);
  Vector k_array_zone_2 = Utils::linspace(k_min_zone_2,k_max_zone_2,n_k_zone_2);
  Vector k_array_zone_3 = Utils::linspace(k_min_zone_3,k_max_zone_3,n_k_zone_3);
  Vector k_array_zone_4 = Utils::linspace(k_min_zone_4,k_max_zone_4,n_k_zone_4);
  Vector k_array_zone_5 = Utils::linspace(k_min_zone_5,k_max_zone_5,n_k_zone_5);


  all_k_arrays.push_back(k_array_zone_0);
  all_k_arrays.push_back(k_array_zone_1);
  all_k_arrays.push_back(k_array_zone_2);
  all_k_arrays.push_back(k_array_zone_3);
  all_k_arrays.push_back(k_array_zone_4);
  all_k_arrays.push_back(k_array_zone_5);
  
  
  //arma::vec k_array_arma = arma::logspace(log10(k_min),log10(k_max),n_k);
  //for (int i=0;i<n_k;i++){
  //  k_array_new.push_back(k_array_arma[i]);
  //}
  line_of_sight_integration(k_array_new);
  
  
  
  //for (int i=0;i<1e6;i++){
  //    std::cout << "her1?" << std::endl;
  //}

  //=========================================================================
  // DONE: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  
  
  auto cell_TT = solve_for_cell(thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //cell_TT_spline.create(k_array_new, k_array_new, "Cell_TT_of_ell");
 
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  double dz = 2.*M_PI/10.;
  int n_z = round(k_max*cosmo->eta_of_x(0.)/dz);
  std::cout << n_z <<std::endl;

  Vector z_array = Utils::linspace(0.,k_max*cosmo->eta_of_x(0.),1e3);//round(dz*n_z));
  std::cout << z_array.size() <<std::endl;
  Vector bessel;
  for (int l=0;l < ells.size();l++){
    int ell = ells[l];
    bessel.clear();
    for (int z_index=0;z_index < z_array.size();z_index++){
      int z = z_array[z_index];
      bessel.push_back(Utils::j_ell(ell, z));
    }
    j_ell_splines[l].create(z_array,bessel,"bessel");
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
  Vector2D result;// = Vector2D(ells.size(), Vector(k_array.size()));
  Vector result_ell;
  const double eta_0 = cosmo->eta_of_x(0.);
  //const double dx = 2.*M_PI/6.;
  //const double n_x = 10./dx
  const double Mpc = Constants.Mpc; 
  //Vector x_array = Utils::linspace(-10,0,100)
 
  //k*(eta_0-eta)
 
  //double eta_p0;
  //double eta_p1;
  Vector k_array_in_use;
  double dx;
  double area;
  double k;
  double npts;
  for (int l = 0;l<ells.size();l++){
    double ell = ells[l];
    std::cout << "LOS" << std::endl;
    std::cout << ell << std::endl;
    if (ell <=10){
      k_array_in_use = all_k_arrays[0];
      npts = n_k_zone_0;
    }
    if (ell > 10 && ell<=30){
      k_array_in_use = all_k_arrays[1];
      npts = n_k_zone_1;
    }

    if (ell > 30 && ell<=100){
      k_array_in_use = all_k_arrays[2];
      npts = n_k_zone_2;
    }
    if (ell > 100 && ell<=300){
      k_array_in_use = all_k_arrays[3];
      npts = n_k_zone_3;
    }
    if (ell > 300 && ell<=700){
      k_array_in_use = all_k_arrays[4];
      npts = n_k_zone_4;
    }
    if (ell > 700){
      k_array_in_use = all_k_arrays[5];
      npts = n_k_zone_5;
    }

    for(int ik = 0; ik < npts; ik++){
      k = k_array_in_use[ik];
      //std::cout << k*Mpc << std::endl;
    

      
    
      //=============================================================================
      // DONE: Implement to solve for the general line of sight integral 
      // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
      // given value of k
      //=============================================================================
      // ...
      // ...
      
      //integrate from -10 to 0...
      // x_array defined in PoweSpectrum.h
      area = 0.0;
      //trapezoid rule
      double dx;
      double ck = Constants.c*k;
      Vector x_array_source;
      x_array_source.clear();
      x_array_source.push_back(-8.);
      double bessel_0;
      double bessel_1;
      for (int i = 0;i< 1e10;i++){
        bessel_0 = Utils::j_ell(ell, (k*(eta_0-cosmo->eta_of_x(x_array_source[i]))));
        ////bessel_0 = j_ell_splines[l](k*(eta_0-cosmo->eta_of_x(x_array_source[i])));
        dx = 2.*M_PI*cosmo->Hp_of_x(x_array_source[i])/(10.*10.*10.*ck);
        if (ell <= 10){
            dx = 2.*M_PI*cosmo->Hp_of_x(x_array_source[i])/(10.*10.*10.*10.*ck); //*10.*10.
        }
        if (ell > 100){
            dx = 2.*M_PI*cosmo->Hp_of_x(x_array_source[i])/(10.*10.*ck); //*10.*10.
        }
                       
        if (x_array_source[i]+dx >= 0.){
          x_array_source.push_back(0.);
          bessel_1 = Utils::j_ell(ell, (k*(eta_0-cosmo->eta_of_x(x_array_source[i+1]))));
          ////bessel_1 = j_ell_splines[l](k*(eta_0-cosmo->eta_of_x(x_array_source[i+1])));
          dx = 0.- x_array_source[i]; 
          //eta_p0 = cosmo->eta_of_x(x_array_source[i]);
          //eta_p1 = cosmo->eta_of_x(x_array_source[i+1]);
          area += ((source_function(x_array_source[i+1],k)*bessel_1)+(source_function(x_array_source[i],k)*bessel_0))/2.*dx;  
          
          break;
        }
        x_array_source.push_back(x_array_source[i]+dx);
        bessel_1 = Utils::j_ell(ell, (k*(eta_0-cosmo->eta_of_x(x_array_source[i+1]))));
        ////bessel_1 = j_ell_splines[l](k*(eta_0-cosmo->eta_of_x(x_array_source[i+1])));
        //eta_p0 = cosmo->eta_of_x(x_array_source[i]);
        //eta_p1 = cosmo->eta_of_x(x_array_source[i+1]);
        //std::cout << "her?" << std::endl;
        area += ((source_function(x_array_source[i+1],k)*bessel_1)+(source_function(x_array_source[i],k)*bessel_0))/2.*dx;
        }
      if (x_array_source[x_array_source.size()-1] < 0.){
        std::cout << "x did not reach 0 in LOS" << std::endl;
      }
      result_ell.push_back(area);
      //result[l][ik] = area;



//      for (int x_index = 0;x_index<x_array.size()-1;x_index++){
//        //dx = x_array[x_index+1] - x_array[x_index];
//        eta_p0 = cosmo->eta_of_x(x_array[x_index]);
//        eta_p1 = cosmo->eta_of_x(x_array[x_index+1]);
//        area += (source_function(x_array[x_index+1],k)*j_ell_splines[l](k*(eta_0-eta_p1)))+(source_function(x_array[x_index],k)*j_ell_splines[l](k*(eta_0-eta_p0)))/2.*dx;
//
//        //A = (a+b)/2*dx
//
//      }
      
      //F_ell(k) = Int dx jell(k(eta  -eta0)) * S(x,k)
      // Store the result for Source_ell(k) in results[ell][ik]
    }
    result.push_back(result_ell);
    result_ell.clear();

  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  //const int n          = 100;
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
  thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  // ...
  // ...
  // ...
  // ...
  for (int i=0;i<nells;i++){
    int ell = ells[i];
    //std::cout << "285" << std::endl;
    if (ell<=10){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[0],thetaT_ell_of_k[i]); 
    }
    if (ell > 10 && ell <= 30){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[1],thetaT_ell_of_k[i]); 
    }
    if (ell > 30 && ell <= 100){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[2],thetaT_ell_of_k[i]); 
    }
    if (ell > 100 && ell <= 300){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[3],thetaT_ell_of_k[i]); 
    }
    if (ell > 300 && ell <= 700){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[4],thetaT_ell_of_k[i]); 
    }
    if (ell > 700){
    thetaT_ell_of_k_spline[i].create(all_k_arrays[5],thetaT_ell_of_k[i]); 
    }
    //std::cout << i << std::endl;
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
  //Vector & log_k_array,
  std::vector<Spline> & f_ell_spline){
  //std::vector<Spline> & g_ell_spline
  
  //Vector k_array = exp(log_k_array);
  const int nells      = ells.size();
  const double eta_0 = cosmo->eta_of_x(0.);

  Vector result;
  result.clear();
  double area;
  //double du;
  //for (int l =0;l<nells;l++){
  //  
  //  for (int index = 0;index<log_k_array.size()-1;index++){
  //    du = log_k_array[index+1]-log_k_array[index];
  //    area += (primordial_power_spectrum(k_array[index+1])*f_ell_spline[l](k_array[index+1])*f_ell_spline[l](k_array[index+1])
  //          + primordial_power_spectrum(k_array[index])*f_ell_spline[l](k_array[index])*f_ell_spline[l](k_array[index]))/2.*du;
  //  }
  //result.push_back(4.*M_PI*area);
  //}
  //============================================================================
  // DONE: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  //std::cout << "268" << std::endl;
  double dk;
  double dlogk;
                                  ////arma::vec k_array_cell = arma::logspace(log10(k_min),log10(k_max),n_k);
  Vector k_array_in_use;// = Utils::linspace(k_min,k_max,1e3);
  //dk = 2.*M_PI/(10.*eta_0);
  //dk = 2.*M_PI/(10.*10*eta_0);
  //k_array_cell.clear();
  //k_array_cell.push_back(k_min);
  //for(int ik = 0; ik < 1e10; ik++){
  //  
  //  if (k_array_cell[ik]+dk >= k_max){
  //    k_array_cell.push_back(k_max);
  //    break;
  //    //dk = k_max- k_array_cell[ik];
  //  }
  //  k_array_cell.push_back(k_array_cell[ik]+dk);
  //}
  //std::cout << "k_array size in cell" <<std::endl;
  //std::cout << k_array_cell.size() <<std::endl;
  //if (k_array_cell[k_array_cell.size()-1] < k_max){
  //      std::cout << "k did not reach k_max in cell" << std::endl;
  //}
  double k_min_cell;
  double k_max_cell;
  arma::vec log_k_array_cell;// = log(k_array_cell);

  for (int l = 0;l<nells;l++){

    double ell = ells[l];
    if (ell <=10){
      k_min_cell = all_k_arrays[0][0];
      k_max_cell = all_k_arrays[0][all_k_arrays[0].size()-1];
      //k_array_in_use = all_k_arrays[0];
      //log_k_array_cell = log(k_array_in_use);
      //npts = n_k_zone_0;
    }
    if (ell > 10 && ell<=30){
      k_min_cell = all_k_arrays[1][0];
      k_max_cell = all_k_arrays[1][all_k_arrays[1].size()-1];
      //k_array_in_use = all_k_arrays[1];
      //log_k_array_cell = log(k_array_in_use);
      //npts = n_k_zone_1;
    }

    if (ell > 30 && ell<=100){
      k_min_cell = all_k_arrays[2][0];
      k_max_cell = all_k_arrays[2][all_k_arrays[2].size()-1];
      //k_array_in_use = all_k_arrays[2];
      //log_k_array_cell = log(k_array_in_use);
      //npts = n_k_zone_2;
    }
    if (ell > 100 && ell<=300){
      k_min_cell = all_k_arrays[3][0];
      k_max_cell = all_k_arrays[3][all_k_arrays[3].size()-1];
      //k_array_in_use = all_k_arrays[3];
      //log_k_array_cell = log(k_array_in_use);
      //npts = n_k_zone_3;
    }
    if (ell > 300 && ell<=700){
      k_min_cell = all_k_arrays[4][0];
      k_max_cell = all_k_arrays[4][all_k_arrays[4].size()-1];
      //k_array_in_use = all_k_arrays[4];
      //log_k_array_cell = log(k_array_in_use);
      //npts = n_k_zone_4;
    }
    if (ell > 700){
      k_min_cell = all_k_arrays[5][0];
      k_max_cell = all_k_arrays[5][all_k_arrays[5].size()-1];
      //k_array_in_use = all_k_arrays[5];
      //log_k_array_cell = log(k_array_in_use);
      //npts = n_k_zone_5;
    }
    k_array_in_use = Utils::linspace(k_min_cell,k_max_cell,1e5);
    log_k_array_cell = log(k_array_in_use);
    //k_array_cell = all_k_arrays[l];
    std::cout << "cell" << std::endl;
    std::cout << ell << std::endl;
    area = 0.0;
    for (int index = 0;index < k_array_in_use.size()-1;index++){
      dlogk = log_k_array_cell[index+1]-log_k_array_cell[index];
      //dk = k_array_cell[index+1] - k_array_cell[index];
      //area += (primordial_power_spectrum(k_array_cell[index+1])*pow(f_ell_spline[l](k_array_cell[index+1]),2)
      //      + primordial_power_spectrum(k_array_cell[index])*pow(f_ell_spline[l](k_array_cell[index]),2))/2.*dlogk;
      //std::cout << "420" << std::endl;
      area += (primordial_power_spectrum(k_array_in_use[index+1])*pow(f_ell_spline[l](k_array_in_use[index+1]),2)
            + primordial_power_spectrum(k_array_in_use[index])*pow(f_ell_spline[l](k_array_in_use[index]),2))/2.*dlogk;    
      //std::cout << "423" << std::endl;   
    }


  result.push_back(4.*M_PI*area);
    //for(size_t ik = 0; ik < 1e10; ik++){
    //  dk = 2.*M_PI/(10.*eta_0);
    //  if (k_array_cell[ik]+dk >= k_max){
    //    k_array_cell.push_back(k_max);
    //    //dk = k_max- k_array_cell[ik];
    //  }
    //  k_array_cell.push_back(k_array_cell[ik]+dk);
    //}

  
  }
  //std::cout << "311" << std::endl;

      
      //area += (primordial_power_spectrum(k_array_cell[ik+1])*f_ell_spline[l](k_array_cell[ik+1])*f_ell_spline[l](k_array_cell[ik+1])
      //      + primordial_power_spectrum(k_array_cell[ik])*f_ell_spline[l](k_array_cell[ik])*f_ell_spline[l](k_array_cell[ik]))/2.*du;

      

      //double k = k_array[ik];
    
      //=============================================================================
      // DONE: Implement to solve for the general line of sight integral 
      // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
      // given value of k
      //=============================================================================
      
      //trapezoid rule
      //double dx;
      //double ck = Constants.c*k;
      //Vector x_array_source;
      //x_array_source.clear();
      //x_array_source.push_back(-10.);
      //for (int i = 0;i< 1e10;i++){
      //  dx = 2.*M_PI*cosmo->Hp_of_x(x_array_source[i])/(10.*ck);
      //    if (x_array_source[i]+dx >= 0.){
      //      x_array_source.push_back(0.);
      //      dx = 0.- x_array_source[i]; 
      //      eta_p0 = cosmo->eta_of_x(x_array_source[i]);
      //      eta_p1 = cosmo->eta_of_x(x_array_source[i+1]);
      //      area += (source_function(x_array_source[i+1],k)*j_ell_splines[l](k*(eta_0-eta_p1)))+(source_function(x_array_source[i],k)*j_ell_splines[l](k*(eta_0-eta_p0)))/2.*dx;  
      //      
      //      break;
      //    }
      //    else{}
      //  x_array_source.push_back(x_array_source[i]+dx);
      //  eta_p0 = cosmo->eta_of_x(x_array_source[i]);
      //  eta_p1 = cosmo->eta_of_x(x_array_source[i+1]);
      //  area += (source_function(x_array_source[i+1],k)*j_ell_splines[l](k*(eta_0-eta_p1)))+(source_function(x_array_source[i],k)*j_ell_splines[l](k*(eta_0-eta_p0)))/2.*dx;
      //  }
      //
      //result[l][ik] = area;



//      for (int x_index = 0;x_index<x_array.size()-1;x_index++){
//        //dx = x_array[x_index+1] - x_array[x_index];
//        eta_p0 = cosmo->eta_of_x(x_array[x_index]);
//        eta_p1 = cosmo->eta_of_x(x_array[x_index+1]);
//        area += (source_function(x_array[x_index+1],k)*j_ell_splines[l](k*(eta_0-eta_p1)))+(source_function(x_array[x_index],k)*j_ell_splines[l](k*(eta_0-eta_p0)))/2.*dx;
//
//        //A = (a+b)/2*dx
//
//      }
      
      //F_ell(k) = Int dx jell(k(eta  -eta0)) * S(x,k)
      // Store the result for Source_ell(k) in results[ell][ik]
    //}
  

  return result;
}


//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3 # SI
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  //double k_ = k*Constants.Mpc/cosmo->get_h();

  double delta_m = Constants.c*Constants.c*k*k*pert->get_Phi(x,k)*2./(3.*(cosmo->get_OmegaB(0.)+cosmo->get_OmegaCDM(0.))*cosmo->get_H0()*cosmo->get_H0())*exp(x);
  double pofk = pow(abs(delta_m),2)*primordial_power_spectrum(k)*2.*M_PI*M_PI/(k*k*k);

  //=============================================================================
  // DONE: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
//double PowerSpectrum::get_j_ell(const double z,const int ell) const{
//  return j_ell_splines[ell](z);
//}
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


void PowerSpectrum::output_matter(const std::string filename) const{
  const double k_min_ = k_min;
  const double k_max_ =  k_max;
  const int    n_pts =  10000;
  const double Mpc = Constants.Mpc;
  const double conversion_P = pow(cosmo->get_h(),3)/pow(Mpc,3);
  const double conversion = Mpc/cosmo->get_h();
  Vector k_array = Utils::linspace(k_min_, k_max_, n_pts);
  //Vector k_array_scaled;
  //for (int i=0;i<k_array.size();i++){
  //  k_array_scaled.push_back(k_array[i]*Mpc/cosmo->get_h());
  //}
   //= k_array*Mpc/cosmo->get_h();
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double k) {
    fp <<    k*conversion                                     << " ";
    fp << get_matter_power_spectrum(0.,k)*conversion_P        << " ";
    //fp << Hp_of_x(x)                               << " ";
    //fp << dHpdx_of_x(x)                            << " ";
    //fp << get_OmegaB(x)                            << " ";
    //fp << get_OmegaCDM(x)                          << " ";
    //fp << get_OmegaLambda(x)                       << " ";
    //fp << get_OmegaR(x)                            << " ";
    //fp << get_OmegaNu(x)                           << " ";
    //fp << get_OmegaK(x)                            << " ";
    //fp << t_of_x(x)                                << " ";
    //fp << get_luminosity_distance_of_x(x)          << " ";
    //fp << ddHpddx_of_x(x)                          << " ";

    fp <<"\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output_Source(const std::string filename) const{
  const double x_min_ = -8.;
  const double x_max_ =  0.;
  const int    n_pts =  10000;
  double k = 340.*cosmo->get_H0()/Constants.c;
  const double Mpc = Constants.Mpc;
  Vector x_array = Utils::linspace(x_min_,x_max_,n_pts);
  //const double conversion_P = pow(cosmo->get_h(),3)/pow(Mpc,3);
  //const double conversion = Mpc/cosmo->get_h();
  //Vector k_array = Utils::linspace(k_min_, k_max_, n_pts);
  //Vector k_array_scaled;
  //for (int i=0;i<k_array.size();i++){
  //  k_array_scaled.push_back(k_array[i]*Mpc/cosmo->get_h());
  //}
   //= k_array*Mpc/cosmo->get_h();
  double eta_0 = cosmo->eta_of_x(0.);
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp <<    x                                                                           << " ";
    fp << pert->get_Source_T(x,k)*Utils::j_ell(2,(k*(eta_0-cosmo->eta_of_x(x))))* 1e3  << " ";
    //fp << pert->get_Source_T(x,k)*get_j_ell(k*(eta_0-cosmo->eta_of_x(x)),0)* 1e3    << " ";
    fp << pert->get_Source_T(x,k)                                                        << " ";
    //fp << get_j_ell(k*(eta_0-cosmo->eta_of_x(x)),0)                                 << " ";
    //fp << Utils::j_ell(100,(k*(eta_0-cosmo->eta_of_x(x))))                               << " ";
    //fp << j_ell_splines[19](k*(eta_0-cosmo->eta_of_x(x)))* 1e3<< " ";
    //fp << j_ell_splines[19](k*(eta_0-cosmo->eta_of_x(x)))* 1e3<< " ";      
          
    //fp <<        Utils::j_ell(100, (k*(eta_0-cosmo->eta_of_x(x))))*1e3 << " ";

    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

void PowerSpectrum::output_45(const std::string filename) const{
  //const double k_min_ = k_min;
  //const double k_max_ = k_max;
  const int    n_pts =  10000;
  //double k = 340.*cosmo->get_H0()/Constants.c;
  const double Mpc = Constants.Mpc;
  //const double c_H0 = Constants.c/cosmo->get_H0();
  const double c_H0 = Constants.c/cosmo->get_H0();
  const double norm = 1e6*cosmo->get_H0()/Constants.c;
  arma::vec k_array = arma::logspace(log10(k_min),log10(k_max),n_pts);
  //const double conversion_P = pow(cosmo->get_h(),3)/pow(Mpc,3);
  //const double conversion = Mpc/cosmo->get_h();
  //Vector k_array = Utils::linspace(k_min_, k_max_, n_pts);
  //Vector k_array_scaled;
  //for (int i=0;i<k_array.size();i++){
  //  k_array_scaled.push_back(k_array[i]*Mpc/cosmo->get_h());
  //}
   //= k_array*Mpc/cosmo->get_h();
  //double eta_0 = cosmo->eta_of_x(0.);
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double k) {
    fp <<    k*c_H0                                     << " ";
    fp << thetaT_ell_of_k_spline[0](k)*thetaT_ell_of_k_spline[0](k)/k*norm            << " ";
    //fp << thetaT_ell_of_k_spline[2](k)            << " ";
    //fp << thetaT_ell_of_k[2] << " ";
    fp <<"\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output_test45(const std::string filename) const{
  //const double k_min_ = k_min;
  //const double k_max_ = k_max;
  //const int    n_pts =  10000;
  //double k = 340.*cosmo->get_H0()/Constants.c;
  const double Mpc = Constants.Mpc;
  //const double c_H0 = Constants.c/cosmo->get_H0();
  const double c_H0 = Constants.c/cosmo->get_H0();
  const double norm = 1e6*cosmo->get_H0()/Constants.c;
  arma::vec k_array_arma = arma::logspace(log10(k_min),log10(k_max),n_k);
  
  //Vector k_array;
  //////////////////////////////////////////////////////////////7for (int i=0;i<k_array_arma.size();i++){
  //////////////////////////////////////////////////////////////7    k_array.push_back(k_array_arma[i]);
  //////////////////////////////////////////////////////////////7}
  //Vector k_array = Utils::linspace(k_min,k_max,10000);

  //const double conversion_P = pow(cosmo->get_h(),3)/pow(Mpc,3);
  //const double conversion = Mpc/cosmo->get_h();
  //Vector k_array = Utils::linspace(k_min_, k_max_, n_pts);
  //Vector k_array_scaled;
  //for (int i=0;i<k_array.size();i++){
  //  k_array_scaled.push_back(k_array[i]*Mpc/cosmo->get_h());
  //}
   //= k_array*Mpc/cosmo->get_h();
  //double eta_0 = cosmo->eta_of_x(0.);
  int i=0;
  std::cout << "output" << std::endl;
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double k) {
    
    fp <<    k*c_H0                                     << " ";
    fp << pow(thetaT_ell_of_k[0][i],2)/k*norm            << " ";
    fp << pow(thetaT_ell_of_k_spline[0](k),2)/k*norm            << " ";
    //fp << pow(thetaT_ell_of_k[2][i],2)/k*norm            << " ";
    //fp << pow(thetaT_ell_of_k[3][i],2)/k*norm            << " ";
    //fp << pow(thetaT_ell_of_k[4][i],2)/k*norm            << " ";
    //p << thetaT_ell_of_k[0][i]/k*norm            << " ";
    //fp << abs(thetaT_ell_of_k[0][i])/k*norm            << " ";
    //fp << thetaT_ell_of_k_spline[2](k)            << " ";
    //fp << thetaT_ell_of_k[2] << " ";
    fp <<"\n";
    i += 1;
  };
  std::for_each(all_k_arrays[5].begin(), all_k_arrays[5].end(), print_data);
}


