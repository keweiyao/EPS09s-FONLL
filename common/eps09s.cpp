/*******************************************************************************
 * EPS09s.cpp
 *
 * Interface for spatially dependent nuclear modifications determined from
 * the EPS09 global analysis.
 *
 * When using this interface, please refer to the following publication:
 *
 * I. Helenius, K.J. Eskola, H. Honkanen and C.A. Salgado,
 *
 * "Impact-parameter dependent nuclear parton distribution functions: EPS09s 
 * and EKS98s and their applications in nuclear hard-processes"
 *
 * Eprint: arxiv:1205.5359 [hep-ph]
 * To appear in JHEP
 *
 * Details of the framework for determining the spatially dependent nPDFs can 
 * be found in the above paper.
 *
 * Questions & comments to:
 *  ilkka.helenius@jyu.fi
 *  kari.eskola@phys.jyu.fi
 *
 *******************************************************************************
 * Instructions:
 *
 * We provide two possible interfaces for the calculation of the spatially 
 * dependent nPDFs:
 * 
 * (1) This method is based on the Eq. (4.2) in our publication. The geometric
 *  integrals T_{AB}^{nm}(b_1,b_2) are separated from the kinematic parts for
 *  which the fit parameters c^i_n(x,Q^2) are provided by a function
 *
 *   eps09s_c(eps_order, eps_pset, a, x, q, cuv, cdv, cu, cd, cs, cc, cb, cg)
 * 
 *  where the inputs are
 * 
 *   eps_order: 1=LO, 2=NLO     ; int
 *   eps_pset : 1,...,31        ; int
 *              1     = central set               * See the EPS09 publication
 *              2,3   = error sets S{+1}, S{-1}   * (JHEP 0904 (2009) 065)
 *              4,5   = error sets S{+2}, S{-2}   * for the details of the
 *              ...   ...                         * error calculations
 *              30,31 = error sets {S+15}, {S-15} * 
 *   a        : Atomic number   ; int
 *   x        : Bjorken-x       ; double
 *   Q        : Scale in GeV    ; double
 *
 *  and the 4 fit parameters c_j^i  (j=1,2,3,4) for each flavor i are returned 
 *  in the following 4-dimensional tables (in units of [c_j^i] = (fm^2)^j) 
 *  (passed by their addresses)
 *
 *	cuv[4] = up valence
 *	cdv[4] = down valence
 *	cu[4]  = up sea
 *	cd[4]  = down sea
 *	cs[4]  = strange
 *	cc[4]  = charm
 *	cb[4]  = bottom
 *	cg[4]  = gluons
 *
 *  For the geometric integrals T_{AB}^{nm}(b_1,b_2) we provide the nuclear  
 *  thickness function T_A(s) as a function
 *
 *   ta_eps09s(a, s)
 *
 *  which returns the value of the thickness function as a double (in units of
 *  1/fm^2). The required inputs are:
 *
 *   a : the mass number of nucleus    ;int
 *   s : the distance from the center of the nucleus in transverse plane [fm]
 *       ;double
 *    
 *  For A > 2 The two parameter Woods-Saxon distribution is used
 *      A = 2 The deuterium wave function is used
 *      A = 1 return 1
 *
 *
 *  At the first call, the program calculates and saves the values of T_A(s_i) 
 *  for various fixed values of s_i. The value of T(s) at any given s is then 
 *  obtained by interpolation. 
 *  Thickness functions for at most 2 nuclei can appear at each time in the 
 *  memory. If more are called, the program will stop.
 *
 * (2) The nuclear modifications are precomputed using the fit parameters and 
 *  the thickness function. The modifications, already multiplied by the 
 *  thickness function T_A(s), are returned by a function
 * 
 *   eps09s(eps_order, eps_pset, a, x, q, s, ruv, rdv, ru, rd, rs, rc, rb, rg)
 *
 *  where the first inputs are the same as above and the s corresponds again
 *  to the distance from the center of nucleus in transverse plane [fm]. The 
 *  modifications for each flavor are returned in parameters ri (passed by 
 *  reference) with the naming scheme introduced in method (1) above. Using the
 *  notation introduced in the publication above we have
 *
 *    ri = T_A(s)*r_i(x,Q^2,s)
 *
 *  The dimensions are [ri] = 1/fm^2, [T_A(s)] = 1/fm^2 and [r_i(x,Q^2,s)] = 1.
 *
 *  As discussed in our paper, we recommend to use the method (1).
 * 
 *******************************************************************************
 *
 * The kinematic domain considered is the same as in the EPS09 analysis,
 *
 *             1e-6 <= x <= 1
 *              1.3 <= Q <= 1000 GeV.
 * 
 * Outside this domain, the values at the nearest end point is returned.
 *
 * The required tables for the interpolation are stored in files (62 in total)
 *
 *   eps09sORDERSET.dat,
 * 
 * where 
 *  ORDER=LO,NLO
 *  SET=1,...,31
 * 
 * The relevant files (LO/NLO) should be located in the working directory.
 *
 *******************************************************************************
 * TEST PROGRAM:
 *
 *  A test function is also included to this code. The function call is
 *
 *   test_eps09s()
 *
 *  which calculates some values of nuclear modifications and tests whether
 *  these are consistent with precomputed values. If all tests are succesfull,
 *  it is safe to use this program.
 *
 *******************************************************************************
 * FORTRAN USAGE:
 * 
 * 1) Uncomment the last few lines of this file (the block 'extern "C"{...}')
 *
 * 2) Compile and link (with gnu compilers):
 *
 *   1) Compile the .cpp code:
 *    $ g++ -c eps09s.cpp
 *
 *   2) Link the object file (eps09s.o) with your fortran program:
 *    $ gfortran FORTRAN_SOURCE.f eps09s.o -lstdc++
 *
 *  where FORTRAN_SOURCE.f is your FORTRAN program and -lstdc++ is the 
 *  C++ Standard Library which includes the required I/O functions etc.
 *
 *  To compile with your favorite compiler, just replace the g++ or gfortran
 *  with some other compiler.
 * 
 *  If you use f77 compiler, you might need to add an another underscore (_) to 
 *  the eps09s_c_ and ta_eps09s_ function calls at block 'extern "C"{...}'
 *  (eps09s_c_(...) -> eps09s_c__(...) and ta_eps09s_(...) -> ta_eps09s__(...) )
 * 
 * 3) From your FORTRAN code call the C++ subroutines and functions without the
 *    underscores in their names, e.g. call eps09s_c(...).
 *
 *******************************************************************************
 */

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace std;

namespace eps09_s{

/*
 * Class containing the fit parameters c_i for all flavors with an interpolation
 * method.
 *
 * Singleton design pattern is used to avoid global variables and not to read
 * the files over and over...
 */
class fit_parameters{
private:
             //set,flavor,nq,nx,nparam;
  double c_values[31][8][51][51][4];
  int order, pset;
  string get_filename(int eps_order, int eps_pset);
  static fit_parameters* _instance;
protected:
  fit_parameters(){}
  fit_parameters(int eps_order, int eps_pset);
  ~fit_parameters(){}
public:
  static fit_parameters* getInstance(int eps_order, int eps_pset);
  void interpolate_c_values(const double x, const double q, const int pset, 
			    double c[8][4]);
};

//Constructor
fit_parameters::fit_parameters(int eps_order, int eps_pset){
  order = eps_order, pset = eps_pset;
  for(int n = 0;n < 31;n++){
    string filename = get_filename(order,n+1);
    ifstream infile(filename.c_str());
    if(!(infile.good())){
      cout << "No file " << filename << " found!" << endl;
      cout << "Exiting program..." << endl;
      infile.close();
      exit(1);
    }
    string dummy;
    string file_version;
    infile >> dummy;
    infile >> file_version;
    if(file_version != "1.1"){ // Checks that file is valid
      cout << "File: " << filename << " not valid!" << endl;
      cout << "Exiting program..." << endl;
      exit(1);
    }

    //Read parameter values from file:
    for(unsigned int i = 0;i < 8;i++){ //Flavors    
      for(unsigned int j = 0;j < 51;j++){ //Q-grid
	for(unsigned int k = 0;k < 51;k++){ //x-grid
	  for(unsigned int l = 0; l < 4;l++){ //c_i
	    infile >> c_values[n][i][j][k][l];
	  }
	}
      }
    }
    infile.close();
  }
}

/*
 * Returns the static instance of fit parameters (singleton pattern).
 * Reads new files if different order is requested.
 */
fit_parameters* fit_parameters::_instance = 0;
fit_parameters* fit_parameters::getInstance(int eps_order, int eps_pset){
  if(_instance == 0)
    _instance = new fit_parameters(eps_order, eps_pset);
  if(_instance->order == eps_order){
    return _instance;    
  } else {
    _instance = new fit_parameters(eps_order, eps_pset);
    return _instance;
  }
}

/*
 * Constructs the file name for selected order & parameter set
 * using string stream from STD library
 *
 * "eps09s" + ORDER + PSET + ".dat"
 * ORDER = LO, NLO  PSET = 1,...,31
 */
string fit_parameters::get_filename(int eps_order, int eps_pset){
  string filename;
  switch(eps_order){
  case 1:
    filename = "EPS09s-data/eps09sLO";
    break;
  case 2:
    filename = "EPS09s-data/eps09sNLO";
    break;
  default: //Checks that order is valid
    cout << "Invalid EPS09s order!! order = " << eps_order << endl;
    cout << "Should be 1 or 2" << endl;
    cout << "Exiting program..." << endl;
    filename = " ";
    exit(2);
    break;
  }  
  if( (eps_pset < 1) || (eps_pset > 31) ){ //Checks that pset is valid
    cout << "Invalid EPS09s parameters set!! set = " << eps_pset << endl;
    cout << "Should be between 1 and 31" << endl;
    cout << "Exiting program..." << endl;
    exit(3);
  }
  stringstream ss_pset;
  ss_pset << eps_pset;
  filename += ss_pset.str();
  filename += ".dat";
  return filename;
}

/*
 * Polynomial interpolation with Newton's divided difference method.
 */
double pol_int(double f_i[], double x_i[], int n, double x){
  for(int i = 1;i < n;i++){
    for(int j = n-1;j > i - 1;j--){
      f_i[j] = (f_i[j] - f_i[j-1])/(x_i[j] - x_i[j-i]);
    }
  }

  double p = f_i[n-1];
  for(int i = n-2;i > -1;i--){
    p = (x - x_i[i])*p + f_i[i];
  }
  
  return p;
}

/*
 * Interpolation for fit parameters for given x, q for all flavors
 */
void fit_parameters::interpolate_c_values(const double xx, const double q, 
					  const int pset, double c[8][4]){
  double xmin = 1e-6, xmax = 1, q2min = 1.69, q2max = 1000000;
  double xcut = 0.1;
  int nx = 51, nq = 51; 
  int nxlog = 25, nxlin = 25;
  int nparams = 4;

  //Freeze Q^2 if outside the limits
  double q2 = q*q;
  if(q2 > q2max){
    q2 = q2max;
  } else if (q2 < q2min){
    q2 = q2min;
  }

  //Freeze x values if outside the limits
  double x = xx;
  if(x > xmax){
    x = xmax;
  } else if (x < xmin){
    x = xmin;
  }

  //Calculate the position in log(log Q^2) grid:
  double i_qd = (nq-1)*log( log(q2)/log(q2min) )/log( log(q2max)/log(q2min) );
  int i_q = static_cast<int>( i_qd );

  //Set the q-index to interval [1,...,49]
  if(i_q < 1){
    i_q = 1;
  } else if (i_q > (nq-2) ){
    i_q = nq-2;
  }

  //Calculate the three nearest points in log(log Q^2) grid
  double q_i[3];
  for(int i = 0;i < 3;i++){
    q_i[i] = sqrt( exp( pow(log(q2max),(i_q+i-1)*1.0/(nq-1))*
			pow(log(q2min),1-(i_q+i-1)*1.0/(nq-1)) ) );
  }

  //Calculate the position in log(x) or x grid
  int i_x;
  if(x <= xcut){
    i_x = static_cast<int>( nxlog*log(x/xmin)/log(xcut/xmin) );
  } else {
    i_x = static_cast<int>( (x-xcut)*nxlin/(xmax-xcut) + nxlog );
  }

  //Set the x-index to interval [1,...,48]
  if(i_x < 1){
    i_x = 1;
  } else if (i_x > nx - 3){
    i_x = nx - 3;
  }

  //Calculate the three nearest points in log(x) or x grid
  double x_i[4];
  for(int i = 0;i < 4;i++){
    if(i_x-1+i < nxlog){
      x_i[i] = xmin*exp( ((i_x-1+i)*1.0/nxlog)*log(xcut/xmin) );
    } else {
      x_i[i] = ( ( i_x - 1 + i - nxlog)*1.0/nxlin )*(xmax-xcut) + xcut;
    }
  }

  //Interpolate c parameters:
  double ccc[8][nparams][3][4], cc[8][nparams][3];
  for(int i = 0; i < 8;i++){
    for(int j = 0; j < nparams;j++){
      for(int k = 0; k < 3;k++){
		for(int l = 0; l < 4;l++){
	 	 ccc[i][j][k][l] = c_values[pset-1][i][i_q+k-1][i_x+l-1][j];
		}
		cc[i][j][k] = pol_int(ccc[i][j][k], x_i, 4, x);
      }
      c[i][j] = pol_int(cc[i][j], q_i, 3, sqrt(q2));
    }
  }
  return;
}


/*
 * Recursive function for adaptiveBoole function. Splits the interval so that
 * the required accuracy is obtained
 */ 
double adaptiveBooleRec(double (*f)(double, void*), void *p,
			double x1, double x9, double epsilon, double S,
			double f1, double f3, double f5, double f7, 
			double f9, int rec_steps){
  double h = x9 - x1, x5 = x1 + h/2;
  double x2 = x1 + h/8, x4 = x1 + 3*h/8, x6 = x1 + 5*h/8, x8 = x1 + 7*h/8;
  double f2 = f(x2,p), f4 = f(x4,p), f6 = f(x6,p), f8 = f(x8,p);
  double S_left = (h/180)*(7*f1 + 32*f2 + 12*f3 + 32*f4 + 7*f5);
  double S_right = (h/180)*(7*f5 + 32*f6 + 12*f7 + 32*f8 + 7*f9);
  double S_tot = S_left + S_right;
  if (rec_steps <= 0 || fabs(S_tot - S) <= 15*epsilon){
    return S_tot + (S_tot - S)/15;
  }
  return adaptiveBooleRec(f,p,x1,x5,epsilon/2, S_left, f1, f2, f3, f4, f5,
			  rec_steps-1)
    +    adaptiveBooleRec(f,p,x5,x9,epsilon/2, S_right, f5, f6, f7, f8, f9,
			  rec_steps-1);
}         
 
/*
 * Adaptive Boole's Rule
 *
 * Input:
 *  double (*f)(double,void*) = pointer to function with unknown parameters
 *  double x1, x5             = interval [a,b]
 *  double epsilon            = error tolerance
 *  int max_recursion_steps     = maximal number of recursions
 *
 * Output:
 *  the result from the integration (double)
 *
 * Parameters tested to be suitable for this purpose. 
 * For other purposes use with caution (or not at all)!
 */ 
double adaptiveBoole(double (*f)(double, void* ), void *p, 
	       double x1, double x5, double epsilon, int max_recursion_steps){
  double h = x5 - x1, x2 = x1 + h/4, x3 = x1 + h/2, x4 = x1 + 3*h/4;
  double f1 = f(x1,p), f2 = f(x2,p), f3 = f(x3,p), f4 = f(x4,p), f5 = f(x5,p);
  double S = (h/90)*(7*f1 + 32*f2 + 12*f3 + 32*f4 + 7*f5);
  return adaptiveBooleRec(f,p,x1,x5,epsilon, S, f1, f2, f3, f4, f5,
			  max_recursion_steps);
}                   

struct paramsi1d1 { int i1; double d1; };

/*
 * Calculates the Woods-Saxon distribution for given position.
 *
 * Mapped with change of variables as z = (1-t)/t so that the integral from
 * 0 < z < Infinity becomes an integral where 0 < t < 1.
 * Woods-Saxon parameters are 
 *
 *  R_A = 1.12*A^(1/3) - 0.86*A(-1/3) fm
 *    d = 0.54 fm
 *  n_0 = 3*A/(4*Pi*R_A^3)*1/(1 + (Pi*d/R_A)^2) (Correct normalization for A>3)
 */
double mapped_woodSaxon_density(double t, void *params){
  if(t == 0) return 0;
  paramsi1d1 * wsparams = (paramsi1d1*)params;
  int a = wsparams->i1;
  double s = wsparams->d1;
  double d = 0.54;
  double r = 1.12*pow(a,1./3.) - 0.86*pow(a,-1./3.);
  double n0 = 0.75*a/(acos(-1.0)*pow(r,3))/(1+(pow(acos(-1.0)*d,2))/(r*r));
  return 2*n0/( t*t*(1 + exp( (sqrt( s*s + (1-t)*(1-t)/(t*t) ) - r )/d ) ) );
}

/*
 * Calculates T_A^{WS}(s) in given point s for given A
 */
double taWoodsSaxon(int a, double s){
  paramsi1d1 params = { a, s };
  return adaptiveBoole(mapped_woodSaxon_density, &params, 0,1,1e-12, 12);
}

  /* Deuterium thickness function */

//S-wave 
double deuterium_u(double r, double beta, double gamma, double epsilon, 
		   double xc, double alpha, double n){
  if(alpha*r > xc){
    return n*sqrt(1 - epsilon*epsilon)*( 1 - exp(-beta*(alpha*r-xc) ) )
      *exp(-alpha*r)/r;
  } else {
    return 0;
  }
}
//D-wave
double deuterium_w(double r, double beta, double gamma, double epsilon, 
		  double xc, double alpha, double n){
  if(alpha*r > xc){
    return n*epsilon*( 1 - exp( -gamma*(alpha*r-xc) ) )*
      ( 1 - exp( -gamma*(alpha*r-xc) ) )*exp(-alpha*r)*
      ( 1 + 3*(( 1 - exp(-gamma*alpha*r) )/(alpha*r))*
	( 1 + ( 1 - exp(-gamma*alpha*r) )/(alpha*r) ) )/r;
  } else {
    return 0;
  }
}

/*
 * Calculates the square of the wave function (two possible parameter set from
 * Nucl.Phys. A730 (2004) 448-459)
 */
double deuterium_psi2_ds(double r, int d_set){
  double alpha = 1.0/4.316;
  double beta, gamma, epsilon, xc, rho, n;
  switch(d_set){
  case 1:
    beta = 4.680, gamma = 2.494, epsilon = 0.03232, xc = 0;
    rho = -27.944041219;
    n = sqrt( 2*alpha/(1 - alpha*rho) );
    break;
  case 2:
    beta = 9.045, gamma = 4.799, epsilon = 0.02438, xc = 0.13;
    rho = -28.1136;
    n = sqrt( 2*alpha/(1 - alpha*rho) );
    break;
  default:
    cout << "Wrong deuteron wave function parameter set! " << endl;
    exit(4);
    break;
  }
  double w = deuterium_w(r, beta, gamma, epsilon, xc, alpha, n);
  double u = deuterium_u(r, beta, gamma, epsilon, xc, alpha, n);
  return (u*u + w*w);
}

/*
 * Converts the wave function from p-n distance to distance from the CM
 */
double deuterium_psi2_ds_cm(double r, int d_set){
  return 8*deuterium_psi2_ds(2*r, d_set);
}

/*
 * Integrand for t_deuterium
 */
double mapped_t_d_integrand(double rl, void *data){
  if(rl == 0) return 0; 
  paramsi1d1 * parameters = (paramsi1d1*)data;
  int paramset = parameters->i1;
  double rt = parameters->d1;
  return 2*deuterium_psi2_ds_cm(sqrt(rt*rt + (1-rl)*(1-rl)/(rl*rl)), 
				paramset)/(rl*rl);
}

/*
 * Deuterium thickness function from Hulthen wave function
 */
double t_deuterium(double rT){
  paramsi1d1 params = { 1, rT };
  return 2*adaptiveBoole(mapped_t_d_integrand, &params, 0,1,1e-12, 12);
}

/*
 * Class containing the calculated values of T_A(s) for different s values
 * and a method which interpolates the T_A(s) for given s.
 */
class thickness_function{
private:
  int n_points_lins, n_points_linu;
  double ta[200];
  double tail_length;
  int a;
  static thickness_function* _instancea;
  static thickness_function* _instanceb;
protected:
  thickness_function(){}
  thickness_function(int a);
  ~thickness_function(){}
public:
  static thickness_function* getInstance(int a);
  double operator()(double s);
};

/*
 * Constructor:
 *
 * Calculates thickness function values for given A for several s values
 */
thickness_function::thickness_function(int aa){
  a = aa;
  n_points_lins = 150;
  n_points_linu = 50;
  tail_length = 4; //[fm]
  double r = 1.12*pow(a,1./3.) - 0.86*pow(a,-1./3.);
  double s_cut = r + tail_length;
  
  if(a == 2){ //Deuterium wave function
    for(int i = 0; i < n_points_lins; i++){
      double s = i*(s_cut)/(1.0*n_points_lins);
      ta[i] = t_deuterium(s);
    }
    for(int i = 0; i < n_points_linu; i++){
      double s = s_cut + s_cut/n_points_lins*i*i;
      ta[i+n_points_lins] = t_deuterium(s);
    }
  } else if (a > 2){ //Woods-Saxon
    for(int i = 0; i < n_points_lins; i++){
      double s = i*(s_cut)/(1.0*n_points_lins);
      ta[i] = taWoodsSaxon( a, s );
    }    
    for(int i = 0; i < n_points_linu; i++){
      double s = s_cut + s_cut/n_points_lins*i*i;
      ta[i+n_points_lins] = taWoodsSaxon( a, s );
    }
  } else { //No thickness for A < 2, Program stops!
    cout << "Error: no thickness function for A < 2!" << endl;
    exit(5);
  }
}

/*
 * Thickness function interpolation:
 *  Interpolates with linear interval in s from 0 to s = R_A + 4fm and 
 *  from this on with linear interval in u = 1/(1 + s).
 */
double thickness_function::operator()(double s){  
  s = fabs(s);
  int ns = this->n_points_lins;
  int nu = this->n_points_linu;
  double aa = this->a;
  double tl = this->tail_length;
  double r = 1.12*pow(aa,1./3.) - 0.86*pow(aa,-1./3.);
  double s_cut = r + tl;
  int index1, index2;
  double s1, s2, interpolated_ta;
  if( s < r + tl ){
    index1 = min(static_cast<int>(floor( ns*s/s_cut )), ns - 1);
    index2 = index1 + 1;

    s1 = index1*s_cut/ns;
    s2 = index2*s_cut/ns;

    if(index1 != index2){
      interpolated_ta = ta[index1]*exp( log(ta[index2]/ta[index1])
      				*(s-s1)/(s2-s1) );
    } else {
      interpolated_ta = ta[index1];
    }

  } else {
    index1 = min(static_cast<int>(floor( sqrt( (ns/s_cut)*(s - s_cut) ) )),
		 nu - 1);
    index2 = min(index1 + 1, nu - 1 );

    s1 = s_cut + s_cut*index1*index1/ns;
    s2 = s_cut + s_cut*index2*index2/ns;

    if(index1 == index2){
      interpolated_ta = 0;
    } else {
      interpolated_ta = ta[index1 + ns]*exp( log(ta[index2+ns]/ta[index1+ns])
					     *(s-s1)/(s2-s1) );
    }
  }  

  return interpolated_ta;
}

/*
 * Return an instance of static thickness function object (singleton pattern)
 * Can handle thickness functions for two different nucleus.
 */
thickness_function* thickness_function::_instancea = 0;
thickness_function* thickness_function::_instanceb = 0;
thickness_function* thickness_function::getInstance(int aa){
  if(_instancea == 0){
    _instancea = new thickness_function(aa);
  }
  if(_instancea->a == aa){
    return _instancea;
  } else {
    if(_instanceb == 0){
      _instanceb = new thickness_function(aa);
    }
    if(_instanceb->a == aa){
      return _instanceb;
    } else {
      cout << "Support only for 2 nucleus in a run! Exiting program..." << endl;
      exit(6);
    }
  }
}

} //End of eps09_s namespace

/*
 * Calculates the impact parameter dependent nuclear modification factors
 * for all flavors for given x and Q values, nucleus A, and for spatial
 * distance s.
 *
 * Modifications are calculated as
 * 
 * r_i = T_A(s)*[1 + c_1(x,Q)*T_A(s) + c_2(x,Q)*[T_A(s)]^2 
 *                 + c_3(x,Q)*[T_A(s)]^3 + c_4(x,Q)*[T_A(s)]^4] 
 */
void eps09s(const int eps_order, const int eps_pset, const int a, 
	    const double x, const double q, const double s, double &ruv, 
	    double &rdv, double &ru, double &rd, double &rs, double &rc, 
	    double &rb, double &rg){
  double r[8];
  if(a > 2){
    if(a < 16){
      cout << "EPS09s Warning:" << endl;
      cout << "Only nuclei with A >= 16 used for fitting!"
	   << endl;
    }
    eps09_s::thickness_function *ta;
    ta = eps09_s::thickness_function::getInstance(a);
    double tas = (*ta)(s);
    eps09_s::fit_parameters *fit_params;
    fit_params = eps09_s::fit_parameters::getInstance(eps_order, eps_pset);
    double c[8][4];
    fit_params->interpolate_c_values( x, q, eps_pset, c);
    for(int i = 0;i < 8;i++){
      r[i] = tas;
      for (int j = 0;j < 4;j++){
	r[i] += c[i][j]*pow(tas,j+2);
      }
    }
  } else if( a == 2 ){
    eps09_s::thickness_function *ta;
    ta = eps09_s::thickness_function::getInstance(a);
    double tas = (*ta)(s);
    for(int i = 0;i < 8;i++){
      r[i] = tas;
    }
  } else {
    for(int i = 0;i < 8;i++){
      r[i] = 1;
    }
  }
  ruv = r[0], rdv = r[1], ru = r[2], rd = r[3], rs = r[4], rc = r[5], 
    rb = r[6], rg = r[7];
}

/*
 * Returns the fit parameter values for given x and Q value. If A < 3
 * return only zeros as nuclear modifications are assumed to be negligible.
 *
 * This should be used with ta_eps09s(a,s) function so that nuclear modification
 * can be calculated as
 *
 *  1 + c^i_1(x,Q^2)*T_A(s) + c^i_2(x,Q^2)*[T_A(s)]^2 + ... 
 * 
 * The fit parameters are tabulated as
 * 
 *      c_1:   c_2:   c_3:   c_4:
 *  uv: cuv[0] cuv[1] cuv[2] cuv[3] 
 *  dv: cdv[0] cdv[1] cdv[2] cdv[3]
 *  us: cus[0] cus[1] cus[2] cus[3]
 *  ds: cds[0] cds[1] cds[2] cds[3]
 *  s:   cs[0]  cs[1]  cs[2]  cs[3]
 *  c:   cc[0]  cc[1]  cc[2]  cc[3]
 *  b:   cb[0]  cb[1]  cb[2]  cb[3]
 *  g:   cg[0]  cg[1]  cg[2]  cg[3]
 *
 */
void eps09s_c(const int eps_order, const int eps_pset, const int a, 
	      const double x, const double q, double cuv[4], double cdv[4], 
	      double cu[4], double cd[4], double cs[4], double cc[4], 
	      double cb[4], double cg[4]){
  if(a > 2){
    if(a < 16){
      cout << "EPS09s Warning:" << endl;
      cout << "Only nuclei with A >= 16 used for fitting!"
	   << endl;
    }
    eps09_s::fit_parameters *fit_params;
    fit_params = eps09_s::fit_parameters::getInstance(eps_order, eps_pset);
    double c[8][4];
    fit_params->interpolate_c_values( x, q, eps_pset, c);
    for(int i = 0;i < 4;i++){
      cuv[i] = c[0][i], cdv[i] = c[1][i], cu[i] = c[2][i], cd[i] = c[3][i], 
	cs[i] = c[4][i], cc[i] = c[5][i], cb[i] = c[6][i], cg[i] = c[7][i];
    }
  } else {
    for(int i = 0;i < 4;i++){
      cuv[i] = 0, cdv[i] = 0, cu[i] = 0, cd[i] = 0, cs[i] = 0, cc[i] = 0,
	cb[i] = 0, cg[i] = 0;
    }
  }
  return;
}

/*
 * Return the thickness function value T_A(s) for given A and s. This should be 
 * used with eps09s_c so that the spatial d^2 s integration can be done 
 * separately. This can speed up the computation significantly.
 */
double ta_eps09s(const int a, const double s){
  if(a > 1){ 
    eps09_s::thickness_function *ta;
    ta = eps09_s::thickness_function::getInstance(a);
    return (*ta)(s);
  } else {
    return 1;
  }
}

/*
 * Test program for the interpolation codes
 */
void test_eps09s(){
  int order = 1;
  int set = 1;
  int a = 208;
  double x = 0.02;
  double q = 1.4;
  double s = 5;
  double r[8],c[8][4];

  eps09s(order, set, a, x, q, s, r[0], r[1], r[2], r[3], r[4], r[5], r[6],r[7]);
  eps09s_c(order, set, a, x, q, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]);
  double ta = ta_eps09s(a,s);

  double ruv = 1, rus = 1, rg = 1;
  for(int i = 0; i < 4; i++){
    ruv += c[0][i]*pow(ta,i+1);
    rus += c[2][i]*pow(ta,i+1);
    rg += c[7][i]*pow(ta,i+1);
  }

  bool test_ok = true;
  if( !(fabs(r[0]/ta - 0.981036) < 1e-5) ) test_ok = false;
  if( !(fabs(r[2]/ta - 0.779357) < 1e-5) ) test_ok = false;
  if( !(fabs(r[7]/ta - 0.80437) < 1e-5) ) test_ok = false;
  if( !(fabs(ruv - 0.981036) < 1e-5) ) test_ok = false;
  if( !(fabs(rus - 0.779357) < 1e-5) ) test_ok = false;
  if( !(fabs(rg - 0.80437) < 1e-5) ) test_ok = false;

  if(test_ok){
    cout << "All tests were succesful!" << endl;
  } else {
    cout << "Test failed!" << endl;
  }
  
}


//-----Generate shaodowing tables with different s1

/*int main(){

	test_eps09s();
	return 0;
}*/

//For FORTRAN usage: Uncomment the following lines (remove /* and */)


extern "C" {

void eps09s_(int *eps_order, int *eps_pset, int *a, double *x, double *q, 
	     double *s, double *ruv, double *rdv, double *ru, double *rd, 
	     double *rs, double *rc, double *rb, double *rg){
  eps09s( *eps_order, *eps_pset, *a, *x, *q, *s, *ruv, *rdv, *ru, *rd, *rs, 
	  *rc, *rb, *rg);
}

void eps09s_c_(int *eps_order, int *eps_pset, int *a, double *x, double *q, 
	       double *cuv, double *cdv, double *cu, double *cd, 
	       double *cs, double *cc, double *cb, double *cg){
  eps09s_c( *eps_order, *eps_pset, *a, *x, *q, cuv, cdv, cu, cd, cs, 
	    cc, cb, cg);
}

double ta_eps09s_(int *a, double *s){
  return ta_eps09s( *a, *s);
}

void test_eps09s_(){
  test_eps09s();
}

}

