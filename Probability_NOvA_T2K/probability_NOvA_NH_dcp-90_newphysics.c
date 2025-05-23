#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "3x3-C/zhetrd3.h"
#include "3x3-C/zheevq3.h"
#include <Python.h>
#include <globes/globes.h> /* GLoBES library */

// Macros
#define SQR(x) ((x) * (x))                         // x^2
#define SQR_ABS(x) (SQR(creal(x)) + SQR(cimag(x))) // |x|^2

/* New parameters */

#define GLB_M0 6
#define GLB_MU 7   
#define GLB_N 8                  


/***************************************************************************
 *     U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E       *
 ***************************************************************************/

double th12;
double th13;
double th23;


double dcp;

double m0;
double mu;
double N;


double d21;
double d31;



// Initialize Python once globally
void initialize_python_once() {
    static int initialized = 0;
    if (!initialized) {
        Py_Initialize();
        PyRun_SimpleString("import sys; sys.path.insert(0, '/Users/ushak/Documents/CERN/Manuel')");
        PyRun_SimpleString("import site; site.addsitedir('/Users/ushak/venvs/globes/lib/python3.13/site-packages')");

        //PyRun_SimpleString("import sys; sys.path.append('.')");
        initialized = 1;
    }
}

void finalize_python_once() {
    static int finalized = 0;
    if (!finalized) {
        Py_Finalize();
        finalized = 1;
    }
}

// Diagonalize using Python
void diagonalize_in_python(double complex H[6][6], double evals[6], double complex evecs[6][6]) {
    initialize_python_once();

    PyObject *pName = PyUnicode_DecodeFSDefault("diag");
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (!pModule) {
        PyErr_Print();
        fprintf(stderr, "Failed to load Python module\n");
        return;
    }

    PyObject *pFunc = PyObject_GetAttrString(pModule, "diagonalize");
    if (!pFunc || !PyCallable_Check(pFunc)) {
        PyErr_Print();
        fprintf(stderr, "Cannot find function \"diagonalize\"\n");
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
        return;
    }

    PyObject *pList = PyList_New(36);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            int idx = i * 6 + j;
            PyObject *complex_tuple = PyComplex_FromDoubles(creal(H[i][j]), cimag(H[i][j]));
            PyList_SetItem(pList, idx, complex_tuple);
        }
    }

    PyObject *pArgs = PyTuple_Pack(1, pList);
    PyObject *pResult = PyObject_CallObject(pFunc, pArgs);
    Py_DECREF(pArgs);
    Py_DECREF(pList);

    if (!pResult) {
        PyErr_Print();
        fprintf(stderr, "Call to diagonalize failed\n");
        Py_DECREF(pFunc);
        Py_DECREF(pModule);
        return;
    }

    PyObject *pEvals = PyTuple_GetItem(pResult, 0);
    PyObject *pEvecs = PyTuple_GetItem(pResult, 1);

    for (int i = 0; i < 6; i++) {
        evals[i] = PyFloat_AsDouble(PyList_GetItem(pEvals, i));
    }

    for (int i = 0; i < 36; i++) {
        PyObject *z = PyList_GetItem(pEvecs, i);
        double re = PyComplex_RealAsDouble(z);
        double im = PyComplex_ImagAsDouble(z);
        evecs[i % 6][i / 6] = re + I * im;
    }

    Py_DECREF(pResult);
    Py_DECREF(pFunc);
    Py_DECREF(pModule);
}
/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/
int my_set_oscillation_parameters(glb_params p, void *user_data)
{
  th12 = glbGetOscParams(p, GLB_THETA_12);
  th13 = glbGetOscParams(p, GLB_THETA_13);
  th23 = glbGetOscParams(p, GLB_THETA_23);
  dcp = glbGetOscParams(p, GLB_DELTA_CP);
  d21 = glbGetOscParams(p, GLB_DM_21); /* Convert to GeV^2 */
  d31 = glbGetOscParams(p, GLB_DM_31); /* Convert to GeV^2 */

  m0 = glbGetOscParams(p, GLB_M0);
  mu = glbGetOscParams(p, GLB_MU);
  N = glbGetOscParams(p, GLB_N);
  

  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/
int my_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p, th12, GLB_THETA_12);
  glbSetOscParams(p, th13, GLB_THETA_13);
  glbSetOscParams(p, th23, GLB_THETA_23);
  glbSetOscParams(p, dcp, GLB_DELTA_CP);
  glbSetOscParams(p, d21, GLB_DM_21); 
  glbSetOscParams(p, d31, GLB_DM_31); 

  glbSetOscParams(p, m0, GLB_M0);
  glbSetOscParams(p, mu, GLB_MU);
  glbSetOscParams(p, N, GLB_N);
  
  

  return 0;
}



/***************************************************************************
 * Calculate oscillation probabilities.                                    *
 * Since for our setup, only P_ee is required, all other entries of P are  *
 * set to zero for simplicity. Furthermore, we neglect matter effects and  *
 * the filter feature (parameter filter_sigma).                            *
 * The formula for P_ee is Eq. (36) from hep-ph/0502147.                   *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:            The buffer where the probabilities are to be stored     *
 *   cp_sign:      +1 if probalities for neutrinos are requested, -1 for   *
 *                 anti-neutrinos.                                         *
 *   E:            The neutrino energy in GeV                              *
 *   psteps:       Number of constant density layers in the matter profile *
 *   length:       The lengths of these layers in km                       *
 *   density:      The individual densities of these layers in g/cm^3      *
 *   filter_sigma: Width of low-pass filter as given in the AEDL file      *
 ***************************************************************************/

int my_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

  double complex A[3][3];
  double complex Q[6][6];
  double w[6];
  double rho, Ve;
  double sin12, sin13, sin23, cos12, cos13, cos23;
  double sinb12, sinb13, sinb23, cosb12, cosb13, cosb23;

  int z, L, ii, jj;
  /* Set all probabilities to zero initially */
  for (ii = 0; ii < 3; ii++)
  {
    for (jj = 0; jj < 3; jj++)
    {
      P[ii][jj] = 0.0;
    }
  }

  rho = 2.8;

  /* Calculate total baseline */
  L = 0.0;

  double complex phase[6][psteps];
    double complex Smat[6][6][psteps];

  for (z = 0; z < psteps; z++)
  {
    L += length[z];
    // L = GLB_KM_TO_EV(L) * 1.0e9;      /* Convert to GeV^{-1} */

    // double L = 810 * 1.0e9;

    sin12 = sin(th12);
    sin13 = sin(th13);
    sin23 = sin(th23);
    cos12 = cos(th12);
    cos13 = cos(th13);
    cos23 = cos(th23);
/*Define PMNS matrix elements*/
    double complex Ue1 = cos12 * cos13;
    double complex Ue2 = sin12 * cos13;
    double complex Ue3 = sin13 * cexp(-I * cp_sign * dcp);

    double complex Um1 = -sin12 * cos23 - cos12 * sin23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Um2 = cos12 * cos23 - sin12 * sin23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Um3 = sin23 * cos13;

    double complex Ut1 = sin12 * sin23 - cos12 * cos23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Ut2 = -cos12 * sin23 - sin12 * cos23 * sin13 * cexp(I * cp_sign * dcp);
    double complex Ut3 = cos23 * cos13;

// Extended 6x6 mixing matrix
complex double U[6][6] = {
    {sqrt((N-1)/N) * Ue1,sqrt((N-1)/N) * Ue2,sqrt((N-1)/N) * Ue3, sqrt(1/N) * Ue1, sqrt(1/N) * Ue2, sqrt(1/N) * Ue3},     // electron row
    {sqrt((N-1)/N) * Um1,sqrt((N-1)/N) * Um2,sqrt((N-1)/N) * Um3, sqrt(1/N) * Um1, sqrt(1/N) * Um2, sqrt(1/N) * Um3},     // muon row
    {sqrt((N-1)/N) * Ut1,sqrt((N-1)/N) * Ut2,sqrt((N-1)/N) * Ut3, sqrt(1/N) * Ut1, sqrt(1/N) * Ut2, sqrt(1/N) * Ut3},     // tau row
    {-sqrt(1/N) * Ue1,-sqrt(1/N) * Ue2, -sqrt(1/N) * Ue3, sqrt((N-1)/N) * Ue1, sqrt((N-1)/N) * Ue2, sqrt((N-1)/N) * Ue3},     // new sterile 1
    {-sqrt(1/N) * Um1,-sqrt(1/N) * Um2, -sqrt(1/N) * Um3, sqrt((N-1)/N) * Um1, sqrt((N-1)/N) * Um2, sqrt((N-1)/N) * Um3},     // new sterile 2
    {-sqrt(1/N) * Ut1,-sqrt(1/N) * Ut2, -sqrt(1/N) * Ut3, sqrt((N-1)/N) * Ut1, sqrt((N-1)/N) * Ut2, sqrt((N-1)/N) * Ut3}      // new sterile 3
};


    Ve = cp_sign * 0.5 * 0.000076 * density[z]; // sqrt(2) * G_F * N_e

   double m1=m0, m2=sqrt(m1 *m1 + d21), m3=sqrt (m1 * m1 + d31), m1H=mu * m1, m2H=mu * m2, m3H=mu * m3;
// Extended 6x6 mass squared matrix
complex double M2[6][6] = {
    {m1 * m1, 0, 0, 0.0, 0.0, 0.0},     // electron row
    {0, m2 * m2, 0, 0.0, 0.0, 0.0},     // muon row
    {0, 0, m3 * m3, 0.0, 0.0, 0.0},     // tau row
    {0.0, 0.0, 0.0, m1H * m1H, 0.0, 0.0},     // new sterile 1
    {0.0, 0.0, 0.0, 0.0, m2H * m2H, 0.0},     // new sterile 2
    {0.0, 0.0, 0.0, 0.0, 0.0, m3H * m3H}      // new sterile 3
};

// Extended 6x6 matter potential matrix
complex double V[6][6] = {
    {Ve, 0, 0, 0.0, 0.0, 0.0},     // electron row
    {0, 0, 0, 0.0, 0.0, 0.0},     // muon row
    {0, 0, 0, 0.0, 0.0, 0.0},     // tau row
    {0.0, 0.0, 0.0, 0, 0.0, 0.0},     // new sterile 1
    {0.0, 0.0, 0.0, 0.0, 0, 0.0},     // new sterile 2
    {0.0, 0.0, 0.0, 0.0, 0.0, 0}      // new sterile 3
};
double complex H[6][6];

// Compute the Hamiltonian H = 1/(2E) * U * M2 * U^dagger + V
for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
        H[i][j] = 0.0 + 0.0 * I;
        for (int k = 0; k < 6; k++) {
            for (int l = 0; l < 6; l++) {
                H[i][j] += U[i][k] * M2[k][l] * conj(U[j][l]);
            }
        }
        H[i][j] /= (2.0 * E);
        H[i][j] += V[i][j];  // Add matter potential
    }
}
diagonalize_in_python(H, w, Q);  // This will call Python and fill the outputs

    // double complex eps00 = cp_sign * eps_e_e * 1e18;
    // double complex eps01 = cp_sign * eps_e_m * cexp(I * cp_sign * phi_e_m) * 1e18;
    // double complex eps02 = cp_sign * eps_e_t * cexp(I * cp_sign * phi_e_t) * 1e18;

    // double complex eps10 = conj(eps01);
    // b double complex eps11 = 0;
    // double complex eps12 = 0;

    // double complex eps20 = conj(eps02);
    // double complex eps21 = 0;
    // double complex eps22 = 0;

    // fprintf(stdout, "%g %g %g\n", cp_sign*0.5*0.000076*rho, cp_sign*0.5*0.000076*density[z], length[z] );

    

    // Time evolution operator: S = Q × diag(exp(-i λ_j L)) × Q^†
    for (int j = 0; j < 6; j++) {
        phase[j][z] = cexp(-I * w[j] * length[z] * 4 * 1.27);  // 1/GeV to km
    }

    // Smat = Q * diag(phase) * Q†
    for (int alpha = 0; alpha < 6; alpha++) {
        for (int beta = 0; beta < 6; beta++) {
            Smat[alpha][beta][z] = 0.0 + 0.0 * I;
            for (int j = 0; j < 6; j++) {
                Smat[alpha][beta][z] += Q[alpha][j] * phase[j][z] * conj(Q[beta][j]);
            }
        }
    }

  }

  for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            P[a][b] = pow(cabs(Smat[b][a][0]), 2);
        }
    }


  return (0);
}



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/***************************************************************************
 *                            M A I N   P R O G R A M                      *
 ***************************************************************************/

int main(int argc, char *argv[])
{

//   char MYFILE[] = "dchi2_emu_IH_5+5_mu+e.dat"; 
  char MYFILE[] = "Prob_NOvA_NH_dcp-90_matter_newphysics.dat"; 
  // char MYFILE[] = "dchi2_emu_IH_5+5+1+1_mu+tau.dat"; 
  // char MYFILE[] = "dchi2_emu_IH_5+5_mu+e+tau.dat"; 
  // char MYFILE[] = "dchi2_emu_IH_5+5+1+1_mu+e+tau.dat"; 

  /* Initialize libglobes */
  glbInit(argv[0]);
  glbRegisterProbabilityEngine(9, /* Number of parameters */
                               &my_probability_matrix,
                               &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters,
                               NULL);



  glbInitExperiment("2024_nova_app.glb", &glb_experiment_list[0], &glb_num_of_exps);



  /* Intitialize output */
  FILE *outfile = NULL;

  outfile = fopen(MYFILE, "w");
  if (outfile == NULL)
  {
    printf("Error opening output file.\n");
    return -1;
  }


  //double true_deltacp = 0.82*M_PI;

  double true_dcp = -90 * M_PI/180 ;

  // IH true, NH test
  double true_theta12 = 33.82*M_PI/180;
  double true_theta13 = 8.60*M_PI/180;
  double true_theta23 = 48.6*M_PI/180;
  double true_sdm = 7.39e-5; // dm21
  double true_ldm = 2.528e-3; // dm32
  //double true_ldm = true_sdm + true_dm32; // dm31

  double m0= 0.01;
  double N = 10;
  double mu = 5;
 

  glb_params true_values = glbAllocParams();
  //glb_params test_values = glbAllocParams();

  // True values
  glbSetOscParams(true_values, true_theta12, GLB_THETA_12);
  glbSetOscParams(true_values, true_theta13, GLB_THETA_13);
  glbSetOscParams(true_values, true_theta23, GLB_THETA_23);
  glbSetOscParams(true_values, true_sdm, GLB_DM_21);
  glbSetOscParams(true_values, true_ldm, GLB_DM_31);
  
  glbSetOscParams(true_values, true_dcp, GLB_DELTA_CP);
  
  glbSetOscParams(true_values, m0, GLB_M0);
  glbSetOscParams(true_values, N, GLB_N);
  glbSetOscParams(true_values, mu, GLB_MU);
  

  

  glbSetDensityParams(true_values, 1.0, GLB_ALL);
  //glbSetDensityParams(test_values, 1.0, GLB_ALL);


  

    glbSetOscillationParameters(true_values);
    glbSetRates();
double e, p, q, r, s;
  for(e=0.01;e<=5;e+=.01)
  {
   
 
    
    
p=glbProfileProbability(0,2,1, +1, e);

q=glbProfileProbability(0,2,1, -1, e);

fprintf(outfile, "%g %g %g \n", e,p, q);

fprintf(stdout, "%g %g %g \n", e,p, q);
    }
    fprintf(outfile, "\n");
  /* Destroy parameter vector(s) */
  glbFreeParams(true_values);
  //glbFreeParams(test_values);
  
  

  // fprintf(outfile, "\n");

  exit(0);

}
