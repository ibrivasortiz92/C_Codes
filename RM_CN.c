#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <crtdbg.h>
#include <unistd.h>
#define printf __mingw_printf
#define eps 1e-18

//////////////////////////////////////////////////////////////////////////////////// PROTOTYPES

///////////////////////////////////////////////////////// FREE MEMORY REQUESTED
void free_memory(long double **ZON, long double **XDOM, long double **YDOM, int **ZMAP, long double **QMAP, long double **MIU, long double **THETA, long double **FI, long double **W, long double **LIST, long double **XVALS, long double **XVECTS, long double **YVALS, long double **YVECTS, long double **RM, long double **PV, long double **FM0, long double **FM1, long double **SM, long double **MFLUX, long double **MFLOW, long double **XFLOW, long double **YFLOW);

//////////////////////////////////////////////////////// LOAD PROBLEM FROM FILE
int input_by_txt(int *N, int *nz, long double **ZON, int *nxr, long double **XDOM, int *nyr, long double **YDOM, int **ZMAP, long double **QMAP, long double *BC, long double *tol, const char *filename);

/////////////////////////////////////////////////////////// GENERATE QUADRATURE
int quad(int N, long double **MIU, long double **THETA, long double **CHI, long double **W, long double **LIST);

////////////////////////////////////////////////////////////////////// SPECTRUM
void myXFunc(int N, long double x, long double MIU[], long double W[], long double c0, long double *y);

void myXRootFunc(int N, long double a, long double b, long double MIU[], long double W[], long double c0, long double *root);

void myYFunc(int N, long double x, long double THETA[], long double W[], long double c0, long double *y);

void myYRootFunc(int N, long double a, long double b, long double THETA[], long double W[], long double c0, long double *root);

int spectrum(int N, long double MIU[], long double THETA[], long double LIST[], long double W[], int nz, long double ZON[], long double **xvals, long double **xvects, long double **yvals, long double **yvects);

void print_spectrum(int N, int nz, long double XVALS[], long double XVECTS[], long double YVALS[], long double YVECTS[]);

////////////////////////////////////////////////////////////// MATRIX FUNCTIONS
void print_vector(int M, long double Vector[]);

void print_matrix(int M, long double Matrix[]);

long double* zeros(int M, long double **OUTPUT);

long double* eye(int M, long double **OUTPUT);

long double* equal(int M, long double Matrix[], long double **OUTPUT);

long double* vector_neg(int M, long double VECTOR[], long double **OUTPUT);

long double* matrix_neg(int M, long double Matrix[], long double **OUTPUT);

long double* matrix_sum(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT);

long double* matrix_mult1(int M, long double Matrix[], long double X[], long double **OUTPUT);

long double* matrix_mult2(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT);

int inv(int M, long double Matrix[], long double **OUTPUT);

long double* vector_concat(int M, long double V1[], long double V2[], long double **OUTPUT);

long double* matrix_concat(int M, long double Matrix1[], long double Matrix2[], long double Matrix3[], long double Matrix4[], long double **OUTPUT);

////////////////////////////////////////////////////////////// RESPONSE MATRIX
int response_matrix(int N, int nz, long double ZON[], int nxr, long double XDOM[], int nyr, long double YDOM[], int ZMAP[], long double QMAP[], long double MIU[], long double THETA[], long double W[], long double XVALS[], long double XVECTS[], long double YVALS[], long double YVECTS[], long double **RM, long double **PV, long double **FM0, long double **FM1, long double **SM);

long double* get_RM(int M, int nyr, int nxr, int ry, int rx, long double RM[], long double **OUTPUT);

long double* get_PV(int M, int nyr, int nxr, int ry, int rx, long double PV[], long double **OUTPUT);

long double* get_FM0(int M, int nyr, int nxr, int ry, int rx, long double FM0[], long double **OUTPUT);

long double* get_FM1(int M, int nyr, int nxr, int ry, int rx, long double FM1[], long double **OUTPUT);

long double* get_SM(int M, int nyr, int nxr, int ry, int rx, long double SM[], long double **OUTPUT);

//////////////////////////////////////////////////////// RM_CN ITERATIVE SCHEME
int rm_cn (int N, int nz, long double ZON[], int nxr, long double XDOM[], int nyr, long double YDOM[], int ZMAP[], long double QMAP[], long double BC[], long double tol, long double W[], long double RM[], long double PV[], long double FM0[], long double FM1[], long double SM[], long double **MFLUX, long double **MFLOW, long double **XFLOW, long double **YFLOW, int *ITER, long double *cpu_time);

////////////////////////////////////////////////////////// AUXILIARY FUNCTIONS
void print_problem(int N, int nz, long double ZON[], int nxr, long double XDOM[], int nyr, long double YDOM[], int ZMAP[], long double QMAP[], long double BC[], long double tol);

void post_processing(int N, int nz, long double ZON[], int nxr, long double XDOM[], int nyr, long double YDOM[], int ZMAP[], long double QMAP[], long double BC[], long double MIU[], long double THETA[], long double W[], long double MFLUX[], long double MFLOW[], long double XFLOW[], long double YFLOW[]);

void json_output(int N, int nxr, long double XDOM[], int nyr, long double YDOM[], int status, int ITER, long double cpu_time, long double MFLUX[], long double MFLOW[], long double XFLOW[], long double YFLOW[]);
                  
///////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////// MAIN
int main(int argc, char *argv[]){

  int status = 0;

  // INPUT VARIABLES
  int N = 0;                   // Quadrature order
  int nz = 0;                  // Number of zones
  long double *ZON = NULL;     // Zone entries
  int nxr = 0;                 // Number of regions in X
  long double *XDOM = NULL;    // X Region entries
  int nyr = 0;                 // Number of regions in Y
  long double *YDOM = NULL;    // Y Region entries
  int *ZMAP = NULL;            // Zone map
  long double *QMAP = NULL;    // External source map
  long double BC[4] = {0.0};   // Boundary conditions
  long double tol = 0.0;       // Tolerance
  const char *filename = NULL; // File name to load

  // QUADRATURE VARIABLES
  long double *MIU = NULL;   // Ordinates in X
  long double *THETA = NULL; // Ordinates in Y
  long double *FI = NULL;    // Ordinates in Z
  long double *W = NULL;     // Weight
  long double *LIST = NULL;  // Ordinates list

  // SPECTRUM VARIABLES
  long double *XVALS = NULL;  // Eigenvalues in X
  long double *XVECTS = NULL; // Eigenvectors in X
  long double *YVALS = NULL;  // Eigenvalues in Y
  long double *YVECTS = NULL; // Eigenvectors in Y

  // RESPONSE MATRIX VARIABLES
  long double *RM = NULL;  // Response matrix
  long double *PV = NULL;  // Particular vector
  long double *FM0 = NULL; // Average flux matrix 0
  long double *FM1 = NULL; // Avergae flux matrix 1
  long double *SM = NULL;  // Average source vector

  // PRINCIPAL VARIABLES
  long double *MFLUX = NULL;  // Escalar flux in the nodes
  long double *MFLOW = NULL;  // Angular flux in the nodes
  long double *XFLOW = NULL;  // Angular flux at the y edges
  long double *YFLOW = NULL;  // Angular flux at the x edges
  int ITER = -1;              // Iterations
  long double cpu_time = 0.0; // CPU time

  if (argc != 2){
    // Usage: ./RM_CN <input>
    status = 1;
    json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);
    return 0;
  }

  // LOAD PROBLEM
  filename = argv[1];
  status = input_by_txt(&N, &nz, &ZON, &nxr, &XDOM, &nyr, &YDOM, &ZMAP, &QMAP, BC, &tol, filename);
  if (status != 0) {
    json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);
    return 0;
  }
  
  // PRINT PROBLEM
  print_problem(N, nz, ZON, nxr, XDOM, nyr, YDOM, ZMAP, QMAP, BC, tol);

  // GENERATE QUADRATURE
  int M = N * (N + 2) / 2;
  status = quad(N, &MIU, &THETA, &FI, &W, &LIST);
  if (status != 0){
    json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);
    free_memory(&ZON, &XDOM, &YDOM, &ZMAP, &QMAP,
                &MIU, &THETA, &FI, &W, &LIST,
                &XVALS, &XVECTS, &YVALS, &YVECTS,
                &RM, &PV, &FM0, &FM1, &SM,
                &MFLUX, &MFLOW, &XFLOW, &YFLOW);
    return 0;
  };

  // PRINT QUADRATURE
  printf("\n\n2. QUADRATURE SETUP:\n\n");
	printf("M\t MIU\t\t THETA\t\t W\n");
	for (int m = 0; m < M; m++) {
		printf("%d\t% .8Lf\t% .8Lf\t% .8Lf\n", m + 1, MIU[m], THETA[m], W[m]);
	}
	printf("\n");

  // GENERATE SPECTRUM
  status = spectrum(N, MIU, THETA, LIST, W, nz, ZON, &XVALS, &XVECTS, &YVALS, &YVECTS);
  if (status != 0){
    json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);
    free_memory(&ZON, &XDOM, &YDOM, &ZMAP, &QMAP,
                &MIU, &THETA, &FI, &W, &LIST,
                &XVALS, &XVECTS, &YVALS, &YVECTS,
                &RM, &PV, &FM0, &FM1, &SM,
                &MFLUX, &MFLOW, &XFLOW, &YFLOW);
    return 0;
  };

  // PRINT SPECTRUM
  print_spectrum(N, nz, XVALS, XVECTS, YVALS, YVECTS);

  // GENERATE RESPONSE MATRIX
  status = response_matrix(N, nz, ZON, nxr, XDOM, nyr, YDOM, ZMAP, QMAP, MIU, THETA, W, XVALS, XVECTS, YVALS, YVECTS, &RM, &PV, &FM0, &FM1, &SM);
  if (status != 0){
    json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);
    free_memory(&ZON, &XDOM, &YDOM, &ZMAP, &QMAP,
                &MIU, &THETA, &FI, &W, &LIST,
                &XVALS, &XVECTS, &YVALS, &YVECTS,
                &RM, &PV, &FM0, &FM1, &SM,
                &MFLUX, &MFLOW, &XFLOW, &YFLOW);
    return 0;
  };

  // PRINT RESPONSE MATRIZ
  // long double *RAUX = NULL;
  // printf("FM0:\n");
  // for (int ry = 0; ry < nyr; ry++){
  //   for (int rx = 0; rx < nxr; rx++){
  //     printf("RY = %d, RX = %d:\n", ry, rx);
  //     RAUX = get_FM0(M, nyr,nxr,ry,rx,FM0,&RAUX);
  //     print_matrix(M,RAUX);
  //   }
  // }

  // ITERATIVE SCHEME
  status = rm_cn (N, nz, ZON, nxr, XDOM, nyr, YDOM, ZMAP, QMAP, BC, tol, W, RM, PV, FM0, FM1, SM, &MFLUX, &MFLOW, &XFLOW, &YFLOW, &ITER, &cpu_time);
  if (status != 0){
    json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);
    free_memory(&ZON, &XDOM, &YDOM, &ZMAP, &QMAP,
                &MIU, &THETA, &FI, &W, &LIST,
                &XVALS, &XVECTS, &YVALS, &YVECTS,
                &RM, &PV, &FM0, &FM1, &SM,
                &MFLUX, &MFLOW, &XFLOW, &YFLOW);
    return 0;
  };

  // POST PROCESSING
  post_processing(N, nz, ZON, nxr, XDOM, nyr, YDOM, ZMAP, QMAP, BC, MIU, THETA, W, MFLUX, MFLOW, XFLOW, YFLOW);

  // JSON OUTPUTS
  json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);

  // FREE MEMORY
  free_memory(&ZON, &XDOM, &YDOM, &ZMAP, &QMAP,
              &MIU, &THETA, &FI, &W, &LIST,
              &XVALS, &XVECTS, &YVALS, &YVECTS,
              &RM, &PV, &FM0, &FM1, &SM,
              &MFLUX, &MFLOW, &XFLOW, &YFLOW);

  _CrtDumpMemoryLeaks();
  return 0;

}
///////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////// IMPLEMENTATION

///////////////////////////////////////////////////////// FREE MEMORY REQUESTED
void free_memory(long double **ZON, long double **XDOM, long double **YDOM, int **ZMAP, long double **QMAP, long double **MIU, long double **THETA, long double **FI, long double **W, long double **LIST, long double **XVALS, long double **XVECTS, long double **YVALS, long double **YVECTS, long double **RM, long double **PV, long double **FM0, long double **FM1, long double **SM, long double **MFLUX, long double **MFLOW, long double **XFLOW, long double **YFLOW){
  
  // INPUT VARIABLES
  if (*ZON != NULL) free(*ZON); if (*XDOM != NULL) free(*XDOM); if (*YDOM != NULL) free(*YDOM); 
  if (*ZMAP != NULL) free(*ZMAP); if (*QMAP != NULL) free(*QMAP);

  // QUADRATURE VARIABLES
  if (*MIU != NULL) free(*MIU); if (*THETA != NULL) free(*THETA); if (*FI != NULL) free(*FI);
  if (*W != NULL) free(*W); if (*LIST != NULL) free(*LIST);

  // SPECTRUM VARIABLES
  if (*XVALS != NULL) free(*XVALS); if (*XVECTS != NULL) free(*XVECTS);
  if (*YVALS != NULL) free(*YVALS); if (*YVECTS != NULL) free(*YVECTS);

  // RESPONSE MATRIX VARIABLES
  if (*RM != NULL) free(*RM); if (*PV != NULL) free(*PV);
  if (*FM0 != NULL) free(*FM0); if (*FM1 != NULL) free(*FM1); if (*SM != NULL) free(*SM);
  
  // PRINCIPAL VARIABLES
  if (*MFLUX != NULL) free(*MFLUX); if (*MFLOW != NULL) free(*MFLOW);
  if (*XFLOW != NULL) free(*XFLOW); if (*YFLOW != NULL) free(*YFLOW);
  
}
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////// LOAD PROBLEM FROM FILE
int input_by_txt(int *N,                  // Quadrature order
                int *nz,                  // Number of zones
                long double **ZON,        // Zone entries
                int *nxr,                 // Number of regions in X
                long double **XDOM,               // X Region entries
                int *nyr,                 // Number of regions in Y
                long double **YDOM,               // Y Region entries
                int **ZMAP,               // Zone map
                long double **QMAP,       // External source map
                long double *BC,          // Boundary conditions
                long double *tol,         // Tolerance
                const char *filename      // File name to load
                ){

  FILE *cfPtr = NULL;
  if ((cfPtr = fopen(filename, "r")) == NULL){
      //Error opening input file
      return 2;
  }

  // QUADRATURE ORDER
  if (fscanf(cfPtr, "%d", N) == 0){
      // Error reading input file
      return 2;
  }
	if (*N < 2 || *N > 18 || (*N)%2 != 0){
		// Invalid quadrature order (N)
	  // Value must be an even number between 2 and 18
		return 2;
	}

  // NUMBER OF ZONES
  if (fscanf(cfPtr, "%d", nz) == 0){
      // Error reading input file
      return 2;
  }
	if (*nz < 1){
		// Invalid number of zones (NZ)
		// Value must be greater than 1
		return 2;
	}

  // ZONE FILLING
  *ZON = malloc(sizeof(long double) * (*nz * 2)); if(*ZON == NULL) return 3;
  double st, ss;
	int zone_count = 0;
  for (int z = 0; z < *nz; z++){
      if (fscanf(cfPtr, "%lf %lf", &st, &ss) == 0){
          // Error reading input file
          free(*ZON);
          return 2;
      }
      (*ZON)[z * 2] = (long double)st; (*ZON)[z * 2 + 1] = (long double)ss;
      zone_count = zone_count + 1;
  }
	if (zone_count != *nz){
		free(*ZON);
		// Error reading the entries of ZON:
		// Must be provided the same number of rows as the number of zones
		return 2;
	}

  // NUMBER OF REGIONS IN X
  if (fscanf(cfPtr, "%d", nxr) == 0){
    //Error reading input file
    free(*ZON);
    return 2;
  }
	if (*nxr < 1){
		free(*ZON);
		// Invalid number of regions in X (NR_X):
		// Value must be greater than 1
		return 2;
	}

  // REGION FILLING IN X
  *XDOM = malloc(sizeof(long double) * (*nxr * 2)); 
  if(*XDOM == NULL){
    free(*ZON); return 3;
  }
  int rx_count = 0;
  double len, nodes;
  for (int xr = 0; xr < *nxr; xr++){
    if (fscanf(cfPtr, "%lf %lf", &len, &nodes) == 0){
      // Error reading input file
      free(*ZON); free(*XDOM);
      return 2;
    }
    (*XDOM)[xr * 2] = (long double)len; (*XDOM)[xr * 2 + 1] = (long double)nodes;
    rx_count = rx_count + 1;
  }
	if (rx_count != *nxr){
		free(*ZON); free(*XDOM);
		// Error reading the entries of XDOM:
		// Must be provided the same number of rows as the number of regions in X
		return 2;
	}

  // NUMBER OF REGIONS IN Y
  if (fscanf(cfPtr, "%d", nyr) == 0){
    // Error reading input file
    free(*ZON); free(*XDOM);
    return 2;
  }
	if (*nxr < 1){
		free(*ZON); free(*XDOM);
		// Invalid number of regions in Y (NR_Y):
		// Value must be greater than 1.");
		return 2;
	}

  // REGION FILLING IN Y
  *YDOM = malloc(sizeof(long double) * (*nyr * 2));
  if(*YDOM == NULL){
    free(*ZON); free(*XDOM); return 3;
  }
	int ry_count = 0;
  for (int yr = 0; yr < *nyr; yr++){
    if (fscanf(cfPtr, "%lf %lf", &len, &nodes) == 0){
      // Error reading input file
      free(*ZON); free(*XDOM); free(*YDOM);
      return 2;
    }
    (*YDOM)[yr * 2] = (long double)len; (*YDOM)[yr * 2 + 1] = (long double)nodes;
    ry_count = ry_count + 1;
  }
	if (ry_count != *nyr){
		free(*ZON); free(*XDOM); free(*YDOM);
		// Error reading the entries of YDOM:\n\n"
		// Must be provided the same number of rows as the number of regions in Y.");
		return 2;
	}

  // ZONE MAPPING
	*ZMAP = malloc(sizeof(int) * (*nyr) * (*nxr));
  if(*ZMAP == NULL){
    free(*ZON); free(*XDOM); free(*YDOM); return 3;
  }
	int z, entry_z_count = 0;
  for (int yr = 0; yr < *nyr; yr++) {
		for (int xr = 0; xr < *nxr; xr++) {
			if (fscanf(cfPtr, "%d", &z) == 0) {
				//Error reading input file!");
				free(*ZON); free(*XDOM); free(*YDOM);
				free(*ZMAP);
        return 2;
			}
			(*ZMAP)[yr * (*nxr) + xr] = z - 1;
			entry_z_count = entry_z_count + 1;
		}
  }
	if (entry_z_count != (*nyr) * (*nxr)){
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP);
		// Error reading the entries of ZMAP:
		// The entries in ZMAP don't match with the regions of the problem
		return 2;
	}

  // EXTERNAL SOURCE MAPPING
	*QMAP = malloc(sizeof(long double) * (*nyr) * (*nxr)); 
  if(*QMAP == NULL){
    free(*ZON); free(*XDOM); free(*YDOM); free(*ZMAP); return 3;
  }
	double q;
	int entry_q_count = 0;
  for (int yr = 0; yr < *nyr; yr++) {
		for (int xr = 0; xr < *nxr; xr++) {
			if (fscanf(cfPtr, "%lf", &q) == 0) {
				// Error reading input file
				free(*ZON); free(*XDOM); free(*YDOM);
				free(*ZMAP); free(*QMAP);
        return 2;
			}
			(*QMAP)[yr * (*nxr) + xr] = (long double)q;
			entry_q_count = entry_q_count + 1;
		}
  }
	if (entry_q_count != (*nyr) * (*nxr)){
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
		// Error reading the entries of QMAP:
		// The entries in QMAP don't match with the regions of the problem
		return 2;
	}

  // BOUNDARY CONDITIONS
  double cond;
	int bc_count = 0;
	for (int c = 0; c < 4; c++) {
		if (fscanf(cfPtr, "%lf", &cond) == 0) {
			// Error reading input file
			free(*ZON); free(*XDOM); free(*YDOM);
			free(*ZMAP); free(*QMAP);
      return 2;
		}
		BC[c] = (long double)cond;
		bc_count = bc_count + 1;
		if (cond < 0.0 ){
			if (cond != -1.0){
				free(*ZON); free(*XDOM); free(*YDOM);
			  free(*ZMAP); free(*QMAP);
				// Invalid boudary condition:
		    // Value must be -1.0, 0.0, or a positive number
				return 2;
			}
		}
	}
	if (bc_count != 4){
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
		// The number of boundary conditions must be 4
		return 2;
	}

  // TOLERANCE
  double dtol;
	if (fscanf(cfPtr, "%lf", &dtol) == 0) {
		// Error reading input file
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
    return 2;
	}
	if (dtol <= 0.0 || dtol >= 1.0){
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
		// The tolerance (TOL) must be a small number between 0.0 and 1.0
		return 2;
	}
  *tol = (long double)dtol;

	fclose(cfPtr);

  return 0;
}
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////// GENERATE QUADRATURE
int quad(int N,                // Quadrature order
         long double **MIU,    // Ordinates in X
         long double **THETA,  // Ordinates in Y
         long double **FI,     // Ordinate list
         long double **W,      // Weight
         long double **LIST    // Ordinate list
         ){
    
  // DIRECTIONS IN THE XY PLANE
  int M = N * (N + 2) / 2;
  
  // ALLOCATE MEMORY FOR THE ORDINATES
  *MIU = malloc(sizeof(long double) * M);
  *THETA = malloc(sizeof(long double) * M);
  *FI = malloc(sizeof(long double) * M);
  *W = malloc(sizeof(long double) * M);
  *LIST = malloc(sizeof(long double) * N / 2);
  if (*MIU == NULL || *THETA == NULL || *FI == NULL || *W == NULL || *LIST == NULL) return 3;

  // PREFIXED VALUES
  long double chi[10] = { 0.0 }, wlist[10] = { 0.0 };
	if (N == 2) {
		chi[0] = 0.5773502692; wlist[0] = 1.0;
	}
	else if (N == 4) {
		chi[0] = 0.3500212; wlist[0] = (long double)1/3;
		chi[1] = 0.8688903;
	}
	else if (N == 6) {
		chi[0] = 0.2666355; wlist[0] = 0.1761263;
		chi[1] = 0.6815076; wlist[1] = 0.1572071;
		chi[2] = 0.9261808;
	}
	else if (N == 8) {
		chi[0] = 0.2182179; wlist[0] = 0.1209877;
		chi[1] = 0.5773503; wlist[1] = 0.0907407;
		chi[2] = 0.7867958; wlist[2] = 0.0925926;
		chi[3] = 0.9511897;
	}
	else if (N == 10) {
		chi[0] = 0.1893213; wlist[0] = 0.0893031;
		chi[1] = 0.5088818; wlist[1] = 0.0725292;
		chi[2] = 0.6943189; wlist[2] = 0.0450438;
		chi[3] = 0.8397600; wlist[3] = 0.0539281;
		chi[4] = 0.9634910;
	}
	else if (N == 12) {
		chi[0] = 0.1672126; wlist[0] = 0.0707626;
		chi[1] = 0.4595476; wlist[1] = 0.0558811;
		chi[2] = 0.6280191; wlist[2] = 0.0373377;
		chi[3] = 0.7600210; wlist[3] = 0.0502819;
		chi[4] = 0.8722706; wlist[4] = 0.0258513;
		chi[5] = 0.9716377;
	}
	else if (N == 14) {
		chi[0] = 0.1519859; wlist[0] = 0.0579970;
		chi[1] = 0.4221570; wlist[1] = 0.0489008;
		chi[2] = 0.5773503; wlist[2] = 0.0221497;
		chi[3] = 0.6988921; wlist[3] = 0.0407009;
		chi[4] = 0.8022263; wlist[4] = 0.0393867;
		chi[5] = 0.8936911; wlist[5] = 0.0245518;
		chi[6] = 0.9766272; wlist[6] = 0.0121325;
	}
	else if (N == 16) {
		chi[0] = 0.1389568; wlist[0] = 0.0489872;
		chi[1] = 0.3922893; wlist[1] = 0.0413296;
		chi[2] = 0.5370966; wlist[2] = 0.0212326;
		chi[3] = 0.6504264; wlist[3] = 0.0256207;
		chi[4] = 0.7467506; wlist[4] = 0.0360486;
		chi[5] = 0.8319966; wlist[5] = 0.0144589;
		chi[6] = 0.9092855; wlist[6] = 0.0344958;
		chi[7] = 0.9805009; wlist[7] = 0.0085179;
	}
	else if (N == 18) {
		chi[0] = 0.1293445; wlist[0] = 0.0422646;
		chi[1] = 0.3680438; wlist[1] = 0.0376127;
		chi[2] = 0.5041652; wlist[2] = 0.0066907;
		chi[3] = 0.6106625; wlist[3] = 0.0391919;
		chi[4] = 0.7011669; wlist[4] = 0.0042550;
		chi[5] = 0.7812562; wlist[5] = 0.0423662;
		chi[6] = 0.8538662; wlist[6] = 0.0092396;
		chi[7] = 0.9207680; wlist[7] = 0.0156648;
		chi[8] = 0.9831277; wlist[8] = 0.0136576;
		                    wlist[9] = 0.0139903;
  }

  // ORDINATES FILLING
	int d = 0, nlevel = N / 2, aux0;
	for (int n = 0; n < nlevel; n++) {
		aux0 = 0;
		for (int m = (nlevel - 1) - n; m >= 0; m--) {
			(*MIU)[d] = chi[m]; (*THETA)[d] = chi[aux0]; 
      (*FI)[d] = chi[n];  (*W)[d] = 0.0;
			d = d + 1; aux0 = aux0 + 1;
		}
	}

  // WEIGHT FILLING
	int p = 0; long double aux1, aux2, aux3;
	for (int n = 0; n < M / 4; n++) {
		if ((*W)[n] == 0.0) {
			(*W)[n] = wlist[p];
			aux1 = (*MIU)[n]; aux2 = (*THETA)[n]; aux3 = (*FI)[n];
			for (int m = 0; m < M / 4; m++) {
				if (aux1 == (*THETA)[m] && aux2 == (*MIU)[m] && aux3 == (*FI)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*THETA)[m] && aux2 == (*FI)[m] && aux3 == (*MIU)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*FI)[m] && aux2 == (*THETA)[m] && aux3 == (*MIU)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*FI)[m] && aux2 == (*MIU)[m] && aux3 == (*THETA)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*MIU)[m] && aux2 == (*FI)[m] && aux3 == (*THETA)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*MIU)[m] && aux2 == (*THETA)[m] && aux3 == (*FI)[m]) {
					(*W)[m] = wlist[p];
				}
			}
			p = p + 1;
		}
	}

  // QUADRANT FILLING
	int aux4 = 0;
	for (int q = 1; q <= 4; q++) {
		for (int m = 0; m < M / 4; m++) {
			if (q == 2) {
				(*MIU)[aux4] = -(*MIU)[m]; (*THETA)[aux4] = (*THETA)[m];
				(*FI)[aux4] = (*FI)[m]; (*W)[aux4] = (*W)[m];
			}
			else if (q == 3) {
				(*MIU)[aux4] = -(*MIU)[m]; (*THETA)[aux4] = -(*THETA)[m];
				(*FI)[aux4] = (*FI)[m]; (*W)[aux4] = (*W)[m];
			}
			else if (q == 4) {
				(*MIU)[aux4] = (*MIU)[m]; (*THETA)[aux4] = -(*THETA)[m];
				(*FI)[aux4] = (*FI)[m]; (*W)[aux4] = (*W)[m];
			}
			aux4 = aux4 + 1;
		}
	}

  // ORDINATE LIST FILL
  for (int i = 0; i < N / 2; i++){
    (*LIST)[i] = chi[i];
  }

  return 0;
    
}
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////// SPECTRUM FUNCTIONS
void myXFunc(int N,              // Quadrature order
             long double x,      // X value
             long double MIU[],  // Ordinate in X
             long double W[],    // Weight
             long double c0,     // Scattering Ratio
             long double *y)     // Y value
             {

  int M = N * (N + 2) / 2;

  (*y) = 0.0;
  for (int m = 0; m < M; m++){
    long double miu = MIU[m], w = W[m];
    (*y) = (*y) + w / (miu + x);
  }
  (*y) = 0.25 * c0 * x * (*y) - 1;

}

void myXRootFunc(int N,              // Quadrature order
                 long double a,      // Left side
                 long double b,      // Right side
                 long double MIU[],  // Ordinate in X
                 long double W[],    // Weight
                 long double c0,     // Scattering ratio
                 long double *root)  // Root
                 {

  long double yb;
  myXFunc(N, b, MIU, W, c0, &yb);
  while(fabs(b - a) > eps){
    long double yc, c = 0.5 * (a + b);
    myXFunc(N, c, MIU, W, c0, &yc);
    if (yc == 0.0){
      a = c;
      b = c;
    }
    else if (yb * yc > 0.0) {
      b = c;
      yb = yc;
    }
    else a = c;
  }
  (*root) = 0.5 * (a + b);

}

void myYFunc(int N,                // Quadrature order
             long double x,        // X value
             long double THETA[],  // Ordinate in Y
             long double W[],      // Weight
             long double c0,       // Scattering Ratio
             long double *y)       // Y value
             {

  int M = N * (N + 2) / 2;

  (*y) = 0.0;
  for (int m = 0; m < M; m++){
    long double theta = THETA[m], w = W[m];
    (*y) = (*y) + w / (theta + x);
  }
  (*y) = 0.25 * c0 * x * (*y) - 1;

}

void myYRootFunc(int N,                // Quadrature order
                 long double a,        // Left side
                 long double b,        // Right side
                 long double THETA[],  // Ordinate in Y
                 long double W[],      // Weight
                 long double c0,       // Scattering ratio
                 long double *root)    // Root
                 {

  long double yb;
  myYFunc(N, b, THETA, W, c0, &yb);
  while(fabs(b - a) > eps){
    long double yc, c = 0.5 * (a + b);
    myYFunc(N, c, THETA, W, c0, &yc);
    if (yc == 0.0){
      a = c;
      b = c;
    }
    else if (yb * yc > 0.0) {
      b = c;
      yb = yc;
    }
    else a = c;
  }
  (*root) = 0.5 * (a + b);

}

int spectrum(int N,                // Quadrature order
            long double MIU[],    // Ordinate in X
            long double THETA[],  // Ordinate in Y
            long double LIST[],   // Ordinate list
            long double W[],      // Weight
            int nz,               // Number of zones
            long double ZON[],    // Zone entries
            long double **xvals,  // Eigenvalues in X
            long double **xvects, // Eigenvectors in X
            long double **yvals,  // Eigenvalues in Y
            long double **yvects  // Eigenvectors in Y
            ){

  // DIRECTIONS IN THE XY PLANE
  int M = N * (N + 2) / 2;

  // MEMORY ALLOCATION
  *xvals = malloc(sizeof(long double) * M * nz);
  *xvects = malloc(sizeof(long double) * M * M * nz);
  *yvals = malloc(sizeof(long double) * M * nz);
  *yvects = malloc(sizeof(long double) * M * M * nz);
  if (*xvals == NULL || *xvects == NULL || *yvals == NULL || *yvects == NULL) return 3;

  // AUXILIARY EIGENVALUES IN X
  long double *xvals1 = NULL, *xvals2 = NULL;
  xvals1 = malloc(sizeof(long double) * N);
  xvals2 = malloc(sizeof(long double) * (M - N));

  // AUXILIARY EIGENVALUES IN Y
  long double *yvals1 = NULL, *yvals2 = NULL;
  yvals1 = malloc(sizeof(long double) * N);
  yvals2 = malloc(sizeof(long double) * (M - N));

  // AUXILIARY EIGENVECTORS IN X
  long double *xvects1 = NULL, *xvects2 = NULL;
  xvects1 = malloc(sizeof(long double) * M * N);
  xvects2 = malloc(sizeof(long double) * M * (M - N));

  // AUXILIARY EIGENVECTORS IN Y
  long double *yvects1 = NULL, *yvects2 = NULL;
  yvects1 = malloc(sizeof(long double) * M * N);
  yvects2 = malloc(sizeof(long double) * M * (M - N));

  int *aux = NULL;
  aux = malloc(sizeof(int) * M);

  if (xvals1 == NULL || xvals2 == NULL || yvals1 == NULL || yvals2 == NULL || xvects1 == NULL || xvects2 == NULL || yvects1 == NULL || yvects2 == NULL || aux == NULL){
    if(xvals1 != NULL) free(xvals1); if(xvals2 != NULL) free(xvals2);
    if(yvals1 != NULL) free(yvals1); if(yvals2 != NULL) free(yvals2);
    if(xvects1 != NULL) free(xvects1); if(xvects2 != NULL) free(xvects2);
    if(yvects1 != NULL) free(yvects1); if(yvects2 != NULL) free(yvects2);
    if (aux != NULL) free(aux); return 3;
  }

  // BODY
  for (int z = 0; z < nz; z++){

    long double st = ZON[z * 2], ss = ZON[z * 2 + 1];
    long double c0 = ss / st;

    // EIGENVALUES
    if (c0 != 0.0){

      // DISPERSION LAW IN X
      for(int i = 0; i < N / 2; i++){

        long double chi_i = LIST[i], chi_f, h, a, b, r;

        if (i == N / 2 - 1) chi_f = 10;
        else chi_f = LIST[i + 1];

        h = (chi_f - chi_i) / (pow(10, N));
        a = chi_i + h;
        b = chi_f - h;

        myXRootFunc(N, a, b, MIU, W, c0, &r);
        xvals1[i] = r;
        
      }

      // ZERO NORMALIZATION IN X
      int k = 0;
      for (int i = 0; i < N / 2; i++){

        long double val = LIST[i];
        int xmult = 0;

        for (int m = 0; m < M; m++){
          long double miu = MIU[m];
          if (val == miu) xmult = xmult + 1;
        }
        xmult = xmult - 1;

        while (xmult > 0){
          xvals2[k] = val;
          xmult = xmult - 1;
          k = k + 1; 
        }

      }

      // DISPERSION LAW IN Y
      for(int i = 0; i < N / 2; i++){

        long double chi_i = LIST[i], chi_f, h, a, b, r;

        if (i == N / 2 - 1) chi_f = 10;
        else chi_f = LIST[i + 1];

        h = (chi_f - chi_i) / (pow(10, N));
        a = chi_i + h;
        b = chi_f - h;

        myYRootFunc(N, a, b, THETA, W, c0, &r);
        yvals1[i] = r;
        
      }

      // ZERO NORMALIZATION IN Y
      k = 0;
      for (int i = 0; i < N / 2; i++){

        long double val = LIST[i];
        int ymult = 0;

        for (int m = 0; m < M; m++){
          long double theta = THETA[m];
          if (val == theta) ymult = ymult + 1;
        }
        ymult = ymult - 1;

        while (ymult > 0){
          yvals2[k] = val;
          ymult = ymult - 1;
          k = k + 1; 
        }

      }

      // ORDERING EIGENVALUES
      long double temp;
      for (int i = 0; i < N / 2; i++){
        for (int j = 0; j < N / 2; j++){
          if (xvals1[i] > xvals1[j]) {
            temp = xvals1[i];
            xvals1[i] = xvals1[j];
            xvals1[j] = temp;
          }
          if (yvals1[i] > yvals1[j]) {
            temp = yvals1[i];
            yvals1[i] = yvals1[j];
            yvals1[j] = temp;
          }
        }
      }
      for (int i = 0; i < N/2; i++){
        xvals1[N / 2 + i] = - xvals1[i];
        yvals1[N / 2 + i] = - yvals1[i];
      }
      for (int i = 0; i < (M - N) / 2; i++){
        xvals2[(M - N) / 2 + i] = - xvals2[i];
        yvals2[(M - N) / 2 + i] = - yvals2[i];
      }

      // EIGENVALUE ASSIGNMENT
      for (int i = 0; i < M; i++){
        if (i < (M - N)){
          (*xvals)[M * z + i] = xvals2[i];
          (*yvals)[M * z + i] = yvals2[i];
        }
        else {
          (*xvals)[M * z + i] = xvals1[i - (M - N)];
          (*yvals)[M * z + i] = yvals1[i - (M - N)];
        }
      }

    }

    // c0 = 0
    else {

      for (int i = 0; i < M; i++){
        long double miu = MIU[i], theta = THETA[i];
        (*xvals)[M * z + i] = - miu;
        (*yvals)[M * z + i] = - theta;
      }

    }
    
    // EIGENVECTORS
    if (c0 != 0.0){

      // EIGENVECTORS CALCULATION BY DISPERSION LAW
      for (int i = 0; i < N; i++){
        for (int m = 0; m < M; m++){
          long double miu = MIU[m], theta = THETA[m];
          xvects1[M * i + m] = 0.25 * c0 * xvals1[i] / (miu + xvals1[i]);
          yvects1[M * i + m] = 0.25 * c0 * yvals1[i] / (theta + yvals1[i]);
        }
      }

      for (int i = 0; i < (M - N); i++){
        for (int m = 0; m < M; m++){
          xvects2[M * i + m] = 0.0;
          yvects2[M * i + m] = 0.0;
        }
      }

      // EIGENVECTORS CALCULATION BY ZERO NORMALIZATION IN X
      for (int m = 0; m < M; m++) aux[m] = 0;
      for (int i = 0; i < (M - N); i++){
        long double val = xvals2[i];
        for (int m = 0; m < M; m++){
          long double mw = W[m], mmiu = MIU[m];
          if (val == - mmiu && aux[m] == 0){
            aux[m] = 1;
            for(int n = 0; n < M; n++){
              long double nw = W[n], nmiu = MIU[n];
              if (val == - nmiu && aux[n] == 0){
                aux[n] = 1;
                xvects2[M * i + n] = - mw / nw;
                xvects2[M * i + m] = 1.0;
                aux[m] = 0;
                break;
              }
            }
            break;
          }
        }
      }

      // EIGENVECTORS CALCULATION BY ZERO NORMALIZATION IN Y
      for (int m = 0; m < M; m++) aux[m] = 0;
      for (int i = 0; i < (M - N); i++){
        long double val = yvals2[i];
        for (int m = 0; m < M; m++){
          long double mw = W[m], mtheta = THETA[m];
          if (val == - mtheta && aux[m] == 0){
            aux[m] = 1;
            for(int n = 0; n < M; n++){
              long double nw = W[n], ntheta = THETA[n];
              if (val == - ntheta && aux[n] == 0){
                aux[n] = 1;
                yvects2[M * i + n] = - mw / nw;
                yvects2[M * i + m] = 1;
                aux[m] = 0;
                break;
              }
            }
            break;
          }
        }
      }

      // EIGENVECTOR ASSIGNMENT
      for (int i = 0; i < M; i++){
        for (int m = 0; m < M; m++){
          if (i < (M - N)){
            (*xvects)[nz * (M * i + m) + z] = xvects2[M * i + m];
            (*yvects)[nz * (M * i + m) + z] = yvects2[M * i + m];
          }
          else {
            (*xvects)[nz * (M * i + m) + z] = xvects1[M * (i - (M - N)) + m];
            (*yvects)[nz * (M * i + m) + z] = yvects1[M * (i - (M - N)) + m];
          }
        }
      }

    }
      
    // c0 = 0
    else {

      for (int i = 0; i < M; i++){
        for (int m = 0; m < M; m++){
          if (m == i){
            (*xvects)[nz * (M * i + m) + z] = 1.0;
            (*yvects)[nz * (M * i + m) + z] = 1.0;
          }
          else {
            (*xvects)[nz * (M * i + m) + z] = 0.0;
            (*yvects)[nz * (M * i + m) + z] = 0.0;
          }
        }
      }

      printf("OK\n");

    }

  } // end z

  // FREE MEMORY
  free(xvals1); free(xvals2); free(yvals1); free(yvals2);
  free(xvects1); free(xvects2); free(yvects1); free(yvects2); free(aux);

  return 0;

}

void print_spectrum(int N, int nz, long double XVALS[], long double XVECTS[], long double YVALS[], long double YVECTS[]) {

  int M = N * (N + 2) / 2;

  printf ("EIGENVALUES IN X:\n");
  for (int z = 0; z < nz; z++){
    printf("z = %d: ", z + 1);
    for (int i = 0; i < M; i++){
      printf("%.5Lf ", XVALS[M * z + i]);
    }
    printf("\n");
  }
  printf("\n");

  printf("EIGENVECTORS IN X:\n");
  for (int z = 0; z < nz; z++){
    for (int m = 0; m < M; m++){
      for (int i = 0; i < M; i++){
        printf("%.5Lf ", XVECTS[nz * (M * i + m) + z]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");

  printf ("EIGENVALUES IN Y:\n");
  for (int z = 0; z < nz; z++){
    printf("z = %d: ", z + 1);
    for (int i = 0; i < M; i++){
      printf("%.5Lf ", YVALS[M * z + i]);
    }
    printf("\n");
  }
  printf("\n");

  printf("EIGENVECTORS IN Y:\n");
  for (int z = 0; z < nz; z++){
    for (int m = 0; m < M; m++){
      for (int i = 0; i < M; i++){
        printf("%.5Lf ", YVECTS[nz * (M * i + m) + z]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");

}
///////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////// MATRIX FUNCTIONS
void print_vector(int M, long double Vector[]){

  for (int i = 0; i < M; i++){
    printf("%.5Lf\n", Vector[i]);
  }
  printf("\n");

}


void print_matrix(int M, long double Matrix[]){

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      printf("%.5Lf ", Matrix[M * j + i]);
    }
    printf("\n");
  }
  printf("\n");

}


long double* zeros(int M, long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      (*OUTPUT)[M * j + i] = 0.0;
    }
  }

  return *OUTPUT;

}


long double* eye(int M, long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      if (j == i) (*OUTPUT)[M * j + i] = 1.0;
      else (*OUTPUT)[M * j + i] = 0.0;
    }
  }

  return *OUTPUT;

}


long double* equal(int M, long double Matrix[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }
  
  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      (*OUTPUT)[M * j + i] = Matrix[M * j + i];
    }
  }

  return *OUTPUT;

}


long double* vector_neg(int M, long double VECTOR[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M);
  }

  for (int i = 0; i < M; i++){
    if (fabs(VECTOR[i]) < eps) (*OUTPUT)[i] = 0.0;
    else (*OUTPUT)[i] = - VECTOR[i];
  }

  return *OUTPUT;

}


long double* matrix_neg(int M, long double Matrix[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      if (fabs(Matrix[M * j + i]) < eps) (*OUTPUT)[M * j + i] = 0.0;
      else (*OUTPUT)[M * j + i] = - Matrix[M * j + i];
    }
  }

  return *OUTPUT;

}


long double* vector_sum(int M, long double Vector1[], long double Vector2[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M);
  }

  for (int i = 0; i < M; i++){
    (*OUTPUT)[i] = Vector1[i] + Vector2[i];
  }

  return *OUTPUT;

}


long double* matrix_sum(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      (*OUTPUT)[M * j + i] = Matrix1[M * j + i] + Matrix2[M * j + i];
    }
  }

  return *OUTPUT;

}


long double* matrix_mult1(int M, long double Matrix[], long double X[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M);
  }

  long double sum;

  for (int j = 0; j < M; j++){
    sum = 0.0;
    for (int i = 0; i < M; i++){
      sum = sum + Matrix[M * j + i] * X[i];
    }
    (*OUTPUT)[j] = sum;
  }

  return *OUTPUT;

}


long double* matrix_mult2(int M, long double Matrix1[], long double Matrix2[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  long double sum;

  for (int k = 0; k < M; k++) {
    for (int j = 0; j < M; j++){
      sum = 0;
      for (int i = 0; i < M; i++){
        sum = sum + Matrix1[M * j + i] * Matrix2[M * i + k];
      }
      (*OUTPUT)[M * j + k] = sum;
    }
  }

  return *OUTPUT;

}


int inv(int M, long double Matrix[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  int p;
  long double val, m;

  long double *temp = NULL, *IDEN = NULL, *Matrix2 = NULL;
  temp = malloc(sizeof(long double) * M);
  IDEN = malloc(sizeof(long double) * M * M);
  Matrix2 = malloc(sizeof(long double) * M * M);
  if (temp == NULL || IDEN == NULL || Matrix2 == NULL){
    free(temp); free(IDEN); free(Matrix2); return 3;
  }
  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      Matrix2[M * j + i] = Matrix[M * j + i];
      if (i == j) IDEN[M * j + i] = 1.0;
      else IDEN[M * j + i] = 0.0;
    }
  }

  for (int j = 0; j < M; j++){

    // PARTIAL PIVOTING
    val = 0;
    for( int k = j; k < M; k++){
      if (val < fabs(Matrix2[M * k + j])){
        val = fabs(Matrix2[M * k + j]);
        p = k;
      }
    }

    for (int i = 0; i < M; i++){
      temp[i] = Matrix2[M * p + i];
      Matrix2[M * p + i] = Matrix2[M * j + i];
      Matrix2[M * j + i] = temp[i];

      temp[i] = IDEN[M * p + i];
      IDEN[M * p + i] = IDEN[M * j + i];
      IDEN[M * j + i] = temp[i];
    }

    if (fabs(Matrix2[M * j + j]) <= eps){
      free(temp); free(IDEN); free(Matrix2);
      free(*OUTPUT); *OUTPUT = NULL;
      return 4;
    }

    // GAUSS ELIMINATION
    for (int jj = j + 1; jj < M; jj++){
      m = Matrix2[M * jj + j] / Matrix2[M * j + j];
      for(int i = 0; i < M; i++){
        Matrix2[M * jj + i] = Matrix2[M * jj + i] - m * Matrix2[M * j + i];
        IDEN[M * jj + i] = IDEN[M * jj + i] - m * IDEN[M * j + i];
      }
    }
  }

  // BACK SUBSTITUTION
  long double sum;
  for (int k = 0; k < M; k++){
    (*OUTPUT)[M * (M - 1) + k] = IDEN[M * (M - 1) + k] / Matrix2[M * M - 1];
    for (int j = M - 2; j >= 0; j--){
      sum = 0;
      for (int i = j; i < M - 1; i++){
        sum = sum + Matrix2[M * j + i + 1] * (*OUTPUT)[M * (i + 1) + k];
      }
      (*OUTPUT)[M * j + k] = (IDEN[M * j + k] - sum) / Matrix2[M * j + j];
      if (fabs((*OUTPUT)[M * j + k]) < eps) (*OUTPUT)[M * j + k] = 0.0;
    }
  }
  
  free(temp); free(IDEN); free(Matrix2);

  return 0;
}

long double* vector_concat(int M, long double V1[], long double V2[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * 2 * M);
  }

  for (int i = 0; i < M; i++){
    // 1
    (*OUTPUT)[i] = V1[i];

    // 2
    (*OUTPUT)[M + i] = V2[i];
  }

  return *OUTPUT;

}


long double* matrix_concat(int M, long double Matrix1[], long double Matrix2[], long double Matrix3[], long double Matrix4[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * 4 * M * M);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){
      // 1
      (*OUTPUT)[2 * M * j + i] = Matrix1[M * j + i];

      // 2
      (*OUTPUT)[2 * M * j + M + i] = Matrix2[M * j + i];

      // 3
      (*OUTPUT)[2 * M * (M + j) + i] = Matrix3[M * j + i];

      // 4
      // 3
      (*OUTPUT)[2 * M * (M + j) + M + i] = Matrix4[M * j + i];
    }
  }

  return *OUTPUT;

}

///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////// RESPONSE MATRIX FUNCTION
int response_matrix(int N,               // Quadrature order
                     int nz,              // Number of zones
                     long double ZON[],   // Zone entries
                     int nxr,             // Number of regions in X
                     long double XDOM[],  // X region entries
                     int nyr,             // Number of regions in Y
                     long double YDOM[],  // Y regions entries
                     int ZMAP[],          // Zone mapping
                     long double QMAP[],  // External source mapping
                     long double MIU[],   // Ordinates in X
                     long double THETA[], // Ordinates in Y
                     long double W[],     // Weight
                     long double XVALS[], // Eigenvalues of X
                     long double XVECTS[],// Eigenvectors of X 
                     long double YVALS[], // Eigenvalues of Y
                     long double YVECTS[],// Eigenvectors of Y
                     long double **RM,    // Response matrix
                     long double **PV,    // Particular vector
                     long double **FM0,   // Average flux matrix 0
                     long double **FM1,   // Average flux matrix 1
                     long double **SM     // Average Source
                     ){
  
  int M = N * (N + 2) / 2;
  int vcond = 0;

  // OUTPUT MATRICES
  *RM = malloc(sizeof(long double) * 4 * M * M * nyr * nxr);
  *PV = malloc(sizeof(long double) * 2 * M * nyr * nxr);
  *FM0 = malloc(sizeof(long double) * M * M * nyr * nxr);
  *FM1 = malloc(sizeof(long double) * M * M * nyr * nxr);
  *SM = malloc(sizeof(long double) * M * nyr * nxr);
  if(*RM == NULL || *PV == NULL || *FM0 == NULL || *FM1 == NULL || *SM == NULL) return 3;

  // AUXILIARY MATRICES
  long double *XA = NULL, *XE = NULL, *YA = NULL, *YE = NULL;
  XA = malloc(sizeof(long double) * M * M);
  XE = malloc(sizeof(long double) * M * M);
  YA = malloc(sizeof(long double) * M * M); 
  YE = malloc(sizeof(long double) * M * M);

  long double *XB = NULL, *YB = NULL;
  XB = malloc(sizeof(long double) * M * M);
  YB = malloc(sizeof(long double) * M * M);

  long double *H = NULL;
  H = malloc(sizeof(long double) * M);

  long double *XE_INV = NULL, *YE_INV = NULL;
  XE_INV = malloc(sizeof(long double) * M * M);
  YE_INV = malloc(sizeof(long double) * M * M);

  long double *M1 = NULL, *M2 = NULL, *M3 = NULL, *M4 = NULL, *TEMP0 = NULL, *TEMP1 = NULL;
  M1 = malloc(sizeof(long double) * M * M);
  M2 = malloc(sizeof(long double) * M * M);
  M3 = malloc(sizeof(long double) * M * M);
  M4 = malloc(sizeof(long double) * M * M);
  TEMP0 = malloc(sizeof(long double) * M * M);
  TEMP1 = malloc(sizeof(long double) * M * M);

  long double *AUX = NULL, *AUX2 = NULL, *AUX_INV = NULL, *RAUX = NULL, *AUX3 = NULL;
  long double *S0 = NULL, *S1 = NULL;
  AUX = malloc(sizeof(long double) * 4 * M * M);
  AUX2 = malloc(sizeof(long double) * 4 * M * M);
  AUX_INV = malloc(sizeof(long double) * 4 * M * M);
  RAUX = malloc(sizeof(long double) * 4 * M * M);
  AUX3 = malloc(sizeof(long double) * 4 * M * M);
  S0 = malloc(sizeof(long double) * 2 * M);
  S1 = malloc(sizeof(long double) * 2 * M);
  
  if(XA == NULL || XE == NULL || YA == NULL || YE == NULL || XB == NULL || YB == NULL || H == NULL || XE_INV == NULL || YE_INV == NULL || M1 == NULL || M2 == NULL || M3 == NULL || M4 == NULL || TEMP0 == NULL || TEMP1 == NULL || AUX == NULL || AUX2 == NULL || AUX_INV == NULL || RAUX == NULL || AUX3 == NULL || S0 == NULL || S1 == NULL){

    if(XA != NULL) free(XA); if(XE != NULL) free(XE); if(YA != NULL) free(YA);
    if(YE != NULL) free(YE); if(H != NULL) free(H); if(XB != NULL) free(XB);
    if(YB != NULL) free(YB);

    if(XE_INV != NULL) free(XE_INV); if(YE_INV != NULL) free(YE_INV);

    if(M1 != NULL) free(M1); if(M2 != NULL) free(M2); if(M3 != NULL) free(M3); 
    if(M4 != NULL) free(M4); if(TEMP0 != NULL) free(TEMP0); if(TEMP1 != NULL) free(TEMP1);

    if(AUX != NULL) free(AUX); if(AUX2 != NULL) free(AUX2); if(AUX_INV != NULL) free(AUX_INV);
    if(RAUX != NULL) free(RAUX); if(AUX3 != NULL) free(AUX3); if(S0 != NULL) free(S0);
    if(S1 != NULL) free(S1);
    
    return 3;
  }

  for (int ry = 0; ry < nyr; ry++){
    for (int rx = 0; rx < nxr; rx++){

      // AUXILIARY VARIABLES
      long double lenx = XDOM[rx * 2], ntcx = XDOM[rx * 2 + 1], hx;
      hx = lenx / ntcx;
      long double leny = YDOM[ry * 2], ntcy = YDOM[ry * 2 + 1], hy;
      hy = leny / ntcy;
      int z = ZMAP[nxr * ry + rx];
      long double st = ZON[z * 2], ss = ZON[z * 2 + 1], c0;
      c0 = ss / st;
      long double Q = QMAP[nxr * ry + rx];

      for (int m = 0; m < M; m++){

        // AVERAGE FLUX MATRICES
        (*SM)[nyr * (M*rx + m) + ry] = Q * (1 + c0 / (1 - c0)) / st;
        long double k_miu, k_theta, k_w;
        for (int k = 0; k < M; k++){
          k_miu = MIU[k]; k_theta = THETA[k]; k_w = W[k];
          (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry] = 0.25 * c0 * k_w * k_miu / (st * hx * (1 - c0));
          (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry] = 0.25 * c0 * k_w * k_theta / (st * hy * (1 - c0));
          if (m == k){
            (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry] = (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry]
                                                          + k_miu / (st * hx);
            (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry] = (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry]
                                                          + k_theta / (st * hy);
          }
          if (k < M / 4){
            (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry] = (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry];
            (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry] = (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry];
          }
          else if (k >= M / 4 && k < M / 2){
            (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry] = - (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry];
            (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry] = (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry];
          }
          else if (k >= M / 2 && k < 3 * M / 4){
            (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry] = - (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry];
            (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry] = - (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry];
          }
          else if (k >= 3 * M / 4 && k < M){
            (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry] = (*FM0)[nyr * (nxr * (M * m + k) + rx) + ry];
            (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry] = - (*FM1)[nyr * (nxr * (M * m + k) + rx) + ry];
          }
        }


        if (m < M / 4){
          for (int k = 0; k < M; k++){
            XA[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hx / XVALS[M * z + k]);
            XE[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hx / XVALS[M * z + k]);

            YA[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hy / YVALS[M * z + k]);
            YE[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hy / YVALS[M * z + k]);
          }
        }
        else if (m >= M / 4 && m < M / 2){
          for (int k = 0; k < M; k++){
            XA[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hx / XVALS[M * z + k]);
            XE[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hx / XVALS[M * z + k]);

            YA[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hy / YVALS[M * z + k]);
            YE[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hy / YVALS[M * z + k]);
          }
        }
        else if (m >= M / 2 && m < 3 * M / 4){
          for (int k = 0; k < M; k++){
            XA[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hx / XVALS[M * z + k]);
            XE[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hx / XVALS[M * z + k]);

            YA[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hy / YVALS[M * z + k]);
            YE[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hy / YVALS[M * z + k]);
          }
        }
        else if (m >= 3 * M / 4 && m < M){
          for (int k = 0; k < M; k++){
            XA[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hx / XVALS[M * z + k]);
            XE[M * m + k] = XVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hx / XVALS[M * z + k]);

            YA[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(- 0.5 * st * hy / YVALS[M * z + k]);
            YE[M * m + k] = YVECTS[nz * (M * k + m) + z] * exp(0.5 * st * hy / YVALS[M * z + k]);
          }
        }

        for (int k = 0; k < M; k++){
          long double miu = MIU[k], theta = THETA[k], w = W[k];
          if (k < M / 4){
            XB[M * m + k] = - 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m) XB[M * m + k] = XB[M * m + k] - theta / (st * hy);

            YB[M * m + k] = - 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m) YB[M * m + k] = YB[M * m + k] - miu / (st * hx);
          }
          else if (k >= M / 4 && k < M/2){
            XB[M * m + k] = - 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m) XB[M * m + k] = XB[M * m + k] - theta / (st * hy);

            YB[M * m + k] = 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m) YB[M * m + k] = YB[M * m + k] + miu / (st * hx);
          }
          else if (k >= M / 2 && k < 3 * M / 4){
            XB[M * m + k] = 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m) XB[M * m + k] = XB[M * m + k] + theta / (st * hy);

            YB[M * m + k] = 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m) YB[M * m + k] = YB[M * m + k] + miu / (st * hx);
          }
          else if (k >= 3 * M / 4 && k < M){
            XB[M * m + k] = 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m) XB[M * m + k] = XB[M * m + k] + theta / (st * hy);

            YB[M * m + k] = - 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m) YB[M * m + k] = YB[M * m + k] - miu / (st * hx);
          }
          
        }

        H[m] = Q / (st * (1 - c0));

      }

      // RESPONSE MATRIX
      vcond = inv(M, XE, &XE_INV);
      if (vcond != 0) {
        if(XA != NULL) free(XA); if(XE != NULL) free(XE); if(YA != NULL) free(YA);
        if(YE != NULL) free(YE); if(H != NULL) free(H); if(XB != NULL) free(XB);
        if(YB != NULL) free(YB);

        if(XE_INV != NULL) free(XE_INV); if(YE_INV != NULL) free(YE_INV);

        if(M1 != NULL) free(M1); if(M2 != NULL) free(M2); if(M3 != NULL) free(M3); 
        if(M4 != NULL) free(M4); if(TEMP0 != NULL) free(TEMP0); if(TEMP1 != NULL) free(TEMP1);

        if(AUX != NULL) free(AUX); if(AUX2 != NULL) free(AUX2); if(AUX_INV != NULL) free(AUX_INV);
        if(RAUX != NULL) free(RAUX); if(AUX3 != NULL) free(AUX3); if(S0 != NULL) free(S0);
        if(S1 != NULL) free(S1);

        return vcond;
      }

      vcond = inv(M, YE, &YE_INV);
      if (vcond != 0) {
        if(XA != NULL) free(XA); if(XE != NULL) free(XE); if(YA != NULL) free(YA);
        if(YE != NULL) free(YE); if(H != NULL) free(H); if(XB != NULL) free(XB);
        if(YB != NULL) free(YB);

        if(XE_INV != NULL) free(XE_INV); if(YE_INV != NULL) free(YE_INV);

        if(M1 != NULL) free(M1); if(M2 != NULL) free(M2); if(M3 != NULL) free(M3); 
        if(M4 != NULL) free(M4); if(TEMP0 != NULL) free(TEMP0); if(TEMP1 != NULL) free(TEMP1);

        if(AUX != NULL) free(AUX); if(AUX2 != NULL) free(AUX2); if(AUX_INV != NULL) free(AUX_INV);
        if(RAUX != NULL) free(RAUX); if(AUX3 != NULL) free(AUX3); if(S0 != NULL) free(S0);
        if(S1 != NULL) free(S1);

        return vcond;
      }

      M1 = eye(M, &M1);

      TEMP0 = matrix_mult2(M, XE_INV, XB, &TEMP0);
      TEMP1 = matrix_mult2(M, XA, TEMP0, &TEMP1);
      TEMP0 = matrix_neg(M, TEMP1, &TEMP0);
      TEMP1 = matrix_sum(M, XB, TEMP0, &TEMP1);
      M2 = matrix_neg(M, TEMP1, &M2);
      
      TEMP0 = matrix_mult2(M, YE_INV, YB, &TEMP0);
      TEMP1 = matrix_mult2(M, YA, TEMP0, &TEMP1);
      TEMP0 = matrix_neg(M, TEMP1, &TEMP0);
      TEMP1 = matrix_sum(M, YB, TEMP0, &TEMP1);
      M3 = matrix_neg(M, TEMP1, &M3);

      M4 = eye(M, &M4);

      AUX = matrix_concat(M, M1, M2, M3, M4, &AUX);
      

      M1 = matrix_mult2(M, XA, XE_INV, &M1);

      M4 = matrix_mult2(M, YA, YE_INV, &M4);

      AUX2 = matrix_concat(M, M1, M2, M3, M4, &AUX2);
      

      vcond = inv(2*M, AUX, &AUX_INV);
      if (vcond != 0) {
        if(XA != NULL) free(XA); if(XE != NULL) free(XE); if(YA != NULL) free(YA);
        if(YE != NULL) free(YE); if(H != NULL) free(H); if(XB != NULL) free(XB);
        if(YB != NULL) free(YB);

        if(XE_INV != NULL) free(XE_INV); if(YE_INV != NULL) free(YE_INV);

        if(M1 != NULL) free(M1); if(M2 != NULL) free(M2); if(M3 != NULL) free(M3); 
        if(M4 != NULL) free(M4); if(TEMP0 != NULL) free(TEMP0); if(TEMP1 != NULL) free(TEMP1);

        if(AUX != NULL) free(AUX); if(AUX2 != NULL) free(AUX2); if(AUX_INV != NULL) free(AUX_INV);
        if(RAUX != NULL) free(RAUX); if(AUX3 != NULL) free(AUX3); if(S0 != NULL) free(S0);
        if(S1 != NULL) free(S1);

        return vcond;
      }
      
      RAUX = matrix_mult2(2*M, AUX_INV, AUX2, &RAUX);

      for (int j = 0; j < 2 * M; j++){
        for (int i = 0; i < 2 * M; i++){
          (*RM)[nyr * (nxr * (2*M*j + i) + rx) + ry] = RAUX[2*M*j + i];
        }
      }

      // PARTICULAR VECTOR
      TEMP0 = matrix_mult2(M, XA, XE_INV, &TEMP0);
      TEMP1 = matrix_neg(M, TEMP0, &TEMP1);
      TEMP0 = eye(M, &TEMP0);
      M1 = matrix_sum(M, TEMP0, TEMP1, &M1);

      M2 = zeros(M, &M2);

      M3 = zeros(M, &M3);

      TEMP0 = matrix_mult2(M, YA, YE_INV, &TEMP0);
      TEMP1 = matrix_neg(M, TEMP0, &TEMP1);
      TEMP0 = eye(M, &TEMP0);
      M4 = matrix_sum(M, TEMP0, TEMP1, &M4);

      AUX3 = matrix_concat(M, M1, M2, M3, M4, &AUX3);

      S0 = vector_concat(M, H, H, &S0);

      S1 = matrix_mult1(2*M, AUX3, S0, &S1);

      S0 = matrix_mult1(2*M, AUX_INV, S1, &S0);
      
      for (int i = 0; i < 2 * M; i++){
        (*PV)[nyr * (2*M*rx + i) + ry] = S0[i];
      }

      
    }
  }

  // FREE MEMORY
  if(XA != NULL) free(XA); if(XE != NULL) free(XE); if(YA != NULL) free(YA);
  if(YE != NULL) free(YE); if(H != NULL) free(H); if(XB != NULL) free(XB);
  if(YB != NULL) free(YB);

  if(XE_INV != NULL) free(XE_INV); if(YE_INV != NULL) free(YE_INV);

  if(M1 != NULL) free(M1); if(M2 != NULL) free(M2); if(M3 != NULL) free(M3); 
  if(M4 != NULL) free(M4); if(TEMP0 != NULL) free(TEMP0); if(TEMP1 != NULL) free(TEMP1);

  if(AUX != NULL) free(AUX); if(AUX2 != NULL) free(AUX2); if(AUX_INV != NULL) free(AUX_INV);
  if(RAUX != NULL) free(RAUX); if(AUX3 != NULL) free(AUX3); if(S0 != NULL) free(S0);
  if(S1 != NULL) free(S1);

  return 0;

}


long double* get_RM(int M, int nyr, int nxr, int ry, int rx, long double RM[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * 4 * M * M);
  }

  for (int j = 0; j < 2 * M; j++){
    for (int i = 0; i < 2 * M; i++){

      (*OUTPUT)[2*M*j + i] = RM[nyr * (nxr * (2*M*j + i) + rx) + ry];

    }
  }

  return *OUTPUT;

}


long double* get_PV(int M, int nyr, int nxr, int ry, int rx, long double PV[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * 2 * M);
  }

  for (int i = 0; i < 2 * M; i++){

    (*OUTPUT)[i] = PV[nyr * (2*M*rx + i) + ry];

  }

  return *OUTPUT;

}


long double* get_FM0(int M, int nyr, int nxr, int ry, int rx, long double FM0[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){

      (*OUTPUT)[M*j + i] = FM0[nyr * (nxr * (M*j + i) + rx) + ry];

    }
  }

  return *OUTPUT;

}


long double* get_FM1(int M, int nyr, int nxr, int ry, int rx, long double FM1[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M * M);
  }

  for (int j = 0; j < M; j++){
    for (int i = 0; i < M; i++){

      (*OUTPUT)[M*j + i] = FM1[nyr * (nxr * (M*j + i) + rx) + ry];

    }
  }

  return *OUTPUT;

}

long double* get_SM(int M, int nyr, int nxr, int ry, int rx, long double SM[], long double **OUTPUT){

  if ((*OUTPUT) == NULL){
    *OUTPUT = malloc(sizeof(long double) * M);
  }

  for (int i = 0; i < M; i++){

    (*OUTPUT)[i] = SM[nyr * (M*rx + i) + ry];

  }

  return *OUTPUT;

}

///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////// RM_CN ITERATIVE SCHEME
int rm_cn (int N,               // Quadrature order
           int nz,              // Number of zones
           long double ZON[],   // Zone entries
           int nxr,             // Number of regions in X
           long double XDOM[],  // X Region entries
           int nyr,             // Number of regions in Y
           long double YDOM[],  // Y Region entries
           int ZMAP[],          // Zone map
           long double QMAP[],  // External source map
           long double BC[],    // Boundary conditions
           long double tol,     // Tolerance
           long double W[],     // Quadrature weight
           long double RM[],    // Response matrix
           long double PV[],    // Source vector
           long double FM0[],   // Average flux matrix 0
           long double FM1[],   // Average flux matrix 1
           long double SM[],    // Average source vector
           long double **MFLUX, // Escalar flux in the nodes
           long double **MFLOW, // Angular flux in the nodes
           long double **XFLOW, // Angular flux at the y edges
           long double **YFLOW, // Angular flux at the x edges
           int *ITER,            // Iteration
           long double *cpu_time // CPU time
           ){

  // INITIALIZATION
  int ntc_x = 0, ntc_y = 0;
  for (int rx = 0; rx < nxr; rx++) {
		ntc_x = ntc_x + (int)XDOM[2*rx + 1];
	}
	for (int ry = 0; ry < nyr; ry++) {
		ntc_y = ntc_y + (int)YDOM[2 * ry + 1];
	}
  int M = N * (N + 2) / 2;
  *MFLUX = malloc(sizeof(long double) * (ntc_y * ntc_x));
	*MFLOW = malloc(sizeof(long double) * (ntc_y * ntc_x) * M);
	*XFLOW = malloc(sizeof(long double) * (ntc_y * (ntc_x + 1)) * M);
	*YFLOW = malloc(sizeof(long double) * ((ntc_y + 1) * ntc_x) * M);
  if (*MFLUX == NULL || *MFLOW == NULL || *XFLOW == NULL || *YFLOW == NULL) return 3;
  for (int j = 0; j < ntc_y; j++) {
		for (int i = 0; i < ntc_x; i++) {
			(*MFLUX)[ntc_x * j + i] = 0.0;
			for (int m = 0; m < M; m++) {
				(*MFLOW)[M * (ntc_x * j + i) + m] = 0.0;
			}
		}
	}
  for (int j = 0; j < ntc_y; j++) {
		for (int i = 0; i < ntc_x + 1; i++) {
			for (int m = 0; m < M; m++) {
				(*XFLOW)[M * ((ntc_x + 1) * j + i) + m] = 0.0;
				if ((i == 0) && (BC[0] != -1.0)) {
					(*XFLOW)[M * ((ntc_x + 1) * j + i) + m] = BC[0];
				}
				if ((i == ntc_x) && (BC[2] != -1.0)) {
					(*XFLOW)[M * ((ntc_x + 1) * j + i) + m] = BC[2];
				}
			}
		}
	}
  for (int j = 0; j < ntc_y + 1; j++) {
		for (int i = 0; i < ntc_x; i++) {
			for (int m = 0; m < M; m++) {
				(*YFLOW)[M * (ntc_x * j + i) + m] = 0.0;
				if (j == 0 && BC[1] != -1.0) {
					(*YFLOW)[M * (ntc_x * j + i) + m] = BC[1];
				}
				if (j == ntc_y && BC[3] != -1.0) {
					(*YFLOW)[M * (ntc_x * j + i) + m] = BC[3];
				}
			}
		}
	}

  // VARIABLES
  clock_t start, end;
  long double ERR = 1.0;
  int j_b, j_f, i_b, i_f;
  int nc_y, nc_x;
  long double *IN = NULL, *OUT = NULL; 
  long double *RESP_MATRIX = NULL, *PART_VECTOR = NULL, *V_AUX = NULL;
  IN = malloc(sizeof(long double) * 2 * M);
  OUT = malloc(sizeof(long double) * 2 * M);
  RESP_MATRIX = malloc(sizeof(long double) * 4 * M * M);
  PART_VECTOR = malloc(sizeof(long double) * 2 * M);
  V_AUX = malloc(sizeof(long double) * 2 * M);
  long double *X_VECTOR = NULL, *Y_VECTOR = NULL, *M_VECTOR = NULL, *V_AUX2 = NULL;
  long double *FM_0 = NULL, *FM_1 = NULL, *S_VECTOR = NULL;
  X_VECTOR = malloc(sizeof(long double) * M);
  Y_VECTOR = malloc(sizeof(long double) * M);
  M_VECTOR = malloc(sizeof(long double) * M);
  V_AUX2 = malloc(sizeof(long double) * M);
  FM_0 = malloc(sizeof(long double) * M * M);
  FM_1 = malloc(sizeof(long double) * M * M);
  S_VECTOR = malloc(sizeof(long double) * M);
  long double m_sum, w, flux0, flux;
  
  if (IN == NULL || OUT == NULL || RESP_MATRIX == NULL || PART_VECTOR == NULL || V_AUX == NULL || X_VECTOR == NULL || Y_VECTOR == NULL || M_VECTOR == NULL || V_AUX2 == NULL || FM_0 == NULL || FM_1 == NULL || S_VECTOR == NULL){
    if (IN != NULL) free(IN); if (OUT != NULL) free(OUT); if (RESP_MATRIX != NULL) free(RESP_MATRIX);
    if (PART_VECTOR != NULL) free(PART_VECTOR); if (V_AUX != NULL) free(V_AUX);

    if (X_VECTOR != NULL) free(X_VECTOR); if (Y_VECTOR != NULL) free(Y_VECTOR);
    if (M_VECTOR != NULL) free(M_VECTOR); if (V_AUX2 != NULL) free(V_AUX2);
    if (FM_0 != NULL) free(FM_0); if (FM_1 != NULL) free(FM_1); if (S_VECTOR != NULL) free(S_VECTOR);

    return 3;
  }

  // ITERATIVE PROCESS
  printf("\n3. ITERATIVE PROCESS:\n\n");
	printf("ITER\tMAX_ERR (MFLUX)\n");
  start = clock(); *ITER = -1;
  sleep(1);
  while (ERR > tol && *ITER < 10000){
    ERR = 0.0; *ITER = *ITER + 1;

    // 1. SW - > NE SWEEP
    j_b = 0; j_f = j_b + 1;
    for (int ry = 0; ry < nyr; ry++) {
      nc_y = (int)YDOM[2 * ry + 1]; 
      for (int j = 0; j < nc_y; j++) {
        i_b = 0; i_f = i_b + 1;
        for (int rx = 0; rx < nxr; rx++) {
          RESP_MATRIX = get_RM(M, nyr, nxr, ry, rx, RM, &RESP_MATRIX);
          PART_VECTOR = get_PV(M, nyr, nxr, ry, rx, PV, &PART_VECTOR);
          nc_x = (int)XDOM[2 * rx + 1];
          for (int i = 0; i < nc_x; i++) {

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 4 && m < M / 2){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
            }

            V_AUX = matrix_mult1(2*M, RESP_MATRIX, IN, &V_AUX);
            OUT = vector_sum(2*M, V_AUX, PART_VECTOR, &OUT);

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 4 && m < M / 2){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
            }

            i_b = i_b + 1; i_f = i_b + 1;

          }
        }
        j_b = j_b + 1; j_f = j_b + 1;
      }
    }

    // 2. SE - > NW SWEEP
    j_b = 0; j_f = j_b + 1;
    for (int ry = 0; ry < nyr; ry++) {
      nc_y = (int)YDOM[2 * ry + 1]; 
      for (int j = 0; j < nc_y; j++) {
        i_b = ntc_x - 1; i_f = i_b + 1;
        for (int rx = nxr-1; rx >= 0; rx--) {
          RESP_MATRIX = get_RM(M, nyr, nxr, ry, rx, RM, &RESP_MATRIX);
          PART_VECTOR = get_PV(M, nyr, nxr, ry, rx, PV, &PART_VECTOR);
          nc_x = (int)XDOM[2 * rx + 1];
          for (int i = 0; i < nc_x; i++) {

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 4 && m < M / 2){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
            }

            V_AUX = matrix_mult1(2*M, RESP_MATRIX, IN, &V_AUX);
            OUT = vector_sum(2*M, V_AUX, PART_VECTOR, &OUT);

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 4 && m < M / 2){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
            }

            i_b = i_b - 1; i_f = i_b + 1;

          }
        }
        j_b = j_b + 1; j_f = j_b + 1;
      }
    }

    // 3. NE - > SW SWEEP
    j_b = ntc_y - 1; j_f = j_b + 1;
    for (int ry = nyr - 1; ry >= 0; ry--) {
      nc_y = (int)YDOM[2 * ry + 1]; 
      for (int j = 0; j < nc_y; j++) {
        i_b = ntc_x - 1; i_f = i_b + 1;
        for (int rx = nxr - 1; rx >= 0; rx--) {
          RESP_MATRIX = get_RM(M, nyr, nxr, ry, rx, RM, &RESP_MATRIX);
          PART_VECTOR = get_PV(M, nyr, nxr, ry, rx, PV, &PART_VECTOR);
          nc_x = (int)XDOM[2 * rx + 1];
          for (int i = 0; i < nc_x; i++) {

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 4 && m < M / 2){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
            }

            V_AUX = matrix_mult1(2*M, RESP_MATRIX, IN, &V_AUX);
            OUT = vector_sum(2*M, V_AUX, PART_VECTOR, &OUT);

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 4 && m < M / 2){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
            }

            i_b = i_b - 1; i_f = i_b + 1;

          }
        }
        j_b = j_b - 1; j_f = j_b + 1;
      }
    }

    // 4. NW - > SE SWEEP
    j_b = ntc_y - 1; j_f = j_b + 1;
    for (int ry = nyr - 1; ry >= 0; ry--) {
      nc_y = (int)YDOM[2 * ry + 1]; 
      for (int j = 0; j < nc_y; j++) {
        i_b = 0; i_f = i_b + 1;
        for (int rx = 0; rx < nxr; rx++) {
          RESP_MATRIX = get_RM(M, nyr, nxr, ry, rx, RM, &RESP_MATRIX);
          PART_VECTOR = get_PV(M, nyr, nxr, ry, rx, PV, &PART_VECTOR);
          nc_x = (int)XDOM[2 * rx + 1];
          for (int i = 0; i < nc_x; i++) {

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 4 && m < M / 2){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                IN[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                IN[M + m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
            }

            V_AUX = matrix_mult1(2*M, RESP_MATRIX, IN, &V_AUX);
            OUT = vector_sum(2*M, V_AUX, PART_VECTOR, &OUT);

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 4 && m < M / 2){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_f + i_b) + m] = OUT[M + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m] = OUT[m];
                (*YFLOW)[M * (ntc_x * j_b + i_b) + m] = OUT[M + m];
              }
            }

            i_b = i_b + 1; i_f = i_b + 1;

          }
        }
        j_b = j_b - 1; j_f = j_b + 1;
      }
    }

    // REFLECTIVE BOUNDARY CONDITIONS
		j_b = 0;
		for (int ry = 0; ry < nyr; ry++) {
			nc_y = (int)YDOM[2 * ry + 1];
			for (int j = 0; j < nc_y; j++) {
				for (int m = 0; m < M/4; m++) {
					if (BC[0] == -1.0) {
						(*XFLOW)[M * ((ntc_x + 1) * j_b + 0) + m] = (*XFLOW)[M * ((ntc_x + 1) * j_b + 0) + M/4 + m];
						(*XFLOW)[M * ((ntc_x + 1) * j_b + 0) + 3*M/4 + m] = (*XFLOW)[M * ((ntc_x + 1) * j_b + 0) + M/2 + m];
					}
					if (BC[2] == -1.0) {
						(*XFLOW)[M * ((ntc_x + 1) * j_b + ntc_x) + M/4 + m] = (*XFLOW)[M * ((ntc_x + 1) * j_b + ntc_x) + m];
						(*XFLOW)[M * ((ntc_x + 1) * j_b + ntc_x) + M/2 + m] = (*XFLOW)[M * ((ntc_x + 1) * j_b + ntc_x) + 3*M/4 + m];
					}
				}
				j_b = j_b + 1;
			}
		}
		i_b = 0;
		for (int rx = 0; rx < nxr; rx++) {
			nc_x = (int)XDOM[2 * rx + 1];
			for (int i = 0; i < nc_x; i++) {
				for (int m = 0; m < M / 4; m++) {
					if (BC[1] == -1.0) {
						(*YFLOW)[M * (ntc_x * 0 + i_b) + m] = (*YFLOW)[M * (ntc_x * 0 + i_b) + 3*M/4 + m];
						(*YFLOW)[M * (ntc_x * 0 + i_b) + M/4 + m] = (*YFLOW)[M * (ntc_x * 0 + i_b) + M/2 + m];
					}
					if (BC[3] == -1.0) {
						(*YFLOW)[M * (ntc_x * ntc_y + i_b) + 3*M/4 + m] = (*YFLOW)[M * (ntc_x * ntc_y + i_b) + m];
						(*YFLOW)[M * (ntc_x * ntc_y + i_b) + M/2 + m] = (*YFLOW)[M * (ntc_x * ntc_y + i_b) + M/4 + m];
					}
				}
				i_b = i_b + 1;
			}
		}

    // SCALAR FLUX, SOURCE TERM AND STOP CRITERIA
    i_b = 0; i_f = i_b + 1;
    for (int rx = 0; rx < nxr; rx++){
      nc_x = (int)XDOM[2 * rx + 1];
			for (int i = 0; i < nc_x; i++) {
        j_b = 0; j_f = j_b + 1;
        for (int ry = 0; ry < nyr; ry++){
          nc_y = (int)YDOM[2 * ry + 1];
          FM_0 = get_FM0(M, nyr, nxr, ry, rx, FM0, &FM_0);
          FM_1 = get_FM1(M, nyr, nxr, ry, rx, FM1, &FM_1);
          S_VECTOR = get_SM(M, nyr, nxr, ry, rx, SM, &S_VECTOR);
          for (int j = 0; j < nc_y; j++){

            for (int m = 0; m < M; m++){
              if (m < M / 4){
                X_VECTOR[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m]
                              - (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                Y_VECTOR[m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m]
                              - (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 4 && m < M / 2){
                X_VECTOR[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m]
                              - (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                Y_VECTOR[m] = (*YFLOW)[M * (ntc_x * j_f + i_b) + m]
                              - (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
              }
              else if (m >= M / 2 && m < 3 * M / 4){
                X_VECTOR[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m]
                              - (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m];
                Y_VECTOR[m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m]
                              - (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
              else if (m >= 3 * M / 4 && m < M){
                X_VECTOR[m] = (*XFLOW)[M * ((ntc_x + 1)* j_b + i_f) + m]
                              - (*XFLOW)[M * ((ntc_x + 1)* j_b + i_b) + m];
                Y_VECTOR[m] = (*YFLOW)[M * (ntc_x * j_b + i_b) + m]
                              - (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
              }
            }

            V_AUX2 = matrix_mult1(M, FM_0, X_VECTOR, &V_AUX2);
            X_VECTOR = matrix_mult1(M, FM_1, Y_VECTOR, &X_VECTOR);
            Y_VECTOR = vector_sum(M, V_AUX2, X_VECTOR, &Y_VECTOR);
            V_AUX2 = vector_neg(M, Y_VECTOR, &V_AUX2);
            M_VECTOR = vector_sum(M, V_AUX2, S_VECTOR, &M_VECTOR);

            m_sum = 0;
            for (int m = 0; m < M; m++){
              w = W[m];
              (*MFLOW)[M * (ntc_x * j_b + i_b) + m] = M_VECTOR[m];
              m_sum = m_sum + M_VECTOR[m] * w;
            }
            flux = 0.25 * m_sum;

            flux0 = (*MFLUX)[ntc_x * j_b + i_b];
            (*MFLUX)[ntc_x * j_b + i_b] = flux;

            if (fabs(1 - flux0 / flux) > ERR) ERR = fabs(1 - flux0 / flux);

            j_b = j_b + 1; j_f = j_b + 1;
          }
        }

        i_b = i_b + 1; i_f = i_b + 1;
      }
    }

    printf("%d\t%.5Le\n", *ITER, ERR);
  }
  end = clock();
  *cpu_time = (long double)(end - start) / CLOCKS_PER_SEC - 1.0;
  printf("CPU TIME = %.5Le\n\n", *cpu_time);

  free(IN); free(OUT); 
  free(RESP_MATRIX); free(PART_VECTOR); free(V_AUX);
  free(X_VECTOR); free(Y_VECTOR); free(M_VECTOR); free(V_AUX2);
  free(FM_0); free(FM_1); free(S_VECTOR);

  if (*ITER >= 10000) return 5;
  return 0;

}

///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////// AUXILIARY FUNCTIONS
void print_problem(int N,           // Quadrature order
                   int nz,          // Number of zones
                   long double ZON[],    // Zone entries
                   int nxr,         // Number of regions in X
                   long double XDOM[],      // X Region entries
                   int nyr,         // Number of regions in Y
                   long double YDOM[],      // Y Region entries
                   int ZMAP[],      // Zone map
                   long double QMAP[],   // External source map
                   long double BC[],     // Boundary conditions
                   long double tol      // Tolerance
                  ){

  printf("\n1. PROBLEM SETUP:\n\n");

  // QUADRATURE ORDER
  printf("N = %d\n\n", N);

  // NUMBER OF ZONES
  printf("NZ = %d\n\n", nz);

  // ZONE ENTRIES
	printf("ZON:\n");
  printf("NO.\tS_T\tS_S\n");
  for (int z = 0; z < nz; z++){
      printf("%d\t%.4Lf\t%.4Lf\n", z+1, ZON[z*2], ZON[z*2 + 1]);
  }
  printf("\n");

  // NUMBER OF REGIONS IN X
  printf("NR_X = %d\n\n", nxr);

  // REGION ENTRIES IN X
	printf("XDOM:\n");
  printf("No.\tLEN\tNODES\n");
  for (int xr = 0; xr < nxr; xr++){
      printf("%d\t%.2Lf\t%d\n", xr + 1, XDOM[xr * 2], (int)XDOM[xr * 2 + 1]);
  }
  printf("\n");

  // NUMBER OF REGIONS IN Y
  printf("NR_Y = %d\n\n", nyr);

  // REGION ENTRIES IN Y
	printf("YDOM:\n");
  printf("NO.\tLEN\tNODES\n");
  for (int yr = 0; yr < nyr; yr++){
      printf("%d\t%.2Lf\t%d\n", yr + 1, YDOM[yr * 2], (int)YDOM[yr * 2 + 1]);
  }
  printf("\n");

  // ZONE MAPPING
  printf("ZMAP:\n");
  for (int yr = 0; yr < nyr; yr++) {
		for (int xr = 0; xr < nxr; xr++) {
			printf("%d\t", ZMAP[yr * nxr + xr] + 1);
		}
		printf("\n");
  }
  printf("\n");

  // EXTERNAL SOURCE MAPPING
  printf("Q_MAP:\n");
  for (int yr = 0; yr < nyr; yr++) {
		for (int xr = 0; xr < nxr; xr++) {
			printf("%.4Lf\t", QMAP[yr * nxr + xr]);
		}
		printf("\n");
  }
  printf("\n");

  // BOUNDARY CONDITIONS
	printf("BOUNDARY_CONDITIONS =\n");
	for (int c = 0; c < 4; c++) {
		printf("%d. ", c + 1);
		if (BC[c] == 0.0) printf("Vacuum\n");
		else if (BC[c] == -1.0) printf("Reflective\n");
    else printf("Isotropic %.2f\n", BC[c]);
	}
	printf("\n");

  // TOLERANCE
	printf("TOL = %.4Le\n\n", tol);

}


void post_processing(int N,               // Quadrature order
                     int nz,              // Number of zones
                     long double ZON[],   // Zone entries
                     int nxr,             // Number of regions in X
                     long double XDOM[],  // X Region entries
                     int nyr,             // Number of regions in Y
                     long double YDOM[],  // Y Region entries
                     int ZMAP[],          // Zone map
                     long double QMAP[],  // External source map
                     long double BC[],    // Boundary conditions
                     long double MIU[],   // Ordinates in X
                     long double THETA[], // Ordinates in Y
                     long double W[],
	                   long double MFLUX[], // Scalar flux in the nodes
                     long double MFLOW[], // Angular flux in the nodes
                     long double XFLOW[], // Angular flux at the y edges
                     long double YFLOW[]  // Angular flux at the x edges
					 ){
	
	// INICIALIZATION
	int j_b = 0, i_b = 0, j_f, i_f, ntc_x = 0, ntc_y = 0;;
	int nc_y, nc_x, z;
	long double h_y, h_x, area_r, sa, miu, theta, w, len_y, len_x;
	long double *MFLUX_R = NULL, *ABS_R = NULL, FUGA[4] = {0.0};
	MFLUX_R = malloc(sizeof(long double)*nyr*nxr); assert(MFLUX_R != NULL);
	ABS_R = malloc(sizeof(long double)*nyr*nxr); assert(ABS_R != NULL);
	for (int ry = 0; ry < nyr; ry++) {
		for (int rx = 0; rx < nxr; rx++) {
			MFLUX_R[nxr*ry + rx] = 0.0;
			ABS_R[nxr*ry + rx] = 0.0;
		}
	}
  for (int rx = 0; rx < nxr; rx++) {
		ntc_x = ntc_x + (int)XDOM[2*rx + 1];
	}
	for (int ry = 0; ry < nyr; ry++) {
		ntc_y = ntc_y + (int)YDOM[2 * ry + 1];
	}
  int M = N * (N + 2) / 2;

	// CALCULATION OF MFLUX and ABS_R
  j_b = 0;
	for(int ry = 0; ry < nyr; ry++){
		len_y = YDOM[2*ry]; nc_y = (int)YDOM[2*ry + 1]; h_y = (long double)len_y/nc_y;
		for(int j = 0; j < nc_y; j++){
			i_b = 0;
			for(int rx = 0; rx < nxr; rx++){
				len_x = XDOM[2*rx]; nc_x = (int)XDOM[2*rx + 1]; h_x = (long double)len_x/nc_x; 
				area_r = len_y*len_x; z = ZMAP[nxr*ry + rx];
				sa = ZON[2*z] - ZON[2*z + 1];
				for(int i = 0; i < nc_x; i++){
					MFLUX_R[nxr*ry + rx] = MFLUX_R[nxr*ry + rx] + h_y*h_x*MFLUX[ntc_x*j_b + i_b]/area_r;
					ABS_R[nxr*ry + rx] = ABS_R[nxr*ry + rx] + sa*h_y*h_x*MFLUX[ntc_x*j_b + i_b];
					i_b = i_b + 1;
				}
			}
			j_b = j_b + 1;
		}
	}

	// CALCULATION OF THE LEAKAGE FROM THE LEFT AND THE RIGHT BOUNDARIES
	j_b = 0;
	for(int ry = 0; ry < nyr; ry++){
		len_y = YDOM[2*ry]; nc_y = (int)YDOM[2*ry + 1]; h_y = (long double)len_y/nc_y;
		for(int j = 0; j < nc_y; j++){
			// LEFT
			for (int m = M/4; m < 3*M/4; m++){
				miu = MIU[m]; w = W[m];	                          
				FUGA[0] = FUGA[0] + 0.5*h_y*w*miu*XFLOW[M * ((ntc_x + 1) * j_b + 0) + m];
			}
			// RIGHT
			for (int m = 0; m < M/4; m++){
				miu = MIU[m]; w = W[m];		                          
				FUGA[2] = FUGA[2] + 0.5*h_y*w*miu*XFLOW[M * ((ntc_x + 1) * j_b + ntc_x) + m];
			}
			for (int m = 3*M/4; m < M; m++){
				miu = MIU[m]; w = W[m];		                          
				FUGA[2] = FUGA[2] + 0.5*h_y*w*miu*XFLOW[M * ((ntc_x + 1) * j_b + ntc_x) + m];
			}                         
			j_b = j_b + 1;
		}
	}

	// CALCULATION OF THE LEAKAGE FROM THE BOTTOM AND THE TOP BOUNDARIES
	i_b = 0;
	for(int rx = 0; rx < nxr; rx++){
		len_x = XDOM[2*rx]; nc_x = (int)XDOM[2*rx + 1]; h_x = (long double)len_x/nc_x;
		for(int i = 0; i < nc_x; i++){
			for (int m = M/2; m < M; m++){
				theta = THETA[m]; w = W[m];   	
				FUGA[1] = FUGA[1] + 0.5*h_x*w*theta*YFLOW[M * (ntc_x * 0 + i_b) + m];
			}
			for (int m = 0; m < M/2; m++){
				theta = THETA[m]; w = W[m];
				FUGA[3] = FUGA[3] + 0.5*h_x*w*theta*YFLOW[M * (ntc_x * ntc_y + i_b) + m];
			}
			i_b = i_b + 1;
		}
	}

	// PRINT POST PROCESSING
	printf("\n\n4. POST - PROCESSING:\n\n");

	// SCALAR FLUX PER REGION
	printf("SCALAR FLUX PER REGION =\n");
	for (int rx = 0; rx < nxr; rx++) {
		printf("\tRX_%d\t", rx + 1);
	}
	printf("\n");
	for (int ry = 0; ry < nyr; ry++) {
		printf("RY_%d\t", ry + 1);
		for (int rx = 0; rx < nxr; rx++) {
			printf("% .5Le\t", MFLUX_R[nxr * ry + rx]);
		}
		printf("\n");
	}
	printf("\n");

	// ABSORPTION RATE PER REGION
	printf("ABSORPTION RATE PER REGION =\n");
	for (int rx = 0; rx < nxr; rx++) {
		printf("\tRX_%d\t", rx + 1);
	}
	printf("\n");
	for (int ry = 0; ry < nyr; ry++) {
		printf("RY_%d\t", ry + 1);
		for (int rx = 0; rx < nxr; rx++) {
			printf("% .5Le\t", ABS_R[nxr * ry + rx]);
		}
		printf("\n");
	}
	printf("\n");

	// LEAKAGE
  printf("LEAKAGE =\n");
	printf("LEFT\t\tBOTTOM\t\tRIGHT\t\tTOP\n");
  for(int i = 0; i < 4; i++){
    if (BC[i] == -1.0) printf("-\t\t");
		else printf("%.4Le\t", FUGA[i]);
	}
	printf("\n\n");
	
	free(MFLUX_R); free(ABS_R);

}


void json_output(int N,               // Quadrature order
                 int nxr,             // Number of regions in X
                 long double XDOM[],  // X Region entries
                 int nyr,             // Number of regions in Y
                 long double YDOM[],  // Y Region entries
                 int status,          // Status
                 int ITER,            // Iterations
                 long double cpu_time,// CPU time
	               long double MFLUX[], // Scalar flux in the nodes
                 long double MFLOW[], // Angular flux in the nodes
                 long double XFLOW[], // Angular flux at the y edges
                 long double YFLOW[]  // Angular flux at the x edges
					 ){
	
	// INICIALIZATION
	int j_b = 0, i_b = 0, j_f, i_f, ntc_x = 0, ntc_y = 0;;
	int nc_y, nc_x;
  for (int rx = 0; rx < nxr; rx++) {
		ntc_x = ntc_x + (int)XDOM[2*rx + 1];
	}
	for (int ry = 0; ry < nyr; ry++) {
		ntc_y = ntc_y + (int)YDOM[2 * ry + 1];
	}
  int M = N * (N + 2) / 2;

  if (status != 0) {
    printf("{\n\"STATUS\": %d\n}\n", status);
  }
  else {

    // MFLUX VERIFICATION
    j_b = 0;
    for(int ry = 0; ry < nyr; ry++){
      nc_y = (int)YDOM[2*ry + 1];
      for(int j = 0; j < nc_y; j++){
        i_b = 0;
        for(int rx = 0; rx < nxr; rx++){
          nc_x = (int)XDOM[2*rx + 1];
          for(int i = 0; i < nc_x; i++){
            if (isfinite(MFLUX[ntc_x*j_b + i_b]) == 0){
              status = 6;
              break;
            }
          }
          if (status != 0) break;
        }
        if (status != 0) break;
        j_b = j_b + 1;
      }
      if (status != 0) break;
    }

    // STATUS
    if (status != 0) {
      printf("{\n\"STATUS\": %d\n}\n", status);
    }
    else {

      printf("{\n\"STATUS\": %d", status);

      // ITERATIONS
      printf(",\n\"ITER\": %d,\n", ITER);

      // CPU TIME
      printf("\"CPU\": %.10Le,\n", cpu_time);

      // MFLUX
      printf("\"MFLUX\": [\n");
      j_b = 0;
      for(int ry = 0; ry < nyr; ry++){
        nc_y = (int)YDOM[2*ry + 1];
        for(int j = 0; j < nc_y; j++){
          i_b = 0;
          printf("[");
          for(int rx = 0; rx < nxr; rx++){
            nc_x = (int)XDOM[2*rx + 1];
            for(int i = 0; i < nc_x; i++){
              if (i_b == ntc_x - 1) {
                if (j_b == ntc_y - 1) printf(" %.10Le ]\n", MFLUX[ntc_x*j_b + i_b]);
                else printf(" %.10Le ],\n", MFLUX[ntc_x*j_b + i_b]);
              }
              else printf(" %.10Le,", MFLUX[ntc_x*j_b + i_b]);
              i_b = i_b + 1;
            }
          }
          j_b = j_b + 1;
        }
      }
      printf("],\n");

      // XFLOW
      printf("\"XFLOW\": [\n");
      for (int m = 0; m < M; m++){
        j_b = 0;
        printf("[\n");
        for(int ry = 0; ry < nyr; ry++){
          nc_y = (int)YDOM[2*ry + 1];
          for(int j = 0; j < nc_y; j++){
            i_b = 0;
            printf("[");
            for(int rx = 0; rx < nxr; rx++){
              nc_x = (int)XDOM[2*rx + 1];
              for(int i = 0; i < nc_x; i++){
                printf(" %.10Le,", XFLOW[M * ((ntc_x + 1) * j_b + i_b) + m]);
                i_b = i_b + 1;
              }
            }
            if (i_b == ntc_x) {
              if (j_b == ntc_y - 1) printf(" %.10Le ]\n", XFLOW[M * ((ntc_x + 1) * j_b + i_b) + m]);
              else printf(" %.10Le ],\n", XFLOW[M * ((ntc_x + 1) * j_b + i_b) + m]);
            }
            j_b = j_b + 1;
          }
        }
        if (m == M-1) printf("]\n");
        else printf("],\n");
      }
      printf("],\n");

      // YFLOW
      printf("\"YFLOW\": [\n");
      for (int m = 0; m < M; m++){
        j_b = 0;
        printf("[\n");
        for(int ry = 0; ry < nyr; ry++){
          nc_y = (int)YDOM[2*ry + 1];
          for(int j = 0; j < nc_y; j++){
            i_b = 0;
            printf("[");
            for(int rx = 0; rx < nxr; rx++){
              nc_x = (int)XDOM[2*rx + 1];
              for(int i = 0; i < nc_x; i++){
                if (i_b == ntc_x - 1) printf(" %.10Le ],\n", YFLOW[M * (ntc_x * j_b + i_b) + m]);
                else printf(" %.10Le,", YFLOW[M * (ntc_x * j_b + i_b) + m]);
                i_b = i_b + 1;
              }
            }
            j_b = j_b + 1;
          }
        }
        if (j_b == ntc_y){
          i_b = 0;
          printf("[");
          for(int rx = 0; rx < nxr; rx++){
            nc_x = (int)XDOM[2*rx + 1];
            for(int i = 0; i < nc_x; i++){
              if (i_b == ntc_x - 1) printf(" %.10Le ]\n", YFLOW[M * (ntc_x * j_b + i_b) + m]);
              else printf(" %.10Le,", YFLOW[M * (ntc_x * j_b + i_b) + m]);
              i_b = i_b + 1;
            }
          }
        }
        if (m == M-1) printf("]\n");
        else printf("],\n");
      }
      printf("]\n");

      printf("}\n");
    }
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////