#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#define printf __mingw_printf
#define eps 1e-12

//////////////////////////////////////////////////////////////////////////////////// PROTOTYPES

///////////////////////////////////////////////////////// FREE MEMORY REQUESTED
void free_memory(long double **ZON, long double **XDOM, long double **YDOM, int **ZMAP, long double **QMAP, long double **MIU, long double **THETA, long double **FI, long double **W, long double **LIST, long double **MFLUX, long double **MFLOW, long double **XFLOW, long double **YFLOW);

//////////////////////////////////////////////////////// LOAD PROBLEM FROM FILE
int input_by_txt(int *N, int *nz, long double **ZON, int *nxr, long double **XDOM, int *nyr, long double **YDOM, int **ZMAP, long double **QMAP, long double *BC, long double *tol, const char *filename);

/////////////////////////////////////////////////////////// GENERATE QUADRATURE
int quad(int N, long double **MIU, long double **THETA, long double **CHI, long double **W, long double **LIST);

/////////////////////////////////////////////////////////// DD ITERATIVE SCHEME
int dd (int N, int nz, long double ZON[], int nxr, long double XDOM[], int nyr, long double YDOM[], int ZMAP[], long double QMAP[], long double BC[], long double tol, long double MIU[], long double THETA[], long double W[], long double **MFLUX, long double **MFLOW, long double **XFLOW, long double **YFLOW, int *ITER, long double *cpu_time);

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

  // ITERATIVE SCHEME
  status = dd (N, nz, ZON, nxr, XDOM, nyr, YDOM, ZMAP, QMAP, BC, tol, MIU, THETA, W, &MFLUX, &MFLOW, &XFLOW, &YFLOW, &ITER, &cpu_time);
  if (status != 0){
    json_output(N, nxr, XDOM, nyr, YDOM, status, ITER, cpu_time, MFLUX, MFLOW, XFLOW, YFLOW);
    free_memory(&ZON, &XDOM, &YDOM, &ZMAP, &QMAP,
                &MIU, &THETA, &FI, &W, &LIST,
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
              &MFLUX, &MFLOW, &XFLOW, &YFLOW);

  return 0;

}
///////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////// IMPLEMENTATION

///////////////////////////////////////////////////////// FREE MEMORY REQUESTED
void free_memory(long double **ZON, long double **XDOM, long double **YDOM, int **ZMAP, long double **QMAP, long double **MIU, long double **THETA, long double **FI, long double **W, long double **LIST, long double **MFLUX, long double **MFLOW, long double **XFLOW, long double **YFLOW){
  
  // INPUT VARIABLES
  if (*ZON != NULL) free(*ZON); if (*XDOM != NULL) free(*XDOM); if (*YDOM != NULL) free(*YDOM); 
  if (*ZMAP != NULL) free(*ZMAP); if (*QMAP != NULL) free(*QMAP);

  // QUADRATURE VARIABLES
  if (*MIU != NULL) free(*MIU); if (*THETA != NULL) free(*THETA); if (*FI != NULL) free(*FI);
  if (*W != NULL) free(*W); if (*LIST != NULL) free(*LIST);
  
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
                long double **XDOM,       // X Region entries
                int *nyr,                 // Number of regions in Y
                long double **YDOM,       // Y Region entries
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

//////////////////////////////////////////////////////// RM_CN ITERATIVE SCHEME
int dd (int N,               // Quadrature order
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
        long double MIU[],   // Ordinate in X
        long double THETA[], // Ordinate in Y
        long double W[],     // Quadrature weight
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
  long double *S = NULL;
  S = malloc(sizeof(long double) * (ntc_y * ntc_x));
  *MFLUX = malloc(sizeof(long double) * (ntc_y * ntc_x));
	*MFLOW = malloc(sizeof(long double) * (ntc_y * ntc_x) * M);
	*XFLOW = malloc(sizeof(long double) * (ntc_y * (ntc_x + 1)) * M);
	*YFLOW = malloc(sizeof(long double) * ((ntc_y + 1) * ntc_x) * M);
  if (*MFLUX == NULL || *MFLOW == NULL || *XFLOW == NULL || *YFLOW == NULL || S == NULL){
    if (S != NULL) free(S); return 3;
  }
  for (int j = 0; j < ntc_y; j++) {
		for (int i = 0; i < ntc_x; i++) {
			S[ntc_x * j + i] = 0.0; (*MFLUX)[ntc_x * j + i] = 0.0;
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
  clock_t start_t, end_t;
  long double ERR = 1.0;

  // ITERATIVE PROCESS
  printf("\n3. ITERATIVE PROCESS:\n\n");
	printf("ITER\tMAX_ERR (MFLUX)\n");
  start_t = clock(); *ITER = -1;
  while (ERR > tol && *ITER < 10000){
    ERR = 0.0; *ITER = *ITER + 1;

    for (int m = 0; m < M; m++){

      // 1. SW - > NE SWEEP
      if (m < M / 4) {
				int j_b = 0, j_f = j_b + 1, i_b = 0, i_f = i_b + 1, nc_y, nc_x, z;
				long double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q, len_y, len_x;
				long double mflow, xflow, yflow;
				miu = MIU[m]; theta = THETA[m];
				for (int ry = 0; ry < nyr; ry++) {
					len_y = YDOM[2 * ry];                  nc_y = (int)YDOM[2 * ry + 1]; 
					h_y = (long double)(len_y / nc_y);
					for (int j = 0; j < nc_y; j++) {
						i_b = 0; i_f = i_b + 1;
						for (int rx = 0; rx < nxr; rx++) {
							len_x = XDOM[2 * rx];          nc_x = (int)XDOM[2 * rx + 1]; 
							h_x = (long double)(len_x / nc_x);    z = ZMAP[nxr * ry + rx];
							q = QMAP[nxr * ry + rx];       st = ZON[2 * z];
							alfa_x = fabs(h_x * st / miu); alfa_y = fabs(h_y * st / theta);
							for (int i = 0; i < nc_x; i++) {
								xflow = (*XFLOW)[M * ((ntc_x + 1) * j_b + i_b) + m];
								yflow = (*YFLOW)[M*(ntc_x*j_b + i_b) + m];
								mflow = (xflow*2/alfa_x + yflow*2/alfa_y + 
									     S[ntc_x*j_b + i_b] + q/st) / (1 +
									     2/alfa_x + 2/alfa_y);
								(*MFLOW)[M * (ntc_x * j_b + i_b) + m] = mflow;
								(*XFLOW)[M * ((ntc_x + 1) * j_b + i_f) + m] = 2 * mflow - xflow;
								(*YFLOW)[M * (ntc_x * j_f + i_b) + m] = 2 * mflow - yflow;
								i_b = i_b + 1; i_f = i_b + 1;
							}
						}
						j_b = j_b + 1; j_f = j_b + 1;
					}
				}
			}

      // 2. SE -> NW SWEEP
      else if (m >= M / 4 && m < M / 2) {
				int j_b = 0, j_f = j_b + 1, i_b = ntc_x - 1, i_f = i_b + 1, nc_y, nc_x, z;
				long double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q, len_y, len_x;
				long double mflow, xflow, yflow;
				miu = MIU[m]; theta = THETA[m];
				for (int ry = 0; ry < nyr; ry++) {
					len_y = YDOM[2 * ry];                  nc_y = (int)YDOM[2 * ry + 1];
					h_y = (long double)(len_y / nc_y);
					for (int j = 0; j < nc_y; j++) {
						i_b = ntc_x - 1; i_f = i_b + 1;
						for (int rx = nxr-1; rx >= 0; rx--) {
							len_x = XDOM[2 * rx];          nc_x = (int)XDOM[2 * rx + 1];
							h_x = (long double)(len_x / nc_x);    z = ZMAP[nxr * ry + rx];
							q = QMAP[nxr * ry + rx];       st = ZON[2 * z];
							alfa_x = fabs(h_x * st / miu); alfa_y = fabs(h_y * st / theta);
							for (int i = 0; i < nc_x; i++) {
								xflow = (*XFLOW)[M * ((ntc_x + 1) * j_b + i_f) + m];
								yflow = (*YFLOW)[M * (ntc_x * j_b + i_b) + m];
								mflow = (xflow * 2 / alfa_x + yflow * 2 / alfa_y +
									     S[ntc_x * j_b + i_b] + q / st) / (1 +
										 2 / alfa_x + 2 / alfa_y);
								(*MFLOW)[M * (ntc_x * j_b + i_b) + m] = mflow;
								(*XFLOW)[M * ((ntc_x + 1) * j_b + i_b) + m] = 2 * mflow - xflow;
								(*YFLOW)[M * (ntc_x * j_f + i_b) + m] = 2 * mflow - yflow;
								i_b = i_b - 1; i_f = i_b + 1;
							}
						}
						j_b = j_b + 1; j_f = j_b + 1;
					}
				}
			}

      // 3. NE -> SW SWEEP
      else if (m >= M / 2 && m < 3 * M / 4) {
				int j_b = ntc_y - 1, j_f = j_b + 1, i_b = ntc_x - 1, i_f = i_b + 1, nc_y, nc_x, z;
				long double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q, len_y, len_x;
				long double mflow, xflow, yflow;
				miu = MIU[m]; theta = THETA[m];
				for (int ry = nyr-1; ry >= 0; ry--) {
					len_y = YDOM[2 * ry];                  nc_y = (int)YDOM[2 * ry + 1];
					h_y = (long double)(len_y / nc_y);
					for (int j = 0; j < nc_y; j++) {
						i_b = ntc_x - 1; i_f = i_b + 1;
						for (int rx = nxr - 1; rx >= 0; rx--) {
							len_x = XDOM[2 * rx];          nc_x = (int)XDOM[2 * rx + 1];
							h_x = (long double)(len_x / nc_x);    z = ZMAP[nxr * ry + rx];
							q = QMAP[nxr * ry + rx];       st = ZON[2 * z];
							alfa_x = fabs(h_x * st / miu); alfa_y = fabs(h_y * st / theta);
							for (int i = 0; i < nc_x; i++) {
								xflow = (*XFLOW)[M * ((ntc_x + 1) * j_b + i_f) + m];
								yflow = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
								mflow = (xflow * 2 / alfa_x + yflow * 2 / alfa_y +
									     S[ntc_x * j_b + i_b] + q / st) / (1 +
										 2 / alfa_x + 2 / alfa_y);
								(*MFLOW)[M * (ntc_x * j_b + i_b) + m] = mflow;
								(*XFLOW)[M * ((ntc_x + 1) * j_b + i_b) + m] = 2 * mflow - xflow;
								(*YFLOW)[M * (ntc_x * j_b + i_b) + m] = 2 * mflow - yflow;
								i_b = i_b - 1; i_f = i_b + 1;
							}
						}
						j_b = j_b - 1; j_f = j_b + 1;
					}
				}
			}

      // 4. NW -> SE SWEEP
      else if (m >= 3 * M / 4 && m < M) {
			  int j_b = ntc_y - 1, j_f = j_b + 1, i_b = 0, i_f = i_b + 1, nc_y, nc_x, z;
				long double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q, len_y, len_x;
				long double mflow, xflow, yflow;
			  miu = MIU[m]; theta = THETA[m];
			  for (int ry = nyr - 1; ry >= 0; ry--) {
				  len_y = YDOM[2 * ry];                  nc_y = (int)YDOM[2 * ry + 1];
				  h_y = (long double)(len_y / nc_y);
				  for (int j = 0; j < nc_y; j++) {
					  i_b = 0; i_f = i_b + 1;
					  for (int rx = 0; rx < nxr; rx++) {
						  len_x = XDOM[2 * rx];          nc_x = (int)XDOM[2 * rx + 1];
						  h_x = (long double)(len_x / nc_x);    z = ZMAP[nxr * ry + rx];
							q = QMAP[nxr * ry + rx];       st = ZON[2 * z];
						  alfa_x = fabs(h_x * st / miu); alfa_y = fabs(h_y * st / theta);
						  for (int i = 0; i < nc_x; i++) {
							  xflow = (*XFLOW)[M * ((ntc_x + 1) * j_b + i_b) + m];
							  yflow = (*YFLOW)[M * (ntc_x * j_f + i_b) + m];
							  mflow = (xflow * 2 / alfa_x + yflow * 2 / alfa_y +
								         S[ntc_x * j_b + i_b] + q / st) / (1 +
									     2 / alfa_x + 2 / alfa_y);
							  (*MFLOW)[M * (ntc_x * j_b + i_b) + m] = mflow;
							  (*XFLOW)[M * ((ntc_x + 1) * j_b + i_f) + m] = 2 * mflow - xflow;
								(*YFLOW)[M * (ntc_x * j_b + i_b) + m] = 2 * mflow - yflow;
							  i_b = i_b + 1; i_f = i_b + 1;
						  }
					  }
					  j_b = j_b - 1; j_f = j_b + 1;
				  }
			  }
			}
    } // endfor m

    // REFLECTIVE BOUNDARY CONDITIONS
		int j_b = 0, nc_y;
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
		int i_b = 0, nc_x;
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
		long double faux, mflux, w, st, ss;
		int z;
		j_b = 0;
		for (int ry = 0; ry < nyr; ry++) {
			nc_y = (int)YDOM[2 * ry + 1];
			for (int j = 0; j < nc_y; j++) {
				i_b = 0;
				for (int rx = 0; rx < nxr; rx++) {
					nc_x = (int)XDOM[2 * rx + 1];       z = ZMAP[nxr * ry + rx];
					st = ZON[2 * z];               ss = ZON[2 * z + 1];
					for (int i = 0; i < nc_x; i++) {
						faux = (*MFLUX)[ntc_x * j_b + i_b]; mflux = 0.0;
						for (int m = 0; m < M; m++) {
							w = W[m];
							mflux = mflux + w * (*MFLOW)[M * (ntc_x * j_b + i_b) + m];
						}
						mflux = 0.25 * mflux;
						(*MFLUX)[ntc_x * j_b + i_b] = mflux;
						
						if (fabs(1 - faux / mflux) > ERR) ERR = fabs(1 - faux / mflux);
						
						S[ntc_x * j_b + i_b] = ss*mflux/st;
						i_b = i_b + 1;
					}
				}
				j_b = j_b + 1;
			}
		}

    printf("%d\t%.5Le\n", *ITER, ERR);
  }
  end_t = clock();
  printf("END %f\n", (double)clock()/10000000);
  *cpu_time = (long double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("CPU TIME = %.5Le\n\n", *cpu_time);

  if (S != NULL) free(S);

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

  // STATUS
  printf("{\n\"STATUS\": %d", status);

  if (status != 0) printf("\n}\n");
  else {

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

///////////////////////////////////////////////////////////////////////////////////////////////