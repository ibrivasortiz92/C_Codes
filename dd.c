#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

/////////////////////////////////////////////////////////// PROTOTYPES

// LOAD PROBLEM FROM FILE
int input_by_txt(int *N, int *nz, double **ZON, int *nxr, int **XDOM, int *nyr, int **YDOM, int **ZMAP, double **QMAP, double *BC, double *tol, const char *filename);

// FREE MEMORY REQUESTED
void free_memory(double **ZON, int **XDOM, int **YDOM, int **ZMAP, double **QMAP, double **MIU, double **THETA, double **CHI, double **W, double **MFLUX, double **MFLOW, double **XFLOW, double **YFLOW);

// GENERATE QUADRATURE
void quad(int N, double **MIU, double **THETA, double **CHI, double **W);

// DD_2D METHOD
void dd_2d(int N, int nz, double ZON[], int nxr, int XDOM[], int nyr, int YDOM[], int ZMAP[], double QMAP[], double BC[], double tol, double MIU[], double THETA[], double W[], double **MFLUX, double **MFLOW, double **XFLOW, double **YFLOW);

// POST PROCESSING
void post_processing(int N, int nz, double ZON[], int nxr, int XDOM[], int nyr, int YDOM[], int ZMAP[], double QMAP[], double BC[], double MIU[], double THETA[], double W[], double MFLUX[], double MFLOW[], double XFLOW[], double YFLOW[]);

int main(int argc, char *argv[]){
  if (argc != 2){
    puts("Usage: ./dd <input.txt>");
    return 1;
  }

  // INPUT VARIABLES
  int N;                // Quadrature order
  int nz;               // Number of zones
  double *ZON = NULL;   // Zone entries
  int nxr;              // Number of regions in X
  int *XDOM = NULL;     // X Region entries
  int nyr;              // Number of regions in Y
  int *YDOM = NULL;     // Y Region entries
  int *ZMAP = NULL;     // Zone map
  double *QMAP = NULL;  // External source map
  double BC[4];         // Boundary conditions
  double tol;           // Tolerance
  const char *filename; // File name to load

  filename = argv[1];

  // LOAD PROBLEM
  if (input_by_txt(&N, &nz, &ZON, &nxr, &XDOM, &nyr, &YDOM, &ZMAP, &QMAP, BC, &tol, filename) > 0) return 1;

  // QUADRATURE VARIABLES
  double *MIU = NULL;   // Ordinates in X
  double *THETA = NULL; // Ordinates in Y
  double *CHI = NULL;   // Ordinates in Z
  double *W = NULL;     // Weight

  // GENERATE AND PRINT QUADRATURE
  quad(N, &MIU, &THETA, &CHI, &W);

  // PRINCIPAL VARIABLES
  double *MFLUX = NULL; // Escalar flux in the nodes
  double *MFLOW = NULL; // Angular flux in the nodes
  double *XFLOW = NULL; // Angular flux at the y edges
  double *YFLOW = NULL; // Angular flux at the x edges

  // METHOD
  dd_2d(N, nz, ZON, nxr, XDOM, nyr, YDOM, ZMAP, QMAP, BC, tol, MIU, THETA, W, &MFLUX, &MFLOW, &XFLOW, &YFLOW);

  // JSON OUTPUT
  post_processing(N, nz, ZON, nxr, XDOM, nyr, YDOM, ZMAP, QMAP, BC, MIU, THETA, W, MFLUX, MFLOW, XFLOW, YFLOW);

  // FREE MEMORY
  free_memory(&ZON, &XDOM, &YDOM, &ZMAP, &QMAP,
              &MIU, &THETA, &CHI, &W,
              &MFLUX, &MFLOW, &XFLOW, &YFLOW);

  return 0;
}

////////////////////////////////////////////////////////// LOAD PROBLEM FROM FILE
int input_by_txt(int *N,             // Quadrature order
                int *nz,             // Number of zones
                double **ZON,        // Zone entries
                int *nxr,            // Number of regions in X
                int **XDOM,          // X Region entries
                int *nyr,            // Number of regions in Y
                int **YDOM,          // Y Region entries
                int **ZMAP,          // Zone map
                double **QMAP,       // External source map
                double *BC,          // Boundary conditions
                double *tol,         // Tolerance
                const char *filename // File name to load
                ){

  FILE *cfPtr = NULL;
  if ((cfPtr = fopen(filename, "r")) == NULL){
      puts ("Error opening input file!");
      return 1;
  }

  // QUADRATURE ORDER
  if (fscanf(cfPtr, "%d", N) == 0){
      puts("Error reading input file!");
      return 2;
  }
	if (*N < 2 || *N > 18){
		puts("Invalid quadrature order (N)!\n\n"
		     "\tValue must be an even number between 2 and 18.");
		return 2;
	}

  // NUMBER OF ZONES
  if (fscanf(cfPtr, "%d", nz) == 0){
      puts("Error reading input file!");
      return 2;
  }
	if (*nz < 1){
		puts("Invalid number of zones (NZ)!\n\n"
		     "\tValue must be greater than 1.");
		return 2;
	}

  // ZONE FILLING
  *ZON = malloc(sizeof(double) * (*nz * 2)); assert(*ZON != NULL);
  double st, ss;
	int zone_count = 0;
  for (int z = 0; z < *nz; z++){
      if (fscanf(cfPtr, "%lf %lf", &st, &ss) == 0){
          puts("Error reading input file!");
    free(*ZON);
          return 2;
      }
      (*ZON)[z * 2] = st; (*ZON)[z * 2 + 1] = ss;
  zone_count = zone_count + 1;
  }
	if (zone_count != *nz){
		free(*ZON);
		puts("Error reading the entries of ZON:\n\n"
		     "\tMust be provided the same number of rows as the number of zones.");
		return 2;
	}

  // NUMBER OF REGIONS IN X
  if (fscanf(cfPtr, "%d", nxr) == 0){
      puts("Error reading input file!");
  free(*ZON);
      return 2;
  }
	if (*nxr < 1){
		free(*ZON);
		puts("Invalid number of regions in X (NR_X)!\n\n"
		     "\tValue must be greater than 1.");
		return 2;
	}

  // REGION FILLING IN X
  *XDOM = malloc(sizeof(int) * (*nxr * 2)); assert(*XDOM != NULL);
  int len, nodes, rx_count = 0;;
  for (int xr = 0; xr < *nxr; xr++){
      if (fscanf(cfPtr, "%d %d", &len, &nodes) == 0){
          puts("Error reading input file!");
    free(*ZON); free(*XDOM);
          return 2;
      }
      (*XDOM)[xr * 2] = len; (*XDOM)[xr * 2 + 1] = nodes;
  rx_count = rx_count + 1;
  }
	if (rx_count != *nxr){
		free(*ZON); free(*XDOM);
		puts("Error reading the entries of XDOM:\n\n"
		     "\tMust be provided the same number of rows as the number of regions in X.");
		return 2;
	}

  // NUMBER OF REGIONS IN Y
  if (fscanf(cfPtr, "%d", nyr) == 0){
      puts("Error reading input file!");
  free(*ZON); free(*XDOM);
      return 2;
  }
	if (*nxr < 1){
		free(*ZON); free(*XDOM);
		puts("Invalid number of regions in Y (NR_Y)!\n\n"
		     "\tValue must be greater than 1.");
		return 2;
	}

  // REGION FILLING IN Y
  *YDOM = malloc(sizeof(int) * (*nyr * 2)); assert(*YDOM != NULL);
	int ry_count = 0;
  for (int yr = 0; yr < *nyr; yr++){
      if (fscanf(cfPtr, "%d %d", &len, &nodes) == 0){
          puts("Error reading input file!");
    free(*ZON); free(*XDOM); free(*YDOM);
          return 2;
      }
      (*YDOM)[yr * 2] = len; (*YDOM)[yr * 2 + 1] = nodes;
  ry_count = ry_count + 1;
  }
	if (ry_count != *nyr){
		free(*ZON); free(*XDOM); free(*YDOM);
		puts("Error reading the entries of YDOM:\n\n"
		     "\tMust be provided the same number of rows as the number of regions in Y.");
		return 2;
	}

    // ZONE MAPPING
	*ZMAP = malloc(sizeof(int) * (*nyr) * (*nxr)); assert(*ZMAP != NULL);
	int z, entry_z_count = 0;
  for (int yr = 0; yr < *nyr; yr++) {
		for (int xr = 0; xr < *nxr; xr++) {
			if (fscanf(cfPtr, "%d", &z) == 0) {
				puts("Error reading input file!");
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
		puts("Error reading the entries of ZMAP:\n\n"
		     "\tThe entries in ZMAP don't match with the regions of the problem.");
		return 2;
	}

  // EXTERNAL SOURCE MAPPING
	*QMAP = malloc(sizeof(double) * (*nyr) * (*nxr)); assert(*QMAP != NULL);
	double q;
	int entry_q_count = 0;
  for (int yr = 0; yr < *nyr; yr++) {
		for (int xr = 0; xr < *nxr; xr++) {
			if (fscanf(cfPtr, "%lf", &q) == 0) {
				puts("Error reading input file!");
				free(*ZON); free(*XDOM); free(*YDOM);
				free(*ZMAP); free(*QMAP);
                return 2;
			}
			(*QMAP)[yr * (*nxr) + xr] = q;
			entry_q_count = entry_q_count + 1;
		}
  }
	if (entry_q_count != (*nyr) * (*nxr)){
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
		puts("Error reading the entries of QMAP:\n\n"
		     "\tThe entries in QMAP don't match with the regions of the problem.");
		return 2;
	}

  // BOUNDARY CONDITIONS
  double cond;
	int bc_count = 0;
	for (int c = 0; c < 4; c++) {
		if (fscanf(cfPtr, "%lf", &cond) == 0) {
			puts("Error reading input file!");
			free(*ZON); free(*XDOM); free(*YDOM);
			free(*ZMAP); free(*QMAP);
            return 2;
		}
		BC[c] = cond;
		bc_count = bc_count + 1;
		if (cond < 0.0 ){
			if (cond != -1.0){
				free(*ZON); free(*XDOM); free(*YDOM);
			    free(*ZMAP); free(*QMAP);
				puts("Invalid boudary condition!\n\n"
		     		"\tValue must be -1.0, 0.0, or a positive number.");
				return 2;
			}
		}
	}
	if (bc_count != 4){
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
		puts("The number of boundary conditions must be 4.\n\n");
		return 2;
	}

  // TOLERANCE
	if (fscanf(cfPtr, "%lf", tol) == 0) {
		puts("Error reading input file!");
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
        return 2;
	}
	if (*tol <= 0.0 || *tol >= 1.0){
		free(*ZON); free(*XDOM); free(*YDOM);
		free(*ZMAP); free(*QMAP);
		puts("The tolerance (TOL) must be a small number between 0.0 and 1.0.\n\n");
		return 2;
	}

	fclose(cfPtr);

  return 0;
}
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////// FREE MEMORY REQUESTED
void free_memory(double **ZON, int **XDOM, int **YDOM, int **ZMAP, double **QMAP,
                 double **MIU, double **THETA, double **CHI, double **W,
				         double **MFLUX, double **MFLOW, double **XFLOW, double **YFLOW){

  // INPUT VARIABLES
  free(*ZON);     free(*XDOM);    free(*YDOM);    free(*ZMAP);    free(*QMAP);

  // QUADRATURE VARIABLES
  free(*MIU);     free(*THETA);   free(*CHI);     free(*W);

	// PRINCIPAL VARIABLES
	free(*MFLUX);	free(*MFLOW);	free(*XFLOW);	free(*YFLOW);
}

///////////////////////////////////////////////////////////////// GENERATE QUADRATURE
void quad(int N,           // Quadrature order
          double **MIU,    // Ordinates in X
          double **THETA,  // Ordinates in Y
          double **CHI,    // Ordinates in Z
          double **W       // Weight
          ){
    
  // DIRECTIONS IN THE XY PLANE
  int M = N * (N + 2) / 2;
  
  // ALLOCATE MEMORY FOR THE ORDINATES
  *MIU = malloc(sizeof(double) * M); assert(*MIU != NULL);
  *THETA = malloc(sizeof(double) * M); assert(*THETA != NULL);
  *CHI = malloc(sizeof(double) * M); assert(*CHI != NULL);
  *W = malloc(sizeof(double) * M); assert(*W != NULL);

  // PREFIXED VALUES
  double chi[10] = { 0.0 }, wlist[10] = { 0.0 };
	if (N == 2) {
		chi[0] = 0.5773502692; wlist[0] = 1.0;
	}
	else if (N == 4) {
		chi[0] = 0.3500212; wlist[0] = (double)1/3;
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
            (*CHI)[d] = chi[n]; (*W)[d] = 0.0;
			d = d + 1; aux0 = aux0 + 1;
		}
	}

  // WEIGHT FILLING
	int p = 0; double aux1, aux2, aux3;
	for (int n = 0; n < M / 4; n++) {
		if ((*W)[n] == 0.0) {
			(*W)[n] = wlist[p];
			aux1 = (*MIU)[n]; aux2 = (*THETA)[n]; aux3 = (*CHI)[n];
			for (int m = 0; m < M / 4; m++) {
				if (aux1 == (*THETA)[m] && aux2 == (*MIU)[m] && aux3 == (*CHI)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*THETA)[m] && aux2 == (*CHI)[m] && aux3 == (*MIU)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*CHI)[m] && aux2 == (*THETA)[m] && aux3 == (*MIU)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*CHI)[m] && aux2 == (*MIU)[m] && aux3 == (*THETA)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*MIU)[m] && aux2 == (*CHI)[m] && aux3 == (*THETA)[m]) {
					(*W)[m] = wlist[p];
				}
				if (aux1 == (*MIU)[m] && aux2 == (*THETA)[m] && aux3 == (*CHI)[m]) {
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
				(*CHI)[aux4] = (*CHI)[m]; (*W)[aux4] = (*W)[m];
			}
			else if (q == 3) {
				(*MIU)[aux4] = -(*MIU)[m]; (*THETA)[aux4] = -(*THETA)[m];
				(*CHI)[aux4] = (*CHI)[m]; (*W)[aux4] = (*W)[m];
			}
			else if (q == 4) {
				(*MIU)[aux4] = (*MIU)[m]; (*THETA)[aux4] = -(*THETA)[m];
				(*CHI)[aux4] = (*CHI)[m]; (*W)[aux4] = (*W)[m];
			}
			aux4 = aux4 + 1;
		}
	}
    
}

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////   METODO DD_2D

/* Implementation of the Diamond Difference method for monoenergetic neutron transport problems   */ 
/* with fixed source, isotropic scattering and in xy geometry.                                    */

void dd_2d(int N,          // Quadrature order
           int nz,         // Number of zones
           double ZON[],   // Zone entries
           int nxr,        // Number of regions in X
           int XDOM[],     // X Region entries
           int nyr,        // Number of regions in Y
           int YDOM[],     // Y Region entries
           int ZMAP[],     // Zone map
           double QMAP[],  // External source map
           double BC[],    // Boundary conditions
           double tol,     // Tolerance
           double MIU[],   // Ordinates in X
           double THETA[], // Ordinates in Y
           double W[],     // Weight
           double **MFLUX, // Escalar flux in the nodes
           double **MFLOW, // Angular flux in the nodes
           double **XFLOW, // Angular flux at the y edges
           double **YFLOW  // Angular flux at the x edges
           ){
    
  // INITIALIZATION
  int ntc_x = 0, ntc_y = 0;
  for (int rx = 0; rx < nxr; rx++) {
		ntc_x = ntc_x + XDOM[2*rx + 1];
	}
	for (int ry = 0; ry < nyr; ry++) {
		ntc_y = ntc_y + YDOM[2 * ry + 1];
	}
  int M = N * (N + 2) / 2;
  double *S = NULL;      // Source term
  S = malloc(sizeof(double) * (ntc_y * ntc_x)); assert(S != NULL);
  *MFLUX = malloc(sizeof(double) * (ntc_y * ntc_x)); assert(*MFLUX != NULL);
	*MFLOW = malloc(sizeof(double) * (ntc_y * ntc_x) * M); assert(*MFLOW != NULL);
	*XFLOW = malloc(sizeof(double) * (ntc_y * (ntc_x + 1)) * M); assert(*XFLOW != NULL);
	*YFLOW = malloc(sizeof(double) * ((ntc_y + 1) * ntc_x) * M); assert(*YFLOW != NULL);
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

  // ITERATIVE PROCESS
  double ERR = 1.0;
	int ITER = -1;
  while (ERR > tol) {
    ERR = 0.0; ITER = ITER + 1;

    for (int m = 0; m < M; m++){

      // 1. SE - > NW SWEEP
      if (m < M / 4) {
				int j_b = 0, j_f = j_b + 1, i_b = 0, i_f = i_b + 1;
				int nc_y, len_y, nc_x, len_x, z;
				double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q;
				double mflow, xflow, yflow;
				miu = MIU[m]; theta = THETA[m];
				for (int ry = 0; ry < nyr; ry++) {
					len_y = YDOM[2 * ry];                  nc_y = YDOM[2 * ry + 1]; 
					h_y = (double)len_y / nc_y;
					for (int j = 0; j < nc_y; j++) {
						i_b = 0; i_f = i_b + 1;
						for (int rx = 0; rx < nxr; rx++) {
							len_x = XDOM[2 * rx];          nc_x = XDOM[2 * rx + 1]; 
							h_x = (double)len_x / nc_x;    z = ZMAP[nxr * ry + rx];
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
				int j_b = 0, j_f = j_b + 1, i_b = ntc_x - 1, i_f = i_b + 1;
				int nc_y, len_y, nc_x, len_x, z;
				double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q;
				double mflow, xflow, yflow;
				miu = MIU[m]; theta = THETA[m];
				for (int ry = 0; ry < nyr; ry++) {
					len_y = YDOM[2 * ry];                  nc_y = YDOM[2 * ry + 1];
					h_y = (double)len_y / nc_y;
					for (int j = 0; j < nc_y; j++) {
						i_b = ntc_x - 1; i_f = i_b + 1;
						for (int rx = nxr-1; rx >= 0; rx--) {
							len_x = XDOM[2 * rx];          nc_x = XDOM[2 * rx + 1];
							h_x = (double)len_x / nc_x;    z = ZMAP[nxr * ry + rx];
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
				int j_b = ntc_y - 1, j_f = j_b + 1, i_b = ntc_x - 1, i_f = i_b + 1;
				int nc_y, len_y, nc_x, len_x, z;
				double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q;
				double mflow, xflow, yflow;
				miu = MIU[m]; theta = THETA[m];
				for (int ry = nyr-1; ry >= 0; ry--) {
					len_y = YDOM[2 * ry];                  nc_y = YDOM[2 * ry + 1];
					h_y = (double)len_y / nc_y;
					for (int j = 0; j < nc_y; j++) {
						i_b = ntc_x - 1; i_f = i_b + 1;
						for (int rx = nxr - 1; rx >= 0; rx--) {
							len_x = XDOM[2 * rx];          nc_x = XDOM[2 * rx + 1];
							h_x = (double)len_x / nc_x;    z = ZMAP[nxr * ry + rx];
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
			  int j_b = ntc_y - 1, j_f = j_b + 1, i_b = 0, i_f = i_b + 1;
			  int nc_y, len_y, nc_x, len_x, z;
				double h_y, h_x, miu, theta, alfa_x, alfa_y, st, q;
				double mflow, xflow, yflow;
			  miu = MIU[m]; theta = THETA[m];
			  for (int ry = nyr - 1; ry >= 0; ry--) {
				  len_y = YDOM[2 * ry];                  nc_y = YDOM[2 * ry + 1];
				  h_y = (double)len_y / nc_y;
				  for (int j = 0; j < nc_y; j++) {
					  i_b = 0; i_f = i_b + 1;
					  for (int rx = 0; rx < nxr; rx++) {
						  len_x = XDOM[2 * rx];          nc_x = XDOM[2 * rx + 1];
						  h_x = (double)len_x / nc_x;    z = ZMAP[nxr * ry + rx];
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
			nc_y = YDOM[2 * ry + 1];
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
			nc_x = XDOM[2 * rx + 1];
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
		double faux, mflux, w, st, ss;
		int z;
		j_b = 0;
		for (int ry = 0; ry < nyr; ry++) {
			nc_y = YDOM[2 * ry + 1];
			for (int j = 0; j < nc_y; j++) {
				i_b = 0;
				for (int rx = 0; rx < nxr; rx++) {
					nc_x = XDOM[2 * rx + 1];       z = ZMAP[nxr * ry + rx];
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
  }

	// FREE MEMORY
	free(S);
}
////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// POST PROCESSING

void post_processing(int N,          // Quadrature order
                     int nz,         // Number of zones
                     double ZON[],   // Zone entries
                     int nxr,        // Number of regions in X
                     int XDOM[],     // X Region entries
                     int nyr,        // Number of regions in Y
                     int YDOM[],     // Y Region entries
                     int ZMAP[],     // Zone map
                     double QMAP[],  // External source map
                     double BC[],    // Boundary conditions
                     double MIU[],   // Ordinates in X
                     double THETA[], // Ordinates in Y
                     double W[],     // Weight
	                 double MFLUX[], // Escalar flux in the nodes
                     double MFLOW[], // Angular flux in the nodes
                     double XFLOW[], // Angular flux at the y edges
                     double YFLOW[]  // Angular flux at the x edges
					 ){
	
	// INICIALIZATION
	int j_b = 0, i_b = 0, j_f, i_f, ntc_x = 0, ntc_y = 0;;
	int len_y, nc_y, len_x, nc_x, z;
	double h_y, h_x, area_r, sa, miu, theta, w;
  for (int rx = 0; rx < nxr; rx++) {
		ntc_x = ntc_x + XDOM[2*rx + 1];
	}
	for (int ry = 0; ry < nyr; ry++) {
		ntc_y = ntc_y + YDOM[2 * ry + 1];
	}
  int M = N * (N + 2) / 2;

	// MFLUX
  printf("{\n\"MFLUX\": [\n");
	for(int ry = 0; ry < nyr; ry++){
		len_y = YDOM[2*ry]; nc_y = YDOM[2*ry + 1]; h_y = (double)len_y/nc_y;
		for(int j = 0; j < nc_y; j++){
			i_b = 0;
      printf("[");
			for(int rx = 0; rx < nxr; rx++){
				len_x = XDOM[2*rx]; nc_x = XDOM[2*rx + 1]; h_x = (double)len_x/nc_x; 
				area_r = (double)len_y*len_x; z = ZMAP[nxr*ry + rx];
				sa = ZON[2*z] - ZON[2*z + 1];
				for(int i = 0; i < nc_x; i++){
          if (i_b == ntc_x - 1) {
						if (j_b == ntc_y - 1) printf(" %.10e ]\n", MFLUX[ntc_x*j_b + i_b]);
						else printf(" %.10e ],\n", MFLUX[ntc_x*j_b + i_b]);
					}
          else printf(" %.10e,", MFLUX[ntc_x*j_b + i_b]);
          i_b = i_b + 1;
				}
			}
			j_b = j_b + 1;
		}
	}
  printf("],\n");

  // MFLOW
  printf("\"MFLOW\": [\n");
  for (int m = 0; m < M; m++){
    j_b = 0;
    printf("[\n");
    for(int ry = 0; ry < nyr; ry++){
		  len_y = YDOM[2*ry]; nc_y = YDOM[2*ry + 1];
		  for(int j = 0; j < nc_y; j++){
			  i_b = 0;
        printf("[");
			  for(int rx = 0; rx < nxr; rx++){
				  len_x = XDOM[2*rx]; nc_x = XDOM[2*rx + 1];
				  for(int i = 0; i < nc_x; i++){
            if (i_b == ntc_x - 1) {
							if (j_b == ntc_y - 1) printf(" %.10e ]\n", MFLOW[M * (ntc_x * j_b + i_b) + m]);
							else printf(" %.10e ],\n", MFLOW[M * (ntc_x * j_b + i_b) + m]);
						}
            else printf(" %.10e,", MFLOW[M * (ntc_x * j_b + i_b) + m]);
            i_b = i_b + 1;
				  }
			  }
			  j_b = j_b + 1;
		  }
	  }
		if (m == M-1) printf("]\n");
    else printf("],\n");
  }
  printf("],\n");

  // XFLOW
  printf("\"XFLOW\": [\n");
  for (int m = 0; m < M; m++){
    j_b = 0;
    printf("[\n");
    for(int ry = 0; ry < nyr; ry++){
		  len_y = YDOM[2*ry]; nc_y = YDOM[2*ry + 1];
		  for(int j = 0; j < nc_y; j++){
			  i_b = 0;
        printf("[");
			  for(int rx = 0; rx < nxr; rx++){
				  len_x = XDOM[2*rx]; nc_x = XDOM[2*rx + 1];
				  for(int i = 0; i < nc_x; i++){
            printf(" %.10e,", XFLOW[M * ((ntc_x + 1) * j_b + i_b) + m]);
            i_b = i_b + 1;
				  }
			  }
        if (i_b == ntc_x) {
					if (j_b == ntc_y - 1) printf(" %.10e ]\n", XFLOW[M * ((ntc_x + 1) * j_b + i_b) + m]);
					else printf(" %.10e ],\n", XFLOW[M * ((ntc_x + 1) * j_b + i_b) + m]);
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
		  len_y = YDOM[2*ry]; nc_y = YDOM[2*ry + 1];
		  for(int j = 0; j < nc_y; j++){
			  i_b = 0;
        printf("[");
			  for(int rx = 0; rx < nxr; rx++){
				  len_x = XDOM[2*rx]; nc_x = XDOM[2*rx + 1];
				  for(int i = 0; i < nc_x; i++){
            if (i_b == ntc_x - 1) printf(" %.10e ],\n", YFLOW[M * (ntc_x * j_b + i_b) + m]);
            else printf(" %.10e,", YFLOW[M * (ntc_x * j_b + i_b) + m]);
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
        len_x = XDOM[2*rx]; nc_x = XDOM[2*rx + 1];
        for(int i = 0; i < nc_x; i++){
          if (i_b == ntc_x - 1) printf(" %.10e ]\n", YFLOW[M * (ntc_x * j_b + i_b) + m]);
          else printf(" %.10e,", YFLOW[M * (ntc_x * j_b + i_b) + m]);
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
/////////////////////////////////////////////////////////////////////////////////////
