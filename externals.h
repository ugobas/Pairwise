extern int Naa, Na2;
extern int PDB;     // Two modes: PDB (statistics of the PDB) or ALI (MSA)
extern int REM;     // Use 1st (1), 2nd (2), 3rd (3) moment of misfold energy 
extern float TEMP;
extern int *ij_bin; // Range of values of |i-j| for PDB option
extern int IJ_MIN; // Minimal value of |i-j|
extern int **label_ij; // For alignments, pairs corresponding to a.a. ij i<j
extern float **E_cont_gap, **E_cont_T; // Energy matrix
extern float *E_cont_1;
extern int INDIRECT;   // Compute indirect correlations through third residue?
extern int Q_GLOB;  // 
extern float NTHR;  // Threshold for excluding small pairs
//extern float S_DEPEND; // Estimated number of dependent sequences

// Concact statistics
// <Nc>, <Nc^2>, <Nc^3>
extern float *Nc1L, *Nc2L, *Nc3L;
extern int Verbose;
extern char cont_type;
extern float cont_thr, cont_thr2;
