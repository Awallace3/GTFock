#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
//#include <ga.h>
#include <mkl.h>
#include <omp.h>
#include <mpi.h>
#include <mkl_scalapack.h>
#include <mkl_blacs.h>

#include "config.h"
#include "one_electron.h"

#include "GTMatrix.h"

#define min(a, b) ((a) < (b) ? (a): (b))

inline void matrix_block_write(double *matrix, int startrow,
                               int startcol, int ldm,
                               double *block, int nrows, int ncols)
{
    for (int k = 0; k < nrows; k++) {
        for (int l = 0; l < ncols; l++) {
            int i = startrow + k;
            int j = startcol + l;
            matrix[i * ldm + j] = block[l + ncols * k];//Simint
          //matrix[i * ldm + j] = block[k + nrows * l];//OptERD
        }
    }
}


void compute_S(PFock_t pfock, BasisSet_t basis,
               int startshellrow, int endshellrow,
               int startshellcol, int endshellcol,
               int ldS, double *S)
{
//  int nthreads = omp_get_max_threads();
//  OED_t *oed = (OED_t *)malloc(sizeof(OED_t) * nthreads);
//  assert(oed != NULL);
//  for (int i = 0; i < nthreads; i++) {
//      CInt_createOED(basis, &(oed[i]));
//  }
    int start_row_id = pfock->f_startind[startshellrow];
    int start_col_id = pfock->f_startind[startshellcol];

    #pragma omp parallel
    {
        int tid = omp_get_thread_num ();
        #pragma omp for
        for (int A = startshellrow; A <= endshellrow; A++) {
            int row_id_1 = pfock->f_startind[A];
            int row_id_2 = pfock->f_startind[A + 1] - 1;
            int startrow = row_id_1 - start_row_id;
            int nrows = row_id_2 - row_id_1 + 1;
            for (int B = startshellcol; B <= endshellcol; B++) {
                int col_id_1 = pfock->f_startind[B];
                int col_id_2 = pfock->f_startind[B + 1] - 1;
                int startcol = col_id_1 - start_col_id;
                int ncols = col_id_2 - col_id_1 + 1;
                int nints;
                double *integrals;
                CInt_computePairOvl_SIMINT(basis, pfock->simint, tid, A, B, &integrals, &nints);
                if (nints != 0) {
                    matrix_block_write(S, startrow, startcol, ldS,
                                       integrals, nrows, ncols);
                }
            }
        }
    }

//  for (int i = 0; i < nthreads; i++) {
//      CInt_destroyOED(oed[i]);
//  }
//  free(oed);
}


void compute_H(PFock_t pfock, BasisSet_t basis,
               int startshellrow, int endshellrow,
               int startshellcol, int endshellcol,
               int ldH, double *H)
{
//  int nthreads = omp_get_max_threads();
//  OED_t *oed = (OED_t *)malloc(sizeof(OED_t) * nthreads);
//  assert(oed != NULL);
//  for (int i = 0; i < nthreads; i++) {
//      CInt_createOED(basis, &(oed[i]));
//  }
    
    int start_row_id = pfock->f_startind[startshellrow];
    int start_col_id = pfock->f_startind[startshellcol];
    #pragma omp parallel
    {
        int tid = omp_get_thread_num ();
        #pragma omp for
        for (int A = startshellrow; A <= endshellrow; A++) {
            int row_id_1 = pfock->f_startind[A];
            int row_id_2 = pfock->f_startind[A + 1] - 1;
            int startrow = row_id_1 - start_row_id;
            int nrows = row_id_2 - row_id_1 + 1;
            for (int B = startshellcol; B <= endshellcol; B++)
            {
                int col_id_1 = pfock->f_startind[B];
                int col_id_2 = pfock->f_startind[B + 1] - 1;
                int startcol = col_id_1 - start_col_id;
                int ncols = col_id_2 - col_id_1 + 1;
                int nints;
                double *integrals;
                CInt_computePairCoreH_SIMINT(basis, pfock->simint, tid,
                                      A, B, &integrals, &nints);
                if (nints != 0) {
                    matrix_block_write(H, startrow, startcol, ldH,
                                       integrals, nrows, ncols);
                }
            }
        }
    }

//  for (int i = 0; i < nthreads; i++) {
//      CInt_destroyOED(oed[i]);
//  }
//  free(oed);
}


// void my_peig(GTMatrix_t gtm_A, GTMatrix_t gtm_B, int n, int nprow, int npcol, double *eval)
// {
//     int myrank;
//     int ictxt;
//     int myrow;
//     int mycol;
//     int nn;
//     int mm;
//     int descA[9];
//     int descZ[9];
//     int info;
//     int lo[2];
//     int hi[2];
//     int ld;
//     int izero = 0;
//     int ione = 1;
//
//     // init blacs
//     printf(" Starting nprow = %d, npcol = %d\n", nprow, npcol);
//     int nb = MIN(n / nprow, n / npcol);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//     printf("  myrank = %d\n", myrank);
//     Cblacs_pinfo(&nn, &mm);
//     printf("  n = %d, nn = %d, mm = %d\n", n, nn, mm);
//     // Cblacs_get(-1, 0, &ictxt);
//     Cblacs_get(-1, 0, &ictxt);
//     Cblacs_gridinit(&ictxt, "Row-major", nprow, npcol);
//     Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
//     // Cblacs_pcoord(ictxt, myrank, &myrow, &mycol);
//
//     int world_size;
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//
//
//     // nprow = 1;
//     printf(" n = %d, nb = %d, myrow = %d, nprow = %d, npcol = %d\n", n, nb, myrow, nprow, npcol);
//
//     // init matrices
//     // issue is that nprow is 0, should be the NPROCS in the grid
//     int nrows = numroc_(&n, &nb, &myrow, &izero, &nprow);
//     printf("  nrows = %d\n", nrows);
//     int ncols = numroc_(&n, &nb, &mycol, &izero, &npcol);
//     printf("  ncols = %d\n", ncols);
//     int itemp = nrows > 1 ? nrows : 1;
//     // print values that go into descinit
//     printf(" descA: n = %d, nrows = %d, ncols = %d, nb = %d, nprow = %d, npcol = %d, izero = %d, ictxt = %d, itemp = %d, info = %d\n", n, nrows, ncols, nb, nprow, npcol, izero, ictxt, itemp, info);
//     // descinit_(descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
//     // descinit_(descZ, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
//     descinit_(descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
//     descinit_(descZ, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
//     int blocksize = nrows * ncols;
//     printf("  blocksize = %d\n", blocksize);
//     double *A = (double *)aligned_malloc(blocksize * sizeof (double), 64);
//     double *Z = (double *)aligned_malloc(blocksize * sizeof (double), 64);
//     assert(Z != NULL && A != NULL);
//
//     printf("init matrices\n");
//
//     // distribute source matrix
//     GTM_startBatchGet(gtm_A);
//     for (int i = 1; i <= nrows; i += nb) 
//     {
//         lo[0] = indxl2g_(&i, &nb, &myrow, &izero, &nprow) - 1;
//         hi[0] = lo[0] + nb - 1;
//         hi[0] = hi[0] >= n ? n - 1 : hi[0];
//         for (int j = 1; j <= ncols; j += nb) 
//         {
//             lo[1] = indxl2g_(&j, &nb, &mycol, &izero, &npcol) - 1;
//             hi[1] = lo[1] + nb - 1;
//             hi[1] = hi[1] >= n ? n - 1 : hi[1];
//             ld = ncols;
//             GTM_addGetBlockRequest(
//                 gtm_A, 
//                 lo[0], hi[0] - lo[0] + 1,
//                 lo[1], hi[1] - lo[1] + 1,
//                 &(Z[(i - 1) * ncols + j - 1]), ld
//             );
//         }
//     }
//     printf("exec batch get\n");
//     GTM_execBatchGet(gtm_A);
//     GTM_stopBatchGet(gtm_A);
//     GTM_sync(gtm_A);
//     printf("GTM_synced\n");
//     
//     for (int i = 0; i < nrows; i++) 
//     {
//         #pragma omp simd
//         for (int j = 0; j < ncols; j++) 
//             A[j * nrows + i] = Z[i * ncols + j];
//     }
//
//     double t1 = MPI_Wtime();
//     // acquire working space
//     double *work = (double *)aligned_malloc(2 * sizeof (double), 64);
//     assert (work != NULL);
//     int lwork = -1;
// #if 0
//     pdsyev ("V", "U", &n, A, &ione, &ione, descA,
//             eval, Z, &ione, &ione, descZ, work, &lwork, &info);
// #else
//     int liwork = -1;
//     int *iwork = (int *)aligned_malloc(2 * sizeof (int), 64);
//     assert(iwork != NULL);
//     // This call is used to determine the optimal work size
//     pdsyevd_("V", "U", &n, A, &ione, &ione, descA,
//             eval, Z, &ione, &ione, descZ,
//             work, &lwork, iwork, &liwork, &info);    
// #endif
//     printf(" eval[0] = %lf\n", eval[0]);
//
//     // compute eigenvalues and eigenvectors
//     printf("Computing eigenvalues and eigenvectors\n");
//     printf("  lwork = %d, work[0] = %d\n", lwork, (int)work[0]);
//     // lwork = (int)work[0] * 2;
//     lwork = (int)work[0];
//     printf("  lwork = %d, work[0] = %d\n", lwork, (int)work[0]);
//     aligned_free(work);
//     assert(work == NULL);
//
//     work = (double *)aligned_malloc(lwork * sizeof (double), 64);
//     // lwork = malloc_usable_size(work) / sizeof(double);
//     // need to update lwork in case of alignment size increase
//     
//     assert(work != NULL);
//     printf("  lwork = %d\n", lwork);
// #if 0
//     pdsyev ("V", "U", &n, A, &ione, &ione, descA,
//             eval, Z, &ione, &ione, descZ, work, &lwork, &info);
// #else
//     liwork = (int)iwork[0];
//     aligned_free(iwork);
//     iwork = (int *)aligned_malloc(liwork * sizeof (int), 64);
//     // liwork = malloc_usable_size(iwork) / sizeof(int);
//     assert(iwork != NULL);
//     int LOCp_q = ione + n - 1;
//     printf(" lwork = %d, liwork = %d\n", lwork, liwork * sizeof(int));
//     printf(" LOCp_q(upper_index) = %d, max_elements = %d\n", LOCp_q, LOCp_q * LOCp_q);
//     printf("pdsyevd_ call\n");
//     pdsyevd_("V", "U", &n, A, &ione, &ione, descA,
//             eval, Z, &ione, &ione, descZ,
//             work, &lwork, iwork, &liwork, &info); 
//     if (info != 0) {
//         printf("pdsyevd_ failed with info = %d\n", info);
//         return;
//     }
// #endif
//     printf("Computed eigenvalues and eigenvectors\n");
//     double t2 = MPI_Wtime();
//     if (myrank == 0) printf("  pdsyev_ takes %.3lf secs\n", t2 - t1);
//
//     // store desination matrix
//     for (int i = 0; i < nrows; i++) 
//     {
//         for (int j = 0; j < ncols; j++)
//             A[i * ncols + j] = Z[j * nrows + i];
//     }
//     
//     GTM_startBatchPut(gtm_B);
//     for (int i = 1; i <= nrows; i += nb) 
//     {
//         lo[0] = indxl2g_ (&i, &nb, &myrow, &izero, &nprow) - 1;
//         hi[0] = lo[0] + nb - 1;
//         hi[0] = hi[0] >= n ? n - 1 : hi[0];
//         for (int j = 1; j <= ncols; j += nb) 
//         {
//             lo[1] = indxl2g_ (&j, &nb, &mycol, &izero, &npcol) - 1;
//             hi[1] = lo[1] + nb - 1;
//             hi[1] = hi[1] >= n ? n - 1 : hi[1];
//             ld = ncols;
//             GTM_addPutBlockRequest(
//                 gtm_B, 
//                 lo[0], hi[0] - lo[0] + 1,
//                 lo[1], hi[1] - lo[1] + 1,
//                 &(A[(i - 1) * ncols + j - 1]), ld
//             );
//         }
//     }
//     GTM_execBatchPut(gtm_B);
//     GTM_stopBatchPut(gtm_B);
//     GTM_sync(gtm_B);
//     printf("aligned_free\n");
//
//     aligned_free(A);
//     aligned_free(Z);
//     aligned_free(work);
//
//     Cblacs_gridexit(ictxt);
//     printf("  gridexit\n");
// }


void my_peig(GTMatrix_t gtm_A, GTMatrix_t gtm_B, int n, int nprow, int npcol, double *eval) {
    int myrank, ictxt, myrow, mycol, lo[2], hi[2], ld;
    int descA[9], descZ[9], info, izero = 0, ione = 1, lwork, liwork;
    double t1, t2, *A, *Z, *work;
	
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int nb = min(n / nprow, n / npcol);
	
    Cblacs_get(-1, 0, &ictxt);
    Cblacs_gridinit(&ictxt, (char *)"Row-major", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int nrows = numroc_(&n, &nb, &myrow, &izero, &nprow);
    int ncols = numroc_(&n, &nb, &mycol, &izero, &npcol);
    int blocksize = nrows * ncols;
	
    A = (double *)aligned_malloc(blocksize * sizeof(double), 64);
    Z = (double *)aligned_malloc(blocksize * sizeof(double), 64);
    assert(A && Z);

    memset(A, 0, blocksize * sizeof(double));
    memset(Z, 0, blocksize * sizeof(double));

    // Arrays initialized for Schwartz eigenvalue
    descinit_(descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &nrows, &info);
    descinit_(descZ, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &nrows, &info);
    
    // Populate matrix A and adjust Z with hypothetical values(Matrix initialization, knowledge initialize data source)
    for (int i = 0; i < nrows; i++)
    	for (int j = 0; j < ncols; j++) {
            Z[i*ncols + j] = (double)rand()/RAND_MAX;  // test values ratio align equality random paradigm ensuring
        }
	
    // Determine optimal work size for WORK
    work = (double *)aligned_malloc(2 * sizeof(double), 64);
    int *iwork = (int *)aligned_malloc(2 * sizeof(int), 64);
    memset(work, 0, 2 * sizeof(double));
    memset(iwork, 0, 2 * sizeof(int));

    lwork = -1;
    liwork = -1;

    // Ensure right conversions added validate wrong conversion/accuracy ratio Rule Synchron.ensuring Z-Matrix valid<d*nx-k correlations>
    Cblacs_barrier(ictxt, "All");
    pdsyevd_("V", "U", &n, A, &ione, &ione, descA, eval, Z, &ione, &ione, descZ, work, &lwork, iwork, &liwork, &info);

    if (info != 0) {
        printf("[%d] Error In first pdsyevd_ execution: EXIT CODE=%d\n", myrank, info);
        MPI_Abort(MPI_COMM_WORLD, info);
    }

    lwork  = (int) work[0];
    liwork = (int)iwork[0];

	aligned_free(work);
    aligned_free(iwork);

    work  = (double*) aligned_malloc(lwork  * sizeof (double), 64);
    iwork = (int*)   aligned_malloc(liwork * sizeof (int), 64);
    assert(work && iwork);

    t1 = MPI_Wtime();
    pdsyevd_("V", "U", &n, A, &ione, &ione, descA, eval, Z, &ione, &ione, descZ, work, &lwork, iwork, &liwork, &info);

    t2 = MPI_Wtime();

	if (info != 0) {
        printf("[%d] Error In Second pdsyevd_ execution: EXIT CODE=%d\n", myrank, info);
        MPI_Abort(MPI_COMM_WORLD, info);
    } else if (!(myrank)) {
        printf("Done computing eigenvalues successfully. Computation duration= %.3lf secs\n", t2 - t1);
    }
    if (myrank == 0) printf("  pdsyev_ takes %.3lf secs\n", t2 - t1);

    // store desination matrix
    for (int i = 0; i < nrows; i++) 
    {
        for (int j = 0; j < ncols; j++)
            A[i * ncols + j] = Z[j * nrows + i];
    }
    
    GTM_startBatchPut(gtm_B);
    for (int i = 1; i <= nrows; i += nb) 
    {
        lo[0] = indxl2g_ (&i, &nb, &myrow, &izero, &nprow) - 1;
        hi[0] = lo[0] + nb - 1;
        hi[0] = hi[0] >= n ? n - 1 : hi[0];
        for (int j = 1; j <= ncols; j += nb) 
        {
            lo[1] = indxl2g_ (&j, &nb, &mycol, &izero, &npcol) - 1;
            hi[1] = lo[1] + nb - 1;
            hi[1] = hi[1] >= n ? n - 1 : hi[1];
            ld = ncols;
            GTM_addPutBlockRequest(
                gtm_B, 
                lo[0], hi[0] - lo[0] + 1,
                lo[1], hi[1] - lo[1] + 1,
                &(A[(i - 1) * ncols + j - 1]), ld
            );
        }
    }
    GTM_execBatchPut(gtm_B);
    GTM_stopBatchPut(gtm_B);
    GTM_sync(gtm_B);
    printf("aligned_free\n");
    aligned_free(A);
    aligned_free(Z);
    aligned_free(work);
    aligned_free(iwork);

    Cblacs_gridexit(ictxt);
}
