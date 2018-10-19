#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
//#include <ga.h>
//#include <macdecls.h>
#include <sys/time.h>

#include "config.h"
#include "taskq.h"
#include "fock_buf.h"

#include "Buzz_Matrix.h"

// NOTICE: load_full_DenMat() and store_local_bufF() needs that num_dmat2==1

void load_full_DenMat(PFock_t pfock)
{
    Buzz_startBatchGet(pfock->bm_Dmat);
    Buzz_addGetBlockRequest(pfock->bm_Dmat, 0, pfock->nbf, 0, pfock->nbf, pfock->D_mat, pfock->nbf);
    Buzz_execBatchGet(pfock->bm_Dmat);
    Buzz_stopBatchGet(pfock->bm_Dmat);
    Buzz_Sync(pfock->bm_Dmat);
}

void store_local_bufF(PFock_t pfock)
{
    int *loadrow = pfock->loadrow;
    int *loadcol = pfock->loadcol;
    int sizerow = pfock->sizeloadrow;
    int sizecol = pfock->sizeloadcol;
    int myrank;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    int lo[2];
    int hi[2];
    
    Buzz_Matrix_t bm_J = pfock->bm_Fmat;
    #ifdef __SCF__
    Buzz_Matrix_t bm_K = pfock->bm_Fmat;
    #else
    Buzz_Matrix_t bm_K = pfock->bm_Kmat;
    #endif
    
    lo[0] = myrank;
    hi[0] = myrank;
    lo[1] = 0;
    
    // local buffers
    int ldF1 = pfock->ldX1;
    int ldF2 = pfock->ldX2;
    int ldF3 = pfock->ldX3;    
    double *F1 = pfock->bm_F1->mat_block;
    double *F2 = pfock->bm_F2->mat_block;
    double *F3 = pfock->bm_F3->mat_block;
    
    Buzz_startBatchUpdate(bm_J);
    
    // update F1
    lo[0] = pfock->sfunc_row;
    hi[0] = pfock->efunc_row;
    for (int A = 0; A < sizerow; A++) 
    {
        lo[1] = loadrow[PLEN * A + P_LO];
        hi[1] = loadrow[PLEN * A + P_HI];
        int posrow = loadrow[PLEN * A + P_W];
        
        Buzz_addAccumulateBlockRequest(
            bm_J, 
            lo[0], hi[0] - lo[0] + 1,
            lo[1], hi[1] - lo[1] + 1,
            F1 + posrow, ldF1
        );
    }

    // update F2
    lo[0] = pfock->sfunc_col;
    hi[0] = pfock->efunc_col;
    for (int B = 0; B < sizecol; B++) 
    {
        lo[1] = loadcol[PLEN * B + P_LO];
        hi[1] = loadcol[PLEN * B + P_HI];
        int poscol = loadcol[PLEN * B + P_W];
        
        Buzz_addAccumulateBlockRequest(
            bm_J, 
            lo[0], hi[0] - lo[0] + 1,
            lo[1], hi[1] - lo[1] + 1,
            F2 + poscol, ldF2
        );
    }

    Buzz_execBatchUpdate(bm_J);
    Buzz_stopBatchUpdate(bm_J);
    Buzz_Sync(bm_J);
    
    // update F3
    Buzz_startBatchUpdate(bm_K);
    for (int A = 0; A < sizerow; A++) 
    {
        lo[0] = loadrow[PLEN * A + P_LO];
        hi[0] = loadrow[PLEN * A + P_HI];
        int posrow = loadrow[PLEN * A + P_W];
        for (int B = 0; B < sizecol; B++) 
        {
            lo[1] = loadcol[PLEN * B + P_LO];
            hi[1] = loadcol[PLEN * B + P_HI];
            int poscol = loadcol[PLEN * B + P_W];
            
            Buzz_addAccumulateBlockRequest(
                bm_K, 
                lo[0], hi[0] - lo[0] + 1,
                lo[1], hi[1] - lo[1] + 1,
                F3 + posrow * ldF3 + poscol, ldF3
            );
        }
    }
    Buzz_execBatchUpdate(bm_K);
    Buzz_stopBatchUpdate(bm_K);
    Buzz_Sync(bm_K);
}


void compute_FD_ptr(PFock_t pfock, int startM, int endM, int *ptrrow, int *rowsize)
{
    for (int A = 0; A < pfock->nshells; A++) {
        ptrrow[A] = -1;
    }    
    // init row pointers
    for (int A = startM; A <= endM; A++) {
        int start = pfock->shellptr[A];
        int end = pfock->shellptr[A + 1]; 
        for (int i = start; i < end; i++) {
            int B = pfock->shellid[i];
            ptrrow[B] = 1;
        }
    }
    for (int i = 0; i < pfock->natoms; i++)
    {
        int start = pfock->s_startind[i];
        int end = pfock->s_startind[i + 1];
        int flag = -1;
        for (int A = start; A < end; A++)
        {
            if (ptrrow[A] != -1)
                flag = 1;
        }
        for (int A = start; A < end; A++)
        {
            ptrrow[A] = flag;
        }
    }
    *rowsize = 0;
    for (int A = 0; A < pfock->nshells; A++)
    {
        if (ptrrow[A] == 1)
        {
            ptrrow[A] = *rowsize;           
            *rowsize += pfock->f_startind[A + 1] - pfock->f_startind[A];
        }
    }
}


void init_FD_load(PFock_t pfock, int *ptrrow, int **loadrow, int *loadsize)
{    
    int loadcount = 0;
    for (int A = 0; A < pfock->nshells; A++) {
        if (ptrrow[A] != -1) {
            while (A < pfock->nshells && ptrrow[A] != -1) {
                A++;
            }           
            loadcount++;
        }
    }
    *loadrow = (int *)PFOCK_MALLOC(sizeof(int) * PLEN * loadcount);
    assert(NULL != *loadrow);
    *loadsize = loadcount;
    
    loadcount = 0;
    for (int A = 0; A < pfock->nshells; A++) {
        int idx = ptrrow[A];
        if (idx != -1) {
            int lo = pfock->f_startind[A];
            while (A < pfock->nshells && ptrrow[A] != -1) {
                A++;
            }           
            int hi = pfock->f_startind[A] - 1;
            (*loadrow)[loadcount * PLEN + P_LO] = lo;
            (*loadrow)[loadcount * PLEN + P_HI] = hi;
            (*loadrow)[loadcount * PLEN + P_W] = idx;
            loadcount++;
        }
    }
}
