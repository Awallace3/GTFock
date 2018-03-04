#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <macdecls.h>
#include <sys/time.h>

#include "config.h"
#include "taskq.h"
#include "fock_task.h"


static inline void atomic_add_f64(volatile double* global_value, double addend)
{
    uint64_t expected_value, new_value;
    do {
        double old_value = *global_value;
        expected_value = _castf64_u64(old_value);
        new_value = _castf64_u64(old_value + addend);
    } while (!__sync_bool_compare_and_swap((volatile uint64_t*)global_value,
                                           expected_value, new_value));
}

/*
static void update_F(int num_dmat, double *integrals, int dimM, int dimN,
                    int dimP, int dimQ,
                    int flag1, int flag2, int flag3,
                    int iMN, int iPQ, int iMP, int iNP, int iMQ, int iNQ,
                    int iMP0, int iMQ0, int iNP0,
                    double **D1, double **D2, double **D3,
                    double *F_MN, double *F_PQ, double *F_NQ,
                    double *F_MP, double *F_MQ, double *F_NP,
                    int sizeX1, int sizeX2, int sizeX3,
                    int sizeX4, int sizeX5, int sizeX6,
                    int ldMN, int ldPQ, int ldNQ, int ldMP, int ldMQ, int ldNP)
{
    int flag4 = (flag1 == 1 && flag2 == 1) ? 1 : 0;
    int flag5 = (flag1 == 1 && flag3 == 1) ? 1 : 0;
    int flag6 = (flag2 == 1 && flag3 == 1) ? 1 : 0;
    int flag7 = (flag4 == 1 && flag3 == 1) ? 1 : 0;

    for (int i = 0 ; i < num_dmat; i++) {
        double *D_MN = D1[i] + iMN;
        double *D_PQ = D2[i] + iPQ;
        double *D_NQ = D3[i] + iNQ;
        double *D_MP = D3[i] + iMP0;
        double *D_MQ = D3[i] + iMQ0;
        double *D_NP = D3[i] + iNP0;    
        double *J_MN = &F_MN[i * sizeX1] + iMN;
        double *J_PQ = &F_PQ[i * sizeX2] + iPQ;
        double *K_NQ = &F_NQ[i * sizeX3] + iNQ;
        double *K_MP = &F_MP[i * sizeX4] + iMP;
        double *K_MQ = &F_MQ[i * sizeX5] + iMQ;
        double *K_NP = &F_NP[i * sizeX6] + iNP;
    
        for (int iN = 0; iN < dimN; iN++) {
            for (int iQ = 0; iQ < dimQ; iQ++) {
                int inq = iN * ldNQ + iQ;
                double k_NQ = 0;
                for (int iM = 0; iM < dimM; iM++) {
                    int imn = iM * ldMN + iN;
                    int imq = iM * ldMQ + iQ;
                    double j_MN = 0;
                    double k_MQ = 0;
                    for (int iP = 0; iP < dimP; iP++) {
                        int ipq = iP * ldPQ + iQ;
                        int imp = iM * ldMP + iP;
                        int inp = iN * ldNP + iP;
                        double I = 
                            integrals[iQ + dimQ*(iP + dimP * (iN + dimN * iM))];//Simint
                          //integrals[iM + dimM*(iN + dimN * (iP + dimP * iQ))];//OptERD
                        // F(m, n) += D(p, q) * 2 * I(m, n, p, q)
                        // F(n, m) += D(p, q) * 2 * I(n, m, p, q)
                        // F(m, n) += D(q, p) * 2 * I(m, n, q, p)
                        // F(n, m) += D(q, p) * 2 * I(n, m, q, p)
                        double vMN = 2.0 * (1 + flag1 + flag2 + flag4) *
                            D_PQ[iP * ldPQ + iQ] * I;
                        j_MN += vMN;
                        // F(p, q) += D(m, n) * 2 * I(p, q, m, n)
                        // F(p, q) += D(n, m) * 2 * I(p, q, n, m)
                        // F(q, p) += D(m, n) * 2 * I(q, p, m, n)
                        // F(q, p) += D(n, m) * 2 * I(q, p, n, m)
                        double vPQ = 2.0 * (flag3 + flag5 + flag6 + flag7) *
                            D_MN[iM * ldMN + iN] * I;
                        atomic_add_f64(&J_PQ[ipq], vPQ);
                        // F(m, p) -= D(n, q) * I(m, n, p, q)
                        // F(p, m) -= D(q, n) * I(p, q, m, n)
                        double vMP = (1 + flag3) *
                            1.0 * D_NQ[iN * ldNQ + iQ] * I;
                        atomic_add_f64(&K_MP[imp], -vMP);
                        // F(n, p) -= D(m, q) * I(n, m, p, q)
                        // F(p, n) -= D(q, m) * I(p, q, n, m)
                        double vNP = (flag1 + flag5) *
                            1.0 * D_MQ[iM * ldNQ + iQ] * I;
                        atomic_add_f64(&K_NP[inp], -vNP);
                        // F(m, q) -= D(n, p) * I(m, n, q, p)
                        // F(q, m) -= D(p, n) * I(q, p, m, n)
                        double vMQ = (flag2 + flag6) *
                            1.0 * D_NP[iN * ldNQ + iP] * I;
                        k_MQ -= vMQ;
                        // F(n, q) -= D(m, p) * I(n, m, q, p)
                        // F(q, n) -= D(p, m) * I(q, p, n, m)
                        double vNQ = (flag4 + flag7) *
                            1.0 * D_MP[iM * ldNQ + iP] * I;
                        k_NQ -= vNQ;
                    }
                    atomic_add_f64(&J_MN[imn], j_MN);
                    atomic_add_f64(&K_MQ[imq], k_MQ);
                }
                atomic_add_f64(&K_NQ[inq], k_NQ);
            }
        } // for (int iN = 0; iN < dimN; iN++)
    } // for (int i = 0 ; i < num_dmat; i++)
}
*/

double *update_F_buf  = NULL;
int update_F_buf_size = 0;

static void update_F(
    int tid, int num_dmat, double *integrals, int dimM, int dimN,
    int dimP, int dimQ,
    int flag1, int flag2, int flag3,
    int iMN, int iPQ, int iMP, int iNP, int iMQ, int iNQ,
    int iMP0, int iMQ0, int iNP0,
    double **D1, double **D2, double **D3,
    double *F_MN, double *F_PQ, double *F_NQ,
    double *F_MP, double *F_MQ, double *F_NP,
    int sizeX1, int sizeX2, int sizeX3,
    int sizeX4, int sizeX5, int sizeX6,
    int ldMN, int ldPQ, int ldNQ, int ldMP, int ldMQ, int ldNP
)
{
    int flag4 = (flag1 == 1 && flag2 == 1) ? 1 : 0;
    int flag5 = (flag1 == 1 && flag3 == 1) ? 1 : 0;
    int flag6 = (flag2 == 1 && flag3 == 1) ? 1 : 0;
    int flag7 = (flag4 == 1 && flag3 == 1) ? 1 : 0;
    
    double *K_NQ_buf = update_F_buf + 2 * update_F_buf_size * tid;
    double *J_PQ_buf = K_NQ_buf + update_F_buf_size;
    __assume_aligned(K_NQ_buf, 64);
    __assume_aligned(J_PQ_buf, 64);
    
    for (int i = 0 ; i < num_dmat; i++) {
        double *D_MN = D1[i] + iMN;
        double *D_PQ = D2[i] + iPQ;
        double *D_NQ = D3[i] + iNQ;
        double *D_MP = D3[i] + iMP0;
        double *D_MQ = D3[i] + iMQ0;
        double *D_NP = D3[i] + iNP0;    
        double *J_MN = &F_MN[i * sizeX1] + iMN;
        double *J_PQ = &F_PQ[i * sizeX2] + iPQ;
        double *K_NQ = &F_NQ[i * sizeX3] + iNQ;
        double *K_MP = &F_MP[i * sizeX4] + iMP;
        double *K_MQ = &F_MQ[i * sizeX5] + iMQ;
        double *K_NP = &F_NP[i * sizeX6] + iNP;
    
        for (int iN = 0; iN < dimN; iN++) 
        {
            memset(K_NQ_buf, 0, sizeof(double) * dimQ);
            for (int iP = 0; iP < dimP; iP++) 
            {
                memset(J_PQ_buf, 0, sizeof(double) * dimQ);
                int inp = iN * ldNP + iP;
                double k_NP = 0.0;
                double vMQ  = (flag2 + flag6) * D_NP[iN * ldNQ + iP];
                for (int iM = 0; iM < dimM; iM++) 
                {
                    int imn = iM * ldMN + iN;
                    int imp = iM * ldMP + iP;
                    
                    int Ibase  = dimQ * (iP + dimP * (iN + dimN * iM));
                    double vPQ = D_MN[iM * ldMN + iN] * 2.0 * (flag3 + flag5 + flag6 + flag7);
                    double vNQ = D_MP[iM * ldNQ + iP] * (flag4 + flag7);
                    
                    double j_MN = 0.0, k_MP = 0.0, k_NP0 = 0.0;
                    
                    #pragma simd
                    for (int iQ = 0; iQ < dimQ; iQ++) 
                    {
                        double I = integrals[iQ + Ibase];
                        
                        j_MN  += D_PQ[iP * ldPQ + iQ] * I;
                        k_MP  -= D_NQ[iN * ldNQ + iQ] * I;
                        k_NP0 -= D_MQ[iM * ldNQ + iQ] * I;
                        
                        J_PQ_buf[iQ] += vPQ * I;
                        K_NQ_buf[iQ] -= vNQ * I;
                    }
                    j_MN *= 2.0 * (1 + flag1 + flag2 + flag4);
                    k_MP *= (1 + flag3);
                    k_NP += k_NP0 * (flag1 + flag5);
                    
                    for (int iQ = 0; iQ < dimQ; iQ++)
                    {
                        int imq = iM * ldMQ + iQ;
                        double I = integrals[iQ + Ibase];
                        atomic_add_f64(&K_MQ[imq], -vMQ * I);
                    }
                    
                    atomic_add_f64(&J_MN[imn], j_MN);
                    atomic_add_f64(&K_MP[imp], k_MP);
                }
                atomic_add_f64(&K_NP[inp], k_NP);
                for (int iQ = 0; iQ < dimQ; iQ++)
                    atomic_add_f64(&J_PQ[iP * ldPQ + iQ], J_PQ_buf[iQ]);
            } // for (int iP = 0; iP < dimP; iP++) 
            for (int iQ = 0; iQ < dimQ; iQ++)
                atomic_add_f64(&K_NQ[iN * ldNQ + iQ], K_NQ_buf[iQ]);
        } // for (int iN = 0; iN < dimN; iN++)
    } // for (int i = 0 ; i < num_dmat; i++)
}

#include "thread_quartet_buf.h"

void update_F_with_KetShellPairList(
    int tid, int num_dmat, double *batch_integrals, int batch_nints, int npairs, 
    KetShellPairList_s *target_shellpair_list,
    double **D1, double **D2, double **D3,
    double *F_MN, double *F_PQ, double *F_NQ, double *F_MP, double *F_MQ, double *F_NP,
    int sizeX1, int sizeX2, int sizeX3, int sizeX4, int sizeX5, int sizeX6,
    int ldX1, int ldX2, int ldX3, int ldX4, int ldX5, int ldX6
)
{
    for (int ipair = 0; ipair < npairs; ipair++)
    {
        update_F(
            tid, num_dmat, &batch_integrals[ipair * batch_nints], 
            target_shellpair_list->fock_info_list[ipair].dimM, 
            target_shellpair_list->fock_info_list[ipair].dimN, 
            target_shellpair_list->fock_info_list[ipair].dimP, 
            target_shellpair_list->fock_info_list[ipair].dimQ,
            target_shellpair_list->fock_info_list[ipair].flag1,
            target_shellpair_list->fock_info_list[ipair].flag2, 
            target_shellpair_list->fock_info_list[ipair].flag3,
            target_shellpair_list->fock_info_list[ipair].iMN, 
            target_shellpair_list->fock_info_list[ipair].iPQ, 
            target_shellpair_list->fock_info_list[ipair].iMP, 
            target_shellpair_list->fock_info_list[ipair].iNP, 
            target_shellpair_list->fock_info_list[ipair].iMQ, 
            target_shellpair_list->fock_info_list[ipair].iNQ,
            target_shellpair_list->fock_info_list[ipair].iMP0, 
            target_shellpair_list->fock_info_list[ipair].iMQ0, 
            target_shellpair_list->fock_info_list[ipair].iNP0,
            D1, D2, D3,
            F_MN, F_PQ, F_NQ, F_MP, F_MQ, F_NP,
            sizeX1, sizeX2, sizeX3, sizeX4, sizeX5, sizeX6,
            ldX1, ldX2, ldX3, ldX4, ldX5, ldX6
        );
    }
}

// for SCF, J = K
void fock_task(BasisSet_t basis, SIMINT_t simint, int ncpu_f, int num_dmat,
               int *shellptr, double *shellvalue,
               int *shellid, int *shellrid, int *f_startind,
               int *rowpos, int *colpos, int *rowptr, int *colptr,
               double tolscr2, int startrow, int startcol,
               int startM, int endM, int startP, int endP,
               double **D1, double **D2, double **D3,
               double *F1, double *F2, double *F3,
               double *F4, double *F5, double *F6, 
               int ldX1, int ldX2, int ldX3,
               int ldX4, int ldX5, int ldX6,
               int sizeX1, int sizeX2, int sizeX3,
               int sizeX4, int sizeX5, int sizeX6,
               double *nitl, double *nsq)              
{
    int startMN = shellptr[startM];
    int endMN = shellptr[endM + 1];
    int startPQ = shellptr[startP];
    int endPQ = shellptr[endP + 1];
    
    if (update_F_buf_size == 0)
    {
        // startPQ & endPQ & f_startind remains unchanged for each call
        // So just allocate the buffer once
        for (int j = startPQ; j < endPQ; j++) 
        {
            int Q = shellid[j];
            int dimQ = f_startind[Q + 1] - f_startind[Q];
            if (dimQ > update_F_buf_size) update_F_buf_size = dimQ;
        }
        update_F_buf_size = (update_F_buf_size + 8) / 8 * 8;  // Align to 64 bytes
        int nthreads = omp_get_max_threads();
        update_F_buf = _mm_malloc(sizeof(double) * nthreads * 2 * update_F_buf_size, 64);
        assert(update_F_buf != NULL);
    }
    
    #pragma omp parallel
    {
        // init    
        int nt = omp_get_thread_num();
        int nf = nt/ncpu_f;
        double *F_MN = &(F1[nf * sizeX1 * num_dmat]);
        double *F_PQ = &(F2[nf * sizeX2 * num_dmat]);
        double *F_NQ = F3;
        double *F_MP = &(F4[nf * sizeX4 * num_dmat]);
        double *F_MQ = &(F5[nf * sizeX5 * num_dmat]);
        double *F_NP = &(F6[nf * sizeX6 * num_dmat]);
        double mynsq = 0.0;
        double mynitl = 0.0;        
        
        // Pending quartets that need to be computed
        ThreadQuartetLists_s *thread_quartet_lists = (ThreadQuartetLists_s*) malloc(sizeof(ThreadQuartetLists_s));
        init_ThreadQuartetLists(thread_quartet_lists);
        
        // Simint shell and multi_shellpair buffer for batch computation
        void *thread_shell_buf, *thread_multi_shellpairs;
        CInt_SIMINT_createThreadShellBuf(&thread_shell_buf);
        CInt_SIMINT_createThreadMultishellpairs(&thread_multi_shellpairs);
        
        #pragma omp for schedule(dynamic) 
        for (int i = startMN; i < endMN; i++) 
        {
            int M = shellrid[i];
            int N = shellid[i];
            
            reset_ThreadQuartetLists(thread_quartet_lists, M, N);
            
            double value1 = shellvalue[i];            
            int dimM = f_startind[M + 1] - f_startind[M];
            int dimN = f_startind[N + 1] - f_startind[N];
            int iX1M = f_startind[M] - f_startind[startrow];
            int iX3M = rowpos[M]; 
            int iXN  = rowptr[i];
            int iMN  = iX1M * ldX1+ iXN;
            int flag1 = (value1 < 0.0) ? 1 : 0;   
            for (int j = startPQ; j < endPQ; j++) 
            {
                int P = shellrid[j];
                int Q = shellid[j];
                if ((M > P && (M + P) % 2 == 1) || 
                    (M < P && (M + P) % 2 == 0))
                    continue;                
                if ((M == P) &&
                    ((N > Q && (N + Q) % 2 == 1) ||
                    (N < Q && (N + Q) % 2 == 0)))
                    continue;
                double value2 = shellvalue[j];
                int dimP = f_startind[P + 1] - f_startind[P];
                int dimQ = f_startind[Q + 1] - f_startind[Q];
                int iX2P = f_startind[P] - f_startind[startcol];
                int iX3P = colpos[P];
                int iXQ  = colptr[j];               
                int iPQ  = iX2P * ldX2 + iXQ;                             
                int iNQ  = iXN  * ldX3 + iXQ;                
                int iMP  = iX1M * ldX4 + iX2P;
                int iMQ  = iX1M * ldX5 + iXQ;
                int iNP  = iXN  * ldX6 + iX2P;
                int iMP0 = iX3M * ldX3 + iX3P;
                int iMQ0 = iX3M * ldX3 + iXQ;
                int iNP0 = iXN  * ldX3 + iX3P;               
                int flag3 = (M == P && Q == N) ? 0 : 1;                    
                int flag2 = (value2 < 0.0) ? 1 : 0;
                if (fabs(value1 * value2) >= tolscr2) 
                {
                    mynsq += 1.0;
                    mynitl += dimM*dimN*dimP*dimQ;              
                    
                    int am_pair_index = CInt_SIMINT_getShellpairAMIndex(simint, P, Q);
                    
                    KetShellPairList_s *target_shellpair_list = &thread_quartet_lists->ket_shellpair_lists[am_pair_index];
                    
                    // Save this shell pair to the target ket shellpair list
                    int add_KetShellPair_ret = add_KetShellPair(
                        target_shellpair_list, P, Q,
                        dimM, dimN, dimP, dimQ, 
                        flag1, flag2, flag3,
                        iMN, iPQ, iMP, iNP, iMQ, iNQ,
                        iMP0, iMQ0, iNP0
                    );
                    assert(add_KetShellPair_ret == 1);
                    
                    // Target ket shellpair list is full, handles it
                    if (target_shellpair_list->num_shellpairs == _SIMINT_NSHELL_SIMD) 
                    {
                        int npairs = target_shellpair_list->num_shellpairs;
                        double *thread_batch_integrals;
                        int thread_batch_nints;
                        
                        CInt_computeShellQuartetBatch_SIMINT(
                            basis, simint, nt,
                            thread_quartet_lists->M, 
                            thread_quartet_lists->N, 
                            target_shellpair_list->P_list,
                            target_shellpair_list->Q_list,
                            npairs, &thread_batch_integrals, &thread_batch_nints,
                            &thread_shell_buf, 
                            &thread_multi_shellpairs
                        );
                        
                        if (thread_batch_nints > 0)
                        {
                            update_F_with_KetShellPairList(
                                nt, num_dmat, thread_batch_integrals, thread_batch_nints, npairs, 
                                target_shellpair_list,
                                D1, D2, D3,
                                F_MN, F_PQ, F_NQ, F_MP, F_MQ, F_NP,
                                sizeX1, sizeX2, sizeX3, sizeX4, sizeX5, sizeX6,
                                ldX1, ldX2, ldX3, ldX4, ldX5, ldX6
                            );
                        }
                        
                        // Ket shellpair list is processed, reset it
                        reset_KetShellPairList(target_shellpair_list);
                    }
                }  // if (fabs(value1 * value2) >= tolscr2) 
            }  // for (int j = startPQ; j < endPQ; j++)
            
            // Process all the remaining shell pairs in the thread's list
            for (int am_pair_index = 0; am_pair_index < _SIMINT_AM_PAIRS; am_pair_index++)
            {
                KetShellPairList_s *target_shellpair_list = &thread_quartet_lists->ket_shellpair_lists[am_pair_index];
                
                if (target_shellpair_list->num_shellpairs > 0)  // Ket shellpair list is not empty, handles it
                {
                    int npairs = target_shellpair_list->num_shellpairs;
                    double *thread_batch_integrals;
                    int thread_batch_nints;
                    
                    CInt_computeShellQuartetBatch_SIMINT(
                        basis, simint, nt,
                        thread_quartet_lists->M, 
                        thread_quartet_lists->N, 
                        target_shellpair_list->P_list,
                        target_shellpair_list->Q_list,
                        npairs, &thread_batch_integrals, &thread_batch_nints,
                        &thread_shell_buf, 
                        &thread_multi_shellpairs
                    );
                    
                    if (thread_batch_nints > 0)
                    {
                        update_F_with_KetShellPairList(
                            nt, num_dmat, thread_batch_integrals, thread_batch_nints, npairs, 
                            target_shellpair_list,
                            D1, D2, D3,
                            F_MN, F_PQ, F_NQ, F_MP, F_MQ, F_NP,
                            sizeX1, sizeX2, sizeX3, sizeX4, sizeX5, sizeX6,
                            ldX1, ldX2, ldX3, ldX4, ldX5, ldX6
                        );
                    }
                    
                    // Ket shellpair list is processed, reset it
                    reset_KetShellPairList(target_shellpair_list);
                }
            }
        }  // for (int i = startMN; i < endMN; i++)

        #pragma omp critical
        {
            *nitl += mynitl;
            *nsq += mynsq;
        }
        
        CInt_SIMINT_freeThreadShellBuf(&thread_shell_buf);
        CInt_SIMINT_freeThreadMultishellpairs(&thread_multi_shellpairs);
        free_ThreadQuartetLists(thread_quartet_lists);
        
        free(thread_quartet_lists);
    } /* #pragma omp parallel */
}


void reset_F(int numF, int num_dmat, double *F1, double *F2, double *F3,
             double *F4, double *F5, double *F6,
             int sizeX1, int sizeX2, int sizeX3,
             int sizeX4, int sizeX5, int sizeX6)
{
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int k = 0; k < numF * sizeX1 * num_dmat; k++) {
            F1[k] = 0.0;    
        }
        #pragma omp for nowait
        for (int k = 0; k < numF * sizeX2 * num_dmat; k++) {
            F2[k] = 0.0;
        }
        #pragma omp for nowait
        for (int k = 0; k < sizeX3 * num_dmat; k++) {
            F3[k] = 0.0;
        }
        #pragma omp for nowait
        for (int k = 0; k < numF * sizeX4 * num_dmat; k++) {
            F4[k] = 0.0;
        }
        #pragma omp for nowait
        for (int k = 0; k < numF * sizeX5 * num_dmat; k++) {
            F5[k] = 0.0;
        }
        #pragma omp for nowait
        for (int k = 0; k < numF * sizeX6 * num_dmat; k++) {
            F6[k] = 0.0;
        }
    }
}


void reduce_F(int numF, int num_dmat,
              double *F1, double *F2, double *F3,
              double *F4, double *F5, double *F6,
              int sizeX1, int sizeX2, int sizeX3,
              int sizeX4, int sizeX5, int sizeX6,
              int maxrowsize, int maxcolsize,
              int maxrowfuncs, int maxcolfuncs,
              int iX3M, int iX3P,
              int ldX3, int ldX4, int ldX5, int ldX6)
{
    #pragma omp parallel
    {
        #pragma omp for
        for (int k = 0; k < sizeX1 * num_dmat; k++) {
            for (int p = 1; p < numF; p++) {
                F1[k] += F1[k + p * sizeX1 * num_dmat];
            }
        }
        #pragma omp for
        for (int k = 0; k < sizeX2 * num_dmat; k++) {
            for (int p = 1; p < numF; p++) {
                F2[k] += F2[k + p * sizeX2 * num_dmat];
            }
        }
        #pragma omp for
        for (int k = 0; k < sizeX4 * num_dmat; k++) {
            for (int p = 1; p < numF; p++) {
                F4[k] += F4[k + p * sizeX4 * num_dmat];   
            }
        }
        #pragma omp for
        for (int k = 0; k < sizeX5 * num_dmat; k++) {
            for (int p = 1; p < numF; p++) {
                F5[k] += F5[k + p * sizeX5 * num_dmat];   
            }
        }
        #pragma omp for
        for (int k = 0; k < sizeX6 * num_dmat; k++) {
            for (int p = 1; p < numF; p++) {
                F6[k] += F6[k + p * sizeX6 * num_dmat];   
            }
        }

        int iMP = iX3M * ldX3 + iX3P;
        int iMQ = iX3M * ldX3;
        int iNP = iX3P;
        for (int k = 0; k < num_dmat; k++) {
            #pragma omp for
            for (int iM = 0; iM < maxrowfuncs; iM++) {
                for (int iP = 0; iP < maxcolfuncs; iP++) {
                    F3[k * sizeX3 + iMP + iM * ldX3 + iP]
                        += F4[k * sizeX4 + iM * ldX4 + iP];
                }
                for (int iQ = 0; iQ < maxcolsize; iQ++) {
                    F3[k * sizeX3 + iMQ + iM * ldX3 + iQ] +=
                        F5[k * sizeX5 + iM * ldX5 + iQ];    
                }
            }
            #pragma omp for
            for (int iN = 0; iN < maxrowsize; iN++) {
                for (int iP = 0; iP < maxcolfuncs; iP++) {
                    F3[k * sizeX3 + iNP + iN * ldX3 + iP] +=
                        F6[k * sizeX6 + iN * ldX6 + iP];
                }
            }
        }
    }
}
