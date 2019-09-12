/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_KERNEL_WRAPPER_H
#define RSS_RINGOCCS_FRESNEL_KERNEL_WRAPPER_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>

#include "__fresnel_kernel.h"

/*  Functions from __fresnel_kernel.h                                         */
static void float_fresnel_psi(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    float *kD   = (float *)args[0];
    float *rho  = (float *)args[1];
    float *rho0 = (float *)args[2];
    float *phi  = (float *)args[3];
    float *phi0 = (float *)args[4];
    float *B    = (float *)args[5];
    float *D    = (float *)args[6];
    float *out  = (float *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_Psi_Float(kD[kD_i], rho[rho_i], rho0[rho0_i],
                                   phi[phi_i], phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void double_fresnel_psi(char **args, npy_intp *dimensions,
                               npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    double *kD   = (double *)args[0];
    double *rho  = (double *)args[1];
    double *rho0 = (double *)args[2];
    double *phi  = (double *)args[3];
    double *phi0 = (double *)args[4];
    double *B    = (double *)args[5];
    double *D    = (double *)args[6];
    double *out  = (double *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_Psi_Double(kD[kD_i], rho[rho_i], rho0[rho0_i],
                                    phi[phi_i], phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void long_double_fresnel_psi(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    long double *kD   = (long double *)args[0];
    long double *rho  = (long double *)args[1];
    long double *rho0 = (long double *)args[2];
    long double *phi  = (long double *)args[3];
    long double *phi0 = (long double *)args[4];
    long double *B    = (long double *)args[5];
    long double *D    = (long double *)args[6];
    long double *out  = (long double *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_Psi_Long_Double(kD[kD_i], rho[rho_i],
                                         rho0[rho0_i], phi[phi_i],
                                         phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}


static void float_fresnel_dpsi_dphi(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    float *kD   = (float *)args[0];
    float *rho  = (float *)args[1];
    float *rho0 = (float *)args[2];
    float *phi  = (float *)args[3];
    float *phi0 = (float *)args[4];
    float *B    = (float *)args[5];
    float *D    = (float *)args[6];
    float *out  = (float *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_dPsi_dPhi_Float(kD[kD_i], rho[rho_i],
                                         rho0[rho0_i], phi[phi_i],
                                         phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void double_fresnel_dpsi_dphi(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    double *kD   = (double *)args[0];
    double *rho  = (double *)args[1];
    double *rho0 = (double *)args[2];
    double *phi  = (double *)args[3];
    double *phi0 = (double *)args[4];
    double *B    = (double *)args[5];
    double *D    = (double *)args[6];
    double *out  = (double *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_dPsi_dPhi_Double(kD[kD_i], rho[rho_i],
                                          rho0[rho0_i], phi[phi_i],
                                          phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void long_double_fresnel_dpsi_dphi(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    long double *kD   = (long double *)args[0];
    long double *rho  = (long double *)args[1];
    long double *rho0 = (long double *)args[2];
    long double *phi  = (long double *)args[3];
    long double *phi0 = (long double *)args[4];
    long double *B    = (long double *)args[5];
    long double *D    = (long double *)args[6];
    long double *out  = (long double *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_dPsi_dPhi_Long_Double(kD[kD_i], rho[rho_i],
                                               rho0[rho0_i], phi[phi_i],
                                               phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void float_fresnel_d2psi_dphi2(char **args, npy_intp *dimensions,
                                      npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    float *kD   = (float *)args[0];
    float *rho  = (float *)args[1];
    float *rho0 = (float *)args[2];
    float *phi  = (float *)args[3];
    float *phi0 = (float *)args[4];
    float *B    = (float *)args[5];
    float *D    = (float *)args[6];
    float *out  = (float *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_d2Psi_dPhi2_Float(kD[kD_i], rho[rho_i],
                                           rho0[rho0_i], phi[phi_i],
                                           phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void double_fresnel_d2psi_dphi2(char **args, npy_intp *dimensions,
                                       npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    double *kD   = (double *)args[0];
    double *rho  = (double *)args[1];
    double *rho0 = (double *)args[2];
    double *phi  = (double *)args[3];
    double *phi0 = (double *)args[4];
    double *B    = (double *)args[5];
    double *D    = (double *)args[6];
    double *out  = (double *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_d2Psi_dPhi2_Double(kD[kD_i], rho[rho_i],
                                            rho0[rho0_i], phi[phi_i],
                                            phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void long_double_fresnel_d2psi_dphi2(char **args, npy_intp *dimensions,
                                            npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    long double *kD   = (long double *)args[0];
    long double *rho  = (long double *)args[1];
    long double *rho0 = (long double *)args[2];
    long double *phi  = (long double *)args[3];
    long double *phi0 = (long double *)args[4];
    long double *B    = (long double *)args[5];
    long double *D    = (long double *)args[6];
    long double *out  = (long double *)args[7];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_d2Psi_dPhi2_Long_Double(kD[kD_i], rho[rho_i],
                                                 rho0[rho0_i], phi[phi_i],
                                                 phi0[phi0_i], B[B_i], D[D_i]);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void float_fresnel_dpsi_dphi_ellipse(char **args, npy_intp *dimensions,
                                            npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    float *kD   =  (float *)args[0];
    float *rho  =  (float *)args[1];
    float *rho0 =  (float *)args[2];
    float *phi  =  (float *)args[3];
    float *phi0 =  (float *)args[4];
    float *B    =  (float *)args[5];
    float *D    =  (float *)args[6];
    float ecc   = *(float *)args[7];
    float peri  = *(float *)args[8];
    float *out  =  (float *)args[9];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_dPsi_dPhi_Ellipse_Float(kD[kD_i], rho[rho_i],
                                                 rho0[rho0_i], phi[phi_i],
                                                 phi0[phi0_i], B[B_i], D[D_i],
                                                 ecc, peri);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void double_fresnel_dpsi_dphi_ellipse(char **args, npy_intp *dimensions,
                                             npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    double *kD   =  (double *)args[0];
    double *rho  =  (double *)args[1];
    double *rho0 =  (double *)args[2];
    double *phi  =  (double *)args[3];
    double *phi0 =  (double *)args[4];
    double *B    =  (double *)args[5];
    double *D    =  (double *)args[6];
    double ecc   = *(double *)args[7];
    double peri  = *(double *)args[8];
    double *out  =  (double *)args[9];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_dPsi_dPhi_Ellipse_Double(kD[kD_i], rho[rho_i],
                                                  rho0[rho0_i], phi[phi_i],
                                                  phi0[phi0_i], B[B_i], D[D_i],
                                                  ecc, peri);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

static void long_double_fresnel_dpsi_dphi_ellipse(char **args,
                                                  npy_intp *dimensions,
                                                  npy_intp* steps, void* data)
{
    /*  Variables for indexing the various input arguments.                   */
    long i;
    long kD_i   = 0;
    long rho_i  = 0;
    long rho0_i = 0;
    long phi_i  = 0;
    long phi0_i = 0;
    long B_i    = 0;
    long D_i    = 0;

    /*  The number of elements in the output array.                           */
    long n = dimensions[0];

    long double *kD   =  (long double *)args[0];
    long double *rho  =  (long double *)args[1];
    long double *rho0 =  (long double *)args[2];
    long double *phi  =  (long double *)args[3];
    long double *phi0 =  (long double *)args[4];
    long double *B    =  (long double *)args[5];
    long double *D    =  (long double *)args[6];
    long double ecc   = *(long double *)args[7];
    long double peri  = *(long double *)args[8];
    long double *out  =  (long double *)args[9];

    /*  Since both arrays and numbers are allowed as inputs, we must make     *
     *  sure we index the variables appropriately. If an argument is only a   *
     *  single chunk of memory, this variable is a single number and we       *
     *  shouldn't try to index it. Otherwise, the input is an array and we    *
     *  should index over it. The variables below are either zero or one,     *
     *  corresponding to whether the argument is a number or an array.        */
    unsigned char kD_steps   = (steps[0] != 0);
    unsigned char rho_steps  = (steps[1] != 0);
    unsigned char rho0_steps = (steps[2] != 0);
    unsigned char phi_steps  = (steps[3] != 0);
    unsigned char phi0_steps = (steps[4] != 0);
    unsigned char B_steps    = (steps[5] != 0);
    unsigned char D_steps    = (steps[6] != 0);

    /*  Compute the Fresnel Psi function, looping over the variables.         */
    for (i = 0; i < n; i++) {
        out[i] = Fresnel_dPsi_dPhi_Ellipse_Long_Double(kD[kD_i], rho[rho_i],
                                                       rho0[rho0_i], phi[phi_i],
                                                       phi0[phi0_i], B[B_i],
                                                       D[D_i], ecc, peri);

        kD_i   += kD_steps;
        rho_i  += rho_steps;
        rho0_i += rho0_steps;
        phi_i  += phi_steps;
        phi0_i += phi0_steps;
        B_i    += B_steps;
        D_i    += D_steps;
    }
}

#endif