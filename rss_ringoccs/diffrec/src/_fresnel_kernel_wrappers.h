/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_KERNEL_WRAPPER_H
#define RSS_RINGOCCS_FRESNEL_KERNEL_WRAPPER_H

/*  Where npy_intp is defined.                                                */
#include <numpy/ndarraytypes.h>

#include "__fresnel_kernel.h"

/*  Functions from __fresnel_kernel.h                                         */
static void float_fresnel_psi(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data){
    npy_intp i;
    npy_intp n = dimensions[0];

    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((float *)out) = Fresnel_Psi_Float(*(float *)kD,   *(float *)rho,
                                            *(float *)rho0, *(float *)phi,
                                            *(float *)phi0, *(float *)B,
                                            *(float *)D);

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

static void double_fresnel_psi(char **args, npy_intp *dimensions,
                               npy_intp* steps, void* data){
    npy_intp i;
    npy_intp n = dimensions[0];

    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((double *)out) = Fresnel_Psi_Double(*(double *)kD, *(double *)rho,
                                              *(double *)rho0, *(double *)phi,
                                              *(double *)phi0, *(double *)B,
                                              *(double *)D);

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

static void long_double_fresnel_psi(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data){
    npy_intp i;
    npy_intp n = dimensions[0];

    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((long double *)out) = Fresnel_Psi_Double(*(long double *)kD,
                                                   *(long double *)rho,
                                                   *(long double *)rho0,
                                                   *(long double *)phi,
                                                   *(long double *)phi0,
                                                   *(long double *)B,
                                                   *(long double *)D);

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

static void float_fresnel_dpsi_dphi(char **args, npy_intp *dimensions,
                                    npy_intp* steps, void* data){
    npy_intp i;
    npy_intp n = dimensions[0];
    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((float *)out) = Fresnel_dPsi_dPhi_Float(
            *(float *)kD, *(float *)rho, *(float *)rho0, *(float *)phi,
            *(float *)phi0, *(float *)B, *(float *)D
        );

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

static void double_fresnel_dpsi_dphi(char **args, npy_intp *dimensions,
                                     npy_intp* steps, void* data){
    npy_intp i;
    npy_intp n = dimensions[0];
    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((double *)out) = Fresnel_dPsi_dPhi_Double(
            *(double *)kD, *(double *)rho, *(double *)rho0, *(double *)phi,
            *(double *)phi0, *(double *)B, *(double *)D
        );

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

static void long_double_fresnel_dpsi_dphi(char **args, npy_intp *dimensions,
                                          npy_intp* steps, void* data){
    npy_intp i;
    npy_intp n = dimensions[0];
    char *kD   = args[0];
    char *rho  = args[1];
    char *rho0 = args[2];
    char *phi  = args[3];
    char *phi0 = args[4];
    char *B    = args[5];
    char *D    = args[6];
    char *out  = args[7];

    npy_intp kD_steps   = steps[0];
    npy_intp rho_steps  = steps[1];
    npy_intp rho0_steps = steps[2];
    npy_intp phi_steps  = steps[3];
    npy_intp phi0_steps = steps[4];
    npy_intp B_steps    = steps[5];
    npy_intp D_steps    = steps[6];
    npy_intp out_steps  = steps[7];

    for (i = 0; i < n; i++) {
        *((long double *)out) = Fresnel_dPsi_dPhi_Long_Double(
            *(long double *)kD, *(long double *)rho, *(long double *)rho0,
            *(long double *)phi, *(long double *)phi0, *(long double *)B,
            *(long double *)D
        );

        kD   += kD_steps;
        rho  += rho_steps;
        rho0 += rho0_steps;
        phi  += phi_steps;
        phi0 += phi0_steps;
        B    += B_steps;
        D    += D_steps;
        out  += out_steps;
    }
}

#endif