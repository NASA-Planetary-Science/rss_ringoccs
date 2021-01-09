/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************
 *                             Diffrec Module                                 *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains the C-Python API code needed to use the C code     *
 *      contained in the diffrec/ directory within Python. It contains        *
 *      wrappers for several functions and the DiffractionCorrection class,   *
 *      which is the primary class for performing diffraction correction.     *
 *      Upon successful compilation there should be a diffrec.*info*.so file. *
 *      The Python interpreter can directly import diffrec if this file is in *
 *      its path via:                                                         *
 *          import diffrec                                                    *
 *      More likely, the user will be at the top-level of rss_ringoccs, where *
 *      one can simply run the following:                                     *
 *          from rss_ringoccs import diffrec                                  *
 *      The following are defined in this file:                               *
 *          DiffractionCorrection:                                            *
 *              A Python class, written in C using the C-Python API. The      *
 *              equivalent of the __init__ function at the python level is    *
 *              defined via Diffrec_init. This is the main class for          *
 *              processing diffracted data. See the file                      *
 *              _diffraction_correction.c and the associated header file      *
 *              _diffraction_correction.h for more information on the         *
 *              algorithms implemented.                                       *
 *          gap_diffraction:                                                  *
 *              A diffraction modeling tool which returns the diffraction     *
 *              profile of a gap, which has 0 complex amplitude everywhere    *
 *              except the gap where it evaluates to 1.                       *
 *          ringlet_diffraction:                                              *
 *              Similar to gap diffraction with the values reversed.          *
 *          right_straightedge:                                               *
 *              diffraction modeling of a "straight-edge". This has value     *
 *              1 for everything to the right of the straight edge, and 0 to  *
 *              the left.                                                     *
 *          left_straightedge:                                                *
 *              Like right_straightedge with the values reverse.              *
 *          fresnel_psi:                                                      *
 *              The Fresnel kernel function, as defined in MTR86.             *
 *          fresnel_scale:                                                    *
 *              The Fresnel scale, also defined in MTR86.                     *
 *          fresnel_dpsi_dphi:                                                *
 *              The partial derivative of the Fresnel kernel with respect to  *
 *              the phi variable.                                             *
 *          fresnel_dphi_dphi_ellipse:                                        *
 *              The partial derivative of the Fresnel kernel with respect to  *
 *              the phi variable where the kernel has been perturbed by       *
 *              the elliptical geometry of the rings.                         *
 *          fresnel_d2phi_dphi2:                                              *
 *              The second partial derivative of the Fresnel kernel with      *
 *              respect to the phi variable.                                  *
 *          single_slit_diffraction:                                          *
 *              Fraunhofer diffraction modeling of a single slit.             *
 *          double_slit_diffraction:                                          *
 *              Fraunhofer diffraction modeling of a double slit.             *
 *          square_wave_diffraction:                                          *
 *              Fresnel diffraction modeling of a square wave.                *
 *          square_well_phase:                                                *
 *              The modeled phase of the diffraction profile from a ringlet.  *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  The following functions are defined and used throughout, but are not      *
 *  available to the Python interpreter.                                      *
 *  GetNormeqFromString:                                                      *
 *      This acts like a Python dictionary where we input a string and        *
 *      retrieve a predefined value. It gives us the normalized equivalent    *
 *      width corresponding to a legal string.                                *
 *      NOTE:                                                                 *
 *          If you want to add your own window function, and your own         *
 *          corresponding normalized equivalent width, you will also need to  *
 *          add the corresponding string (the name of the window) to several  *
 *          of the error checking routines throughout the code.               *
 *  GetRangeFromString:                                                       *
 *      This also acts like a Python dictionary, giving us two values, the    *
 *      start and end points, corresponding to a given legal string.          *
 *      NOTE:                                                                 *
 *          Currently, all legal strings correspond to features of Saturn.    *
 *          If you wish to add your own string/feature, you will need to      *
 *          append it to this function and the error checking routines.       *
 *  capsule_cleanup:                                                          *
 *      A very important function, even though it is only two lines. This     *
 *      function ensures that memory we have allocated for variables in C     *
 *      that are then passed to the Python interpreter will still be freed    *
 *      when Python no longer needs it. Without this we will have serious     *
 *      memory leaks. It is primiarly used for safely returning numpy arrays. *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  All of the functions are declared "static" and cannot be used outside of  *
 *  this file. The only purpose of this file is generating the diffrec module.*
 *                                                                            *
 *  This code uses complex numbers throughout, and is compatible with the C99 *
 *  standard. To use this code, make sure your compiler supports C99 or more  *
 *  recent standards of the C Programming Language.                           *
 *                                                                            *
 *  In addition, for full use of all the routines one will need the FFTW      *
 *  library. Installation instruction can be found in the rss_ringoccs PDF.   *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  For successful compiling you will need the following NON-STANDARD         *
 *  header files / libraries:                                                 *
 *  libfftw3:                                                                 *
 *      The shared library file for the FFTW package. This is needed to       *
 *      use the FFT options of the DiffractionCorrection class.               *
 *  Python.h:                                                                 *
 *      The main header file for Python. If this code is compiled via the     *
 *      setup.py script provided, then Python.h will automatically be         *
 *      included correctly.                                                   *
 *  numpy/ndarraytypes.h:                                                     *
 *  numpy/ufuncobject.h:                                                      *
 *      Header files that allow us to use the Numpy API functions for         *
 *      accessing and creating numpy arrays at the C level. To get the        *
 *      location of these files you will need to use numpy's get_include()    *
 *      function. The setup.py script does this for you.                      *
 *  librssringoccsdiffractioncorrection.so:                                   *
 *      A non-standard library written and maintained by the rss_ringoccs     *
 *      project, it provide all of the diffraction correction algorithms      *
 *      written at the C level. The setup.py file will automatically link     *
 *      this program to the correct headers. To use this module Python will   *
 *      need to be able to find this .so. The setup scripts will install it   *
 *      in /usr/local/lib which is a standard location on Linux and MacOS.    *
 *  librssringoccsspecialfunctions.so:                                        *
 *      Another non-standard library maintained by rss_ringoccs, it contains  *
 *      a plethora of mathematical functions that are used in diffraction     *
 *      correction such as bessel functions, Legendre polynomials, and many   *
 *      more. The same criterion for librssringoccsdiffractioncorrection must *
 *      be satisfied in order to correctly link and load this library.        *
 *      This file is required implicitly by the                               *
 *      librssringoccsdiffractioncorrection.so library. We do not directly    *
 *      make calls to its functions.                                          *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************
 *                                History                                     *
 ******************************************************************************
 *  2020/09/01 (Ryan Maguire):                                                *
 *      Added DiffractionCorrection class in C. Updated comments.             *
 ******************************************************************************/

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*  The standard library header stdlib contains malloc, calloc, and realloc,  *
 *  free, as well as strtol (string-to-long) which converts a string like     *
 *  "314" to the integer 314.                                                 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/*  The following are NON-STANDARD header files that MUST BE IN YOUR PATH.    *
 *  If you installed python using anaconda then Python.h should automatically *
 *  be included in your path. Also, if you are using the setup.py script      *
 *  provided then inclusion of these files should be done for you.            */
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  The following header files are NON-STANDARD, and are a part of the        *
 *  rss_ringoccs package. The setup scripts will add the correct CFLAGS so    *
 *  compiler should see these without user interaction.                       */

/*  This contains the PyDiffrecObj, using the C-Python API, which allows us   *
 *  to define the DiffractionCorrection class in C.                           */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#include "rss_ringoccs_py_api.h"
#include "rss_ringoccs_Py_DLP_to_C_DLP.c"
#include "rss_ringoccs_C_Tau_to_Py_Tau.c"
#include "rss_ringoccs_Get_Py_Perturb.c"
#include "rss_ringoccs_Get_Py_Range.c"
#include "rss_ringoccs_Get_Py_Vars_From_Self.c"

/*  Deallocating function for the DiffractionCorrection class.                */
static void Diffrec_dealloc(PyDiffrecObj *self)
{
    Py_XDECREF(self->B_rad_vals);
    Py_XDECREF(self->D_km_vals);
    Py_XDECREF(self->F_km_vals);
    Py_XDECREF(self->f_sky_hz_vals);
    Py_XDECREF(self->p_norm_fwd_vals);
    Py_XDECREF(self->p_norm_vals);
    Py_XDECREF(self->phase_fwd_vals);
    Py_XDECREF(self->phase_rad_vals);
    Py_XDECREF(self->phase_vals);
    Py_XDECREF(self->phi_rad_vals);
    Py_XDECREF(self->phi_rl_rad_vals);
    Py_XDECREF(self->power_vals);
    Py_XDECREF(self->raw_tau_threshold_vals);
    Py_XDECREF(self->rev_info);
    Py_XDECREF(self->rho_corr_pole_km_vals);
    Py_XDECREF(self->rho_corr_timing_km_vals);
    Py_XDECREF(self->rho_dot_kms_vals);
    Py_XDECREF(self->rho_km_vals);
    Py_XDECREF(self->t_oet_spm_vals);
    Py_XDECREF(self->t_ret_spm_vals);
    Py_XDECREF(self->t_set_spm_vals);
    Py_XDECREF(self->tau_threshold_vals);
    Py_XDECREF(self->tau_vals);
    Py_XDECREF(self->w_km_vals);
    Py_XDECREF(self->dathist);
    Py_XDECREF(self->history);
    Py_XDECREF(self->T_hat_vals);
    Py_XDECREF(self->T_hat_fwd_vals);
    Py_XDECREF(self->T_vals);
    Py_XDECREF(self->rx_km_vals);
    Py_XDECREF(self->ry_km_vals);
    Py_XDECREF(self->rz_km_vals);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

/*  The init function for the dirrection correction class. This is the        *
 *  equivalent of the __init__ function defined in a normal python class.     */
static int Diffrec_init(PyDiffrecObj *self, PyObject *args, PyObject *kwds)
{
    /*  Declare variables for a DLP and Tau object.                           */
    rssringoccs_DLPObj *dlp;
    rssringoccs_TAUObj *tau;

    /*  For computing the calculation time.                                   */
    clock_t t1 = clock();
    clock_t t2;

    /*  The list of the keywords accepted by the DiffractionCorrection class. *
     *  dlp and res are REQUIRED inputs, the rest are optional. If the user   *
     *  does not provide these optional keywords, we must set them ourselves. */
    static char *kwlist[] = {
        "dlp",
        "res",
        "rng",
        "wtype",
        "use_fwd",
        "use_norm",
        "verbose",
        "bfac",
        "sigma",
        "psitype",
        "write_file",
        "res_factor",
        "ecc",
        "peri",
        "perturb",
        "interp",
        NULL
    };

    /*  Python objects needed throughout the computation.                     */
    PyObject *DLPInst;
    PyObject *tmp;
    PyObject *history;
    PyObject *rngreq;
    PyObject *perturb;

    /*  Set the default keyword options.                                      */

    /*  The kbmd20 is a new window, a modifed Kaiser-Bessel with alpha set to *
     *  two pi. The modification ensures the window goes to zero at its edges *
     *  while evaluating to one at the center, unlike the actual              *
     *  Kaiser-Bessel which is discontinuous at the edge of the window. The   *
     *  two pi factor is mostly guess work since it accurately reproduces the *
     *  PDS results. The real window used for that data is not known too me.  *
     *  The actual code for the window functions is in special_functions/     */
    self->wtype = "kbmd20";

    /*  Fresnel 4 is a new option, not mentioned in any of the papers but     *
     *  documented in our accompanying PDF. It uses Legendre polynomials to   *
     *  approximate the Fresnel kernel. It essentially takes Fresnels         *
     *  quadratic method to the next step, a quartic, hence the name. It is   *
     *  extremely fast (all of Rev007 takes less than a second) and very      *
     *  accurate for all but the most extreme occultations (like Rev133).     */
    self->psitype = "fresnel4";

    /*  Default range is "all", denoting [1.0, 400000.0]. We'll set later.    */
    rngreq = PyBytes_FromString("all");;

    /*  By default, forward computations are not run, FFTs are not used, and  *
     *  the run is silent (verbose is off).                                   */
    self->use_fwd = rssringoccs_False;
    self->verbose = rssringoccs_False;

    /*  Using the bfac guarantees accurate window sizes in the case of a poor *
     *  Allen deviation. Window normalization is also recommended since the   *
     *  integral is scaled by the width of the window, and hence for small    *
     *  window sizes the result might return close to zero.                   */
    self->bfac = rssringoccs_True;
    self->use_norm = rssringoccs_True;

    /*  The default sigma value is the one for Cassini.                       */
    self->sigma = 2.0e-13;

    /*  If res_factor was not set, set to 0.75. This value was specified by   *
     *  Essam Marouf as necessary to ensure the reconstruction matches the    *
     *  PDS results. No justification is known to me.                         */
    self->res_factor = 0.75;

    /*  The write files option allows the user to export the end result to a  *
     *  CSV. Default is set to false.                                         */
    self->write_file = rssringoccs_False;

    /*  The default geometry assumes the rings are circular, so we set both   *
     *  the eccentricity and the periapse to zero.                            */
    self->ecc = 0.0;
    self->peri = 0.0;

    /*  By default, polynomial interpolation is off.                          */
    self->interp = 0;

    /*  Default polynomial perturbation is off.                               */
    perturb = NULL;

    /*  Extract the inputs and keywords supplied by the user. If the data     *
     *  cannot be extracted, raise a type error and return to caller. A short *
     *  explaination of PyArg_ParseTupleAndKeywords. The inputs args and kwds *
     *  are somewhat straight-forward, they're the arguments and keywords     *
     *  passed by the string. The cryptic string is not straight-forward. The *
     *  | symbol means everything after need not be positional, and we can    *
     *  specify arguments and keywords by name when calling                   *
     *  DiffractionCorrection, for example                                    *
     *  DiffractionCorrect(..., wtype="blah"). O indicates a Python object,   *
     *  and d is a Python float. This is the DLP and res variables. The $     *
     *  symbold means everything after is optional. s is a string, p is a     *
     *  Boolean (p for "predicate"). b is an integer, and the colon : denotes *
     *  that the input list has ended.                                        */
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Od$OsppppdspdddOb:", kwlist,
                                     &DLPInst,          &self->input_res,
                                     &rngreq,           &self->wtype,
                                     &self->use_fwd,    &self->use_norm,
                                     &self->verbose,    &self->bfac,
                                     &self->sigma,      &self->psitype,
                                     &self->write_file, &self->res_factor,
                                     &self->ecc,        &self->peri,
                                     &perturb,          &self->interp))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rCould not parse input variables.\n\n"
            "\rInputs:\n"
            "\r\tDLPInst:  \tAn instance of the DLP Class.\n"
            "\r\tres:      \tRequested resolution in km (float).\n\n"
            "\rKeywords:\n"
            "\r\trng       \tThe requested range (str or list).\n"
            "\r\twtype     \tThe requested window type (str).\n"
            "\r\tuse_fwd   \tForward computation (bool).\n"
            "\r\tuse_norm  \tWindow normalization (bool).\n"
            "\r\tverbose   \tPrint status updates (bool).\n"
            "\r\tbfac      \tUse b-factor in window width (bool).\n"
            "\r\tsigma     \tThe Allen deviation (float).\n"
            "\r\tpsitype   \tRequested Frensel kernel approximation (str).\n"
            "\r\twrite_file\tWrite output to file (bool).\n"
            "\r\tres_factor\tScaling factor for resolution (float).\n"
            "\r\tecc       \tEccentricity of rings (bool).\n"
            "\r\tperi      \tPeriapse of rings (bool).\n"
            "\r\tperturb   \tRequested perturbation to Fresnel kernel (list).\n"
            "\r\tinterp    \tPolynomial interpolation to kernel (int).\n"
        );
        return -1;
    }

    if (self->verbose)
    {
        puts("Diffraction Correction:");
        puts("\tDiffraction Correction: Retrieving history from DLP...");
    }

    if (!PyObject_HasAttrString(DLPInst, "history"))
    {
        PyErr_Format(
            PyExc_AttributeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\thistory\n\n"
        );
        return -1;
    }
    else
        history = PyObject_GetAttrString(DLPInst, "history");

    /*  Store the dlp history inside of the DiffractionCorrection class. This *
     *  tmp, Py_INCREF, Py_XDECREF method is recommended in the Python C-API  *
     *  documentation as a means of safely storing the variable.              */
    tmp = self->dathist;
    Py_INCREF(history);
    self->dathist = history;
    Py_XDECREF(tmp);

    /*  If verbose was set, print a status update.                            */
    if (self->verbose)
        puts("\tDiffraction Correction: Converting Py DLP to C DLP...");

    dlp = rssringoccs_Py_DLP_to_C_DLP(DLPInst);

    if (dlp == NULL)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rFailed to pass variables to C. rssringoccs_Py_DLP_to_C_DLP\n"
            "\rreturned NULL. Returning.\n\n"
        );
        return -1;
    }

    if (dlp->error_occurred)
    {
        if (dlp->error_message == NULL)
        {
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rFailed to pass variables to C. rssringoccs_Py_DLP_to_C_DLP\n"
                "\rreturned a dlp with error_occurred set to True. No\n"
                "\rerror message was set. Returning.\n\n"
            );
        }
        else
        {
            PyErr_Format(PyExc_RuntimeError, "%s", dlp->error_message);
            free(dlp->error_message);
        }
        free(dlp);
        return -1;
    }

    /*  If verbose was set, print a status update.                            */
    if (self->verbose)
        puts("\tDiffraction Correction: Creating C Tau object...");

    tau = rssringoccs_Create_TAUObj(dlp, self->input_res * self->res_factor);

    if (self->verbose)
        puts("\tDiffraction Correction: Passing Py variables to tau...");

    rssringoccs_Get_Py_Vars_From_Self(tau, self);
    rssringoccs_Get_Py_Perturb(tau, perturb);
    rssringoccs_Get_Py_Range(tau, rngreq);
    rssringoccs_Tau_Set_WType(self->wtype, tau);
    rssringoccs_Tau_Set_Psitype(self->psitype, tau);

    if (self->verbose)
        puts("\tDiffraction Correction: Running reconstruction...");

    rssringoccs_Reconstruction(tau);

    if (self->verbose)
        puts("\tDiffraction Correction: Converting C tau to Py tau...");

    rssringoccs_C_Tau_to_Py_Tau(self, tau);

    if (tau == NULL)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrssringoccs_Create_TAUObj returned NULL for tau. Returning.\n\n"
        );

        return -1;
    }

    if (tau->error_occurred)
    {
        if (tau->error_message == NULL)
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rtau returned with error_occurred set to true but no\n"
                "\rerror message. Returning.\n\n"
            );
        else
            PyErr_Format(PyExc_RuntimeError, "%s\n", tau->error_message);

        rssringoccs_Destroy_Tau(&tau);
        return -1;
    }

    free(tau->psitype);
    free(tau->wtype);

    /*  We are now freeing the C tau object. The data pointers are still      *
     *  accessible via the self PyObject. Note, we are freeing the pointer to *
     *  the rssringoccs_TAUObj and NOT the pointers inside the object. The    *
     *  data is still available in self.                                      */
    free(tau);

    /*  Similarly, we free the DLP. This does not free the data from the      *
     *  input DLP PyObject. Those are also still available.                   */
    free(dlp);

    t2 = clock();

    if (self->verbose)
        printf(
            "\tDiffraction Correction: Computation Time %f\n",
            (double)(t2 - t1)/CLOCKS_PER_SEC
        );

    return 1;
}

static PyObject *
DiffractionCorrection(PyDiffrecObj *self, PyObject *Py_UNUSED(ignored))
{
    return PyUnicode_FromFormat("DiffractionCorrection");
}

static PyMethodDef DiffractionCorrection_methods[] =
{
    {
        "DiffractionCorrection",
        (PyCFunction) DiffractionCorrection,
        METH_NOARGS,
        "Diffraction correction class."
    },
    {NULL}
};

static PyMemberDef Custom_members[] = {
    {
        "rho_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, rho_km_vals),
        0,
        "Ring radius."
    },
    {
        "phase_rad_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, phase_rad_vals),
        0,
        "Raw diffracted phase."
    },
    {
        "B_rad_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, B_rad_vals),
        0,
        "Ring inclination angle."
    },
    {
        "D_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, D_km_vals),
        0,
        "Spacecraft to ring-intercept point distance."
    },
    {
        "rx_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, rx_km_vals),
        0,
        "x coordinate of the spacecraft in planetocentric frame."
    },
    {
        "ry_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, ry_km_vals),
        0,
        "y coordinate of the spacecraft in planetocentric frame."
    },
    {
        "rz_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, rz_km_vals),
        0,
        "z coordinate of the spacecraft in planetocentric frame."
    },
    {
        "F_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, F_km_vals),
        0,
        "Fresnel scale."
    },
    {
        "T_hat_fwd_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, T_hat_fwd_vals),
        0,
        "Complex transmittance of the forward modeling data."
    },
    {
        "T_hat_vals",
        T_OBJECT_EX, offsetof(PyDiffrecObj, T_hat_vals),
        0,
        "Raw diffraction data"
    },
    {
        "T_vals",
        T_OBJECT_EX, offsetof(PyDiffrecObj, T_vals),
        0,
        "Diffraction corrected data"
    },
    {
        "f_sky_hz_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, f_sky_hz_vals),
        0,
        "Frequency of the input signal"
    },
    {
        "p_norm_fwd_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, p_norm_fwd_vals),
        0,
        "Forward modeling of power"
    },
    {
        "p_norm_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, p_norm_vals),
        0,
        "Raw power data"
    },
    {
        "phase_fwd_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, phase_fwd_vals),
        0,
        "Forward modeling of phase"
    },
    {
        "phase_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, phase_vals),
        0,
        "Diffraction corrected phase"
    },
    {
        "phi_rad_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, phi_rad_vals),
        0,
        "Ring azimuth angle"
    },
    {
        "phi_rl_rad_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, phi_rl_rad_vals),
        0,
        "Ring longitude angle"
    },
    {
        "power_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, power_vals),
        0,
        "Diffraction corrected power"
    },
    {
        "raw_tau_threshold_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, raw_tau_threshold_vals),
        0,
        "Raw threshold optical depth"
    },
    {
        "rev_info",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, rev_info),
        0,
        "Information about the occultation"
    },
    {
        "rho_corr_pole_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, rho_corr_pole_km_vals),
        0,
        "Ring radius with pole correction."
    },
    {
        "rho_corr_timing_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, rho_corr_timing_km_vals),
        0,
        "Ring radius with timing correction."
    },
    {
        "rho_dot_kms_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, rho_dot_kms_vals),
        0,
        "Time derivative of the ring radius."
    },
    {
        "t_oet_spm_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, t_oet_spm_vals),
        0,
        "Observed event time in seconds past midnight"
    },
    {
        "t_ret_spm_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, t_ret_spm_vals),
        0,
        "Ring event time in seconds past midnight"
    },
    {
        "t_set_spm_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, t_set_spm_vals),
        0,
        "Spacecraft event time in seconds past midnight"
    },
    {
        "tau_threshold_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, tau_threshold_vals),
        0,
        "Diffraction corrected threshold optical depth"
    },
    {
        "tau_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, tau_vals),
        0,
        "Optical depth"
    },
    {
        "w_km_vals",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, w_km_vals),
        0,
        "window width as a function of ring radius"
    },
    {
        "dathist",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, dathist),
        0,
        "History of input dlp instance"
    },
    {
        "history",
        T_OBJECT_EX,
        offsetof(PyDiffrecObj, history),
        0,
        "History of this tau instance"
    },
    {
        "bfac",
        T_BOOL,
        offsetof(PyDiffrecObj, bfac),
        0,
        "Use of b-factor in window width"
    },
    {
        "verbose",
        T_BOOL,
        offsetof(PyDiffrecObj, verbose),
        0,
        "Print status updates"
    },
    {
        "use_norm",
        T_BOOL,
        offsetof(PyDiffrecObj, use_norm),
        0,
        "Use of window normalization"
    },
    {
        "use_fwd",
        T_BOOL,
        offsetof(PyDiffrecObj, use_fwd),
        0,
        "Forward modeling Boolean"
    },
    {
        NULL
    }  /* Sentinel */
};

static PyTypeObject DiffrecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "diffrec.DiffractionCorrection",
    .tp_doc = "Diffraction Correction class.",
    .tp_basicsize = sizeof(PyDiffrecObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) Diffrec_init,
    .tp_dealloc = (destructor) Diffrec_dealloc,
    .tp_members = Custom_members,
    .tp_methods = DiffractionCorrection_methods,
};

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "custom",
    .m_doc = "Module containing the rss_ringoccs class.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_diffrec(void)
{
    PyObject *m;
    int pymod_addobj;
    if (PyType_Ready(&DiffrecType) < 0)
        return NULL;

    m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    Py_INCREF(&DiffrecType);
    pymod_addobj = PyModule_AddObject(m, "DiffractionCorrection",
                                      (PyObject *) &DiffrecType);

    if (pymod_addobj < 0)
    {
        Py_DECREF(&DiffrecType);
        Py_DECREF(m);
        return NULL;
    }

    import_array();
    return m;
}
