/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
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
 ******************************************************************************/

/*  NULL and free are defined here.                                           */
#include <stdlib.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_strdup function declared here.                                       */
#include <libtmpl/include/tmpl_string.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  Macro for safely creating None objects.                                   */
#define MAKE_NONE(var)                                                         \
    do {                                                                       \
        PyObject *tmp = self->var;                                             \
        Py_INCREF(Py_None);                                                    \
        self->var = Py_None;                                                   \
        Py_XDECREF(tmp);                                                       \
    } while(0)

#define DESTROY_VAR(var) if (var) {free(var); var = NULL;}

/*  The init function for the dirrection correction class. This is the        *
 *  equivalent of the __init__ function defined in a normal python class.     */
int Diffrec_init(PyDiffrecObj *self, PyObject *args, PyObject *kwds)
{
    /*  Declare variables for a DLP and Tau object.                           */
    rssringoccs_DLPObj *dlp;
    rssringoccs_TAUObj *tau;

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
        "res_factor",
        "ecc",
        "peri",
        "perturb",
        NULL
    };

    /*  Python objects needed throughout the computation.                     */
    PyObject *DLPInst;
    PyObject *tmp;
    PyObject *dlp_tmp;

    /*  Set the default keyword options.                                      */

    /*  Default polynomial perturbation is off.                               */
    PyObject *perturb = NULL;

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
    PyObject *rngreq = PyUnicode_FromString("all");

    /*  By default, forward computations are not run, FFTs are not used, and  *
     *  the run is silent (verbose is off).                                   */
    self->use_fwd = tmpl_False;
    self->verbose = tmpl_False;

    /*  Using the bfac guarantees accurate window sizes in the case of a poor *
     *  Allen deviation. Window normalization is also recommended since the   *
     *  integral is scaled by the width of the window, and hence for small    *
     *  window sizes the result might return close to zero.                   */
    self->bfac = tmpl_True;
    self->use_norm = tmpl_True;

    /*  The default sigma value is the one for Cassini.                       */
    self->sigma = 2.0e-13;

    /*  If res_factor was not set, set to 0.75. This value was specified by   *
     *  Essam Marouf as necessary to ensure the reconstruction matches the    *
     *  PDS results. No justification is known to me.                         */
    self->res_factor = 0.75;

    /*  The default geometry assumes the rings are circular, so we set both   *
     *  the eccentricity and the periapse to zero.                            */
    self->ecc = 0.0;
    self->peri = 0.0;

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
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Od$OsppppdsdddO:", kwlist,
                                     &DLPInst,          &self->input_res,
                                     &rngreq,           &self->wtype,
                                     &self->use_fwd,    &self->use_norm,
                                     &self->verbose,    &self->bfac,
                                     &self->sigma,      &self->psitype,
                                     &self->res_factor, &self->ecc,
                                     &self->peri,       &perturb))
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
            "\r\tres_factor\tScaling factor for resolution (float).\n"
            "\r\tecc       \tEccentricity of rings (bool).\n"
            "\r\tperi      \tPeriapse of rings (bool).\n"
            "\r\tperturb   \tRequested perturbation to Fresnel kernel (list).\n"
        );
        return -1;
    }

    if (self->verbose)
    {
        puts("Diffraction Correction:");
        puts("\tDiffraction Correction: Retrieving history from DLP...");
    }

    /*  If verbose was set, print a status update.                            */
    if (self->verbose)
        puts("\tDiffraction Correction: Converting Py DLP to C DLP...");

    dlp = crssringoccs_Py_DLP_To_C_DLP(DLPInst);

    if (dlp == NULL)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rFailed to pass variables to C. rssringoccs_Py_DLP_To_C_DLP\n"
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
                "\rFailed to pass variables to C. rssringoccs_Py_DLP_To_C_DLP\n"
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

    tau = rssringoccs_Tau_Create_From_DLP(dlp, self->input_res * self->res_factor);

    if (self->verbose)
        puts("\tDiffraction Correction: Passing Py variables to tau...");

    crssringoccs_Get_Py_Vars_From_Tau_Self(tau, self);
    crssringoccs_Get_Py_Perturb(tau, perturb);
    crssringoccs_Get_Py_Range(tau, rngreq);

    rssringoccs_Tau_Set_Window_Type(self->wtype, tau);
    rssringoccs_Tau_Set_Psi_Type(self->psitype, tau);

    if (self->verbose)
        puts("\tDiffraction Correction: Running reconstruction...");

    rssringoccs_Reconstruction(tau);

    if (self->verbose)
        puts("\tDiffraction Correction: Converting C tau to Py tau...");

    crssringoccs_C_Tau_To_Py_Tau(self, tau);

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

        rssringoccs_Tau_Destroy(&tau);
        return -1;
    }

    /*  We are now freeing the C tau object. The data pointers are still      *
     *  accessible via the self PyObject. Note, we are freeing the pointer to *
     *  the rssringoccs_TAUObj and NOT the pointers inside the object. The    *
     *  data is still available in self.                                      */
    free(tau);

    /*  Similarly, we free the DLP. This does not free the data from the      *
     *  input DLP PyObject. Those are also still available.                   */
    free(dlp);

    if (self->verbose)
        puts("\tDiffraction Correction: Building arguments dictionary...");

    dlp_tmp = Py_BuildValue(
        "{s:O,s:d}",
        "dlp_inst", PyObject_GetAttrString(DLPInst, "history"),
        "res",      self->input_res
    );

    tmp = self->input_vars;
    Py_INCREF(dlp_tmp);
    self->input_vars = dlp_tmp;
    Py_XDECREF(tmp);

    if (self->verbose)
        puts("\tDiffraction Correction: Building keywords dictionary...");

    dlp_tmp = Py_BuildValue(
        "{s:O,s:s,s:s,s:d,s:d,s:d,s:d,s:O,s:O}",
        "rng",        rngreq,
        "wtype",      self->wtype,
        "psitype",    self->psitype,
        "sigma",      self->sigma,
        "ecc",        self->ecc,
        "peri",       self->peri,
        "res_factor", self->res_factor,
        "use_norm",   PyBool_FromLong(self->use_norm),
        "bfac",       PyBool_FromLong(self->bfac)
    );

    tmp = self->input_kwds;
    Py_INCREF(dlp_tmp);
    self->input_kwds = dlp_tmp;
    Py_XDECREF(tmp);

    self->outfiles = NULL;

    return 1;
}

