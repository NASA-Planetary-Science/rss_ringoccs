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

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  The init function for the dirrection correction class. This is the        *
 *  equivalent of the __init__ function defined in a normal python class.     */
int
crssringoccs_DiffractionCorrection_Init(crssringoccs_PyDiffrecObj *self,
                                        PyObject *args,
                                        PyObject *kwds)
{
    rssringoccs_DLPObj *dlp;

    /*  Python objects needed throughout the computation.                     */
    PyObject *obj;

    /*  Declare variables for a DLP and Tau object.                           */
    double resolution;

    /*  The list of the keywords accepted by the DiffractionCorrection class. *
     *  dlp and res are REQUIRED inputs, the rest are optional. If the user   *
     *  does not provide these optional keywords, we must set them ourselves. */
    char * kwlist[] = {
        "dlp",
        "resolution",
        "rng",
        "wtype",
        "use_fwd",
        "use_norm",
        "verbose",
        "bfac",
        "sigma",
        "psitype",
        "resolution_factor",
        "eccentricity",
        "periapse",
        "perturb",
        NULL
    };

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
    const int success = PyArg_ParseTupleAndKeywords(
        args, kwds, "|Od$OsppppdsdddO:", kwlist,
        &obj,                     &self->input_resolution_km,
        &self->rngreq,            &self->wtype,
        &self->use_fwd,           &self->use_norm,
        &self->verbose,           &self->bfac,
        &self->sigma,             &self->psitype,
        &self->resolution_factor, &self->eccentricity,
        &self->periapse,          &self->perturb
    );

    if (!success)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tDiffractionCorrection\n\n"
            "\rCould not parse input variables.\n\n"
            "\rInputs:\n"
            "\r\tobj:               A Python object with all of the DLP data.\n"
            "\r\tresolution:        Requested resolution in km (float).\n\n"
            "\rKeywords:\n"
            "\r\trng:               Requested range (str or list).\n"
            "\r\twtype:             Requested window type (str).\n"
            "\r\tuse_fwd:           Forward computation (bool).\n"
            "\r\tuse_norm:          Window normalization (bool).\n"
            "\r\tverbose:           Print status updates (bool).\n"
            "\r\tbfac:              Use b-factor in window width (bool).\n"
            "\r\tsigma:             Allen deviation (float).\n"
            "\r\tpsitype:           Frensel kernel approxmiation (str).\n"
            "\r\tresolution_factor: Scale factor for resolution (float).\n"
            "\r\teccentricity:      Eccentricity of rings (float).\n"
            "\r\tperiapse:          Periapse of rings (float).\n"
            "\r\tperturb:           Fresnel kernel perturbation (list).\n"
        );

        return -1;
    }

    if (self->verbose)
        puts(
            "\rDiffraction Correction:\n"
            "\r\tDiffractionCorrection: Passing Python object to C..."
        );

    dlp = crssringoccs_PyObject_To_DLP(obj);

    if (!dlp)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tDiffractionCorrection\n\n"
            "\rFailed to pass variables to C. rssringoccs_Py_DLP_To_C_DLP\n"
            "\rreturned NULL. Returning.\n\n"
        );

        return -1;
    }

    if (dlp->error_occurred)
    {
        PyErr_Format(PyExc_RuntimeError, "%s", dlp->error_message);
        free(dlp);
        return -1;
    }

    /*  Lastly, copy the verbose Boolean.                                     */
    dlp->verbose = self->verbose;

    resolution = self->input_resolution_km * self->resolution_factor;
    self->tau = rssringoccs_Tau_Create_From_DLP(dlp, resolution);

    if (!self->tau)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tDiffractionCorrection\n\n"
            "\rrssringoccs_Tau_Create_From_DLP returned NULL.\n"
        );

        free(dlp);
        return -1;
    }

    crssringoccs_DiffractionCorrection_Set_Keywords(self);
    crssringoccs_DiffractionCorrection_Set_Perturb(self);
    crssringoccs_DiffractionCorrection_Set_Range(self);

    rssringoccs_Tau_Set_Window_Type(self->wtype, self->tau);
    rssringoccs_Tau_Set_Psi_Type(self->psitype, self->tau);
    rssringoccs_Reconstruction(self->tau);

    crssringoccs_DiffractionCorrection_Create_Argument_Dictionary(self, obj);
    crssringoccs_DiffractionCorrection_Create_Keyword_Dictionary(self);
    crssringoccs_DiffractionCorrection_Finish(self, obj);

    if (self->tau->error_occurred)
    {
        PyErr_Format(PyExc_RuntimeError, "%s\n", self->tau->error_message);
        rssringoccs_Tau_Destroy(&self->tau);
        free(dlp);
        return -1;
    }

    return 1;
}
