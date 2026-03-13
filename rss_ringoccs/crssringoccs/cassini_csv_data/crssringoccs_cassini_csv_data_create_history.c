#include "../crssringoccs.h"

/*  Functions for creating the history (user name, OS, date, etc.).           */
#include <libtmpl/include/tmpl_calendar.h>
#include <libtmpl/include/tmpl_utility.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>


#define RSSRINGOCCS_TO_STRING(x) #x
#define RSSRINGOCCS_MAKE_STRING(x) RSSRINGOCCS_TO_STRING(x)

static const char * const crssringoccs_python_version =
    RSSRINGOCCS_MAKE_STRING(PY_MAJOR_VERSION)
    "."
    RSSRINGOCCS_MAKE_STRING(PY_MINOR_VERSION);

#ifndef __STDC_VERSION__
    const char * const crssringoccs_c_version = "Unknown";
#else
    const char * const crssringoccs_c_version =
        RSSRINGOCCS_MAKE_STRING(__STDC_VERSION__);
#endif

void
crssringoccs_CassiniCSVData_Create_History(crssringoccs_PyCSVObj *self,
                                           const char *geo_str,
                                           const char *cal_str,
                                           const char *dlp_str,
                                           const char *tau_str,
                                           tmpl_Bool use_deprecated)
{
    /*  Python objects needed throughout the computation.                     */
    PyObject *input_variables = NULL;
    PyObject *input_keywords = NULL;
    PyObject *history = NULL;
    PyObject *py_use_deprecated = NULL;

    /*  If the tau variable was not set, set the string to "None" to prevent  *
     *  segfaults from trying to access a NULL pointer.                       */
    if (!tau_str)
        tau_str = "None";

    /*  Python booleans are a type of PyObject. They can be created by        *
     *  casting our Boolean to a long int as follows.                         */
    py_use_deprecated = PyBool_FromLong(use_deprecated);

    /*  Check for errors, this object should no longer be NULL.               */
    if (!py_use_deprecated)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_CassiniCSVData_Create_History\n\n"
            "\rPyBool_FromLong returned NULL.\n\n"
        );

        return;
    }

    /*  Create a dictionary (Python object) with the input arguments.         */
    input_variables = Py_BuildValue(
        "{s:s,s:s,s:s}",
        "geo", geo_str,
        "cal", cal_str,
        "dlp", dlp_str
    );

    /*  Check that the dictionary was successfully created.                   */
    if (!input_variables)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_CassiniCSVData_Create_History\n\n"
            "\rPy_BuildValue returned NULL for input_variables.\n\n"
        );

        goto CLEANUP;
    }

    /*  Create a dictionary with the optional keywords.                       */
    input_keywords = Py_BuildValue(
        "{s:s,s:O}",
        "tau", tau_str,
        "use_deprecated", py_use_deprecated
    );

    /*  Check that the dictionary was successfully created.                   */
    if (!input_keywords)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_CassiniCSVData_Create_History\n\n"
            "\rPy_BuildValue returned NULL for input_keywords.\n\n"
        );

        goto CLEANUP;
    }

    /*  Create the history object, which is a Python dictionary.              */
    history = Py_BuildValue(
        "{s:s,s:s,s:s,s:s,s:s,s:s,s:s,s:s,s:O,s:O}",
        "rss_ringoccs Version", rssringoccs_Version(),
        "libtmpl Version", tmpl_Version(),
        "Python Version", crssringoccs_python_version,
        "C Version", crssringoccs_c_version,
        "Host Name", tmpl_Host_Name(),
        "User Name", tmpl_User_Name(),
        "Run Date", tmpl_Local_Calendar_Date_And_Time(),
        "Operating System", tmpl_Operating_System(),
        "Positional Args", input_variables,
        "Keyword Args", input_keywords
    );

    /*  Check that the dictionary was successfully created.                   */
    if (!history)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcrssringoccs_CassiniCSVData_Create_History\n\n"
            "\rPy_BuildValue returned NULL for history.\n\n"
        );

        goto CLEANUP;
    }

    /*  Begin reference counting for the new history object.                  */
    Py_XSETREF(self->history, history);

    CLEANUP:
        Py_CLEAR(py_use_deprecated);
        Py_CLEAR(input_variables);
        Py_CLEAR(input_keywords);
}
/*  End of crssringoccs_CassiniCSVData_Create_History.                        */

#undef RSSRINGOCCS_TO_STRING
#undef RSSRINGOCCS_MAKE_STRING
