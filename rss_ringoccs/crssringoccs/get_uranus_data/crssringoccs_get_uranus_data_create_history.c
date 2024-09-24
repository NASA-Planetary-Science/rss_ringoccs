#include "../crssringoccs.h"

/*  sprintf routine provided here.                                            */
#include <stdio.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  Functions for creating the history (user name, OS, date, etc.).           */
#include <libtmpl/include/tmpl_calendar.h>
#include <libtmpl/include/tmpl_utility.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>

void
crssringoccs_GetUranusData_Create_History(crssringoccs_PyCSVObj *self,
                                          const char *geo_str,
                                          const char *dlp_str,
                                          const char *tau_str,
                                          tmpl_Bool dlp_in_radians)
{
    /*  Python objects needed throughout the computation.                     */
    PyObject *tmp, *input_variables, *input_keywords, *history, *py_bool;

    /*  Character array for the Python version. We'll use sprintf on this     *
     *  later with the macros provided in Python.h.                           */
    char python_version_string[16];

    /*  If the tau variable was not set, set the string to "None" to prevent  *
     *  segfaults from trying to access a NULL pointer.                       */
    if (!tau_str)
        tau_str = "None";

    /*  Python.h provided macros for the major and minor versioning of        *
     *  Python. To create a string out of this we use sprintf.                */
    sprintf(python_version_string, "%d.%d", PY_MAJOR_VERSION, PY_MINOR_VERSION);

    /*  Python booleans are a type of PyObject. They can be created by        *
     *  casting our Boolean to a long int as follows.                         */
    py_bool = PyBool_FromLong((long int)dlp_in_radians);

    /*  Create a dictionary (Python object) with the input arguments.         */
    input_variables = Py_BuildValue(
        "{s:s,s:s}",
        "geo", geo_str,
        "dlp", dlp_str
    );

    /*  Create a dictionary with the optional keywords. The Boolean is a      *
     *  PyOBject. "N" means we do not create a new reference to it, as        *
     *  opposed to "O" which does. This means when we destroy the history     *
     *  object, the py_bool object is deleted as well.                        */
    input_keywords = Py_BuildValue(
        "{s:s,s:N}",
        "tau", tau_str,
        "dlp_in_radians", py_bool
    );

    /*  Create the history object, which is a Python dictionary.              */
    history = Py_BuildValue(
        "{s:s,s:s,s:s,s:s,s:s,s:s,s:s,s:N,s:N}",
        "rss_ringoccs Version", rssringoccs_Get_Version(),
        "libtmpl Version", tmpl_Version(),
        "Python Version", python_version_string,
        "Host Name", tmpl_Host_Name(),
        "User Name", tmpl_User_Name(),
        "Run Date", tmpl_Local_Calendar_Date_And_Time(),
        "Operating System", tmpl_Operating_System(),
        "Positional Args", input_variables,
        "Keyword Args", input_keywords
    );

    /*  Begin reference counting for the new history object.                  */
    tmp = self->history;
    Py_INCREF(history);
    self->history = history;
    Py_XDECREF(tmp);
}
/*  End of crssringoccs_GetUranusData_Create_History.                         */
