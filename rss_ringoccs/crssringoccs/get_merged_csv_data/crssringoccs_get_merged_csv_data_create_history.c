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
crssringoccs_GetMergedCSVData_Create_History(crssringoccs_PyCSVObj *self,
                                             const char *dlpm_str)
{
    /*  Python objects needed throughout the computation.                     */
    PyObject *tmp, *input_variables, *input_keywords, *history;

    /*  Character array for the Python version. We'll use sprintf on this     *
     *  later with the macros provided in Python.h.                           */
    char python_version_string[16];

    /*  Python.h provided macros for the major and minor versioning of        *
     *  Python. To create a string out of this we use sprintf.                */
    sprintf(python_version_string, "%d.%d", PY_MAJOR_VERSION, PY_MINOR_VERSION);

    /*  Create a dictionary (Python object) with the input arguments.         */
    input_variables = Py_BuildValue(
        "{s:s}",
        "dlpm", dlpm_str
    );

    /*  There are no optional keywords, so set this to None.                  */
    Py_INCREF(Py_None);
    input_keywords = Py_None;

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
