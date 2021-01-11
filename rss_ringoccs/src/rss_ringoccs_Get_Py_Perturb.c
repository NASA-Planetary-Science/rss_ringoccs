

void rssringoccs_Get_Py_Perturb(rssringoccs_TAUObj *tau, PyObject *perturb)
{
    PyObject *iter;
    PyObject *next;
    unsigned int n;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    /*  Check that the input perturb is a list with 5 elements.               */
    if (!perturb)
    {
        for (n=0; n<5; ++n)
            tau->perturb[n] = 0.0;
    }

    /*  If the user supplied a perturb list, parse it and extract values.     */
    else if (PyList_Check(perturb))
    {
        /*  If the list is not the correct size, raise an error.              */
        if (PyList_Size(perturb) != 5)
        {
            tau->error_occurred = rssringoccs_True;
            tau->error_message = rssringoccs_strdup(
                "\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Get_Py_Perturb\n\n"
                "\rInput perturb is a list but does not have 5 entries.\n"
                "\rperturb must be a list of five real numbers.\n"
            );
            return;
        }

        iter = PyObject_GetIter(perturb);

        /*  Loop over the elements of the list, see if they can be converted  *
         *  to doubles, and store them in the tau->perturb variable.          */
        for (n = 0; n < 5; ++n)
        {
            next = PyIter_Next(iter);

            /*  If the element is an integer, convert to double and save it.  */
            if (PyLong_Check(next))
                tau->perturb[n] = PyLong_AsDouble(next);

            /*  Convert from Python float to C double with PyFloat_AsDouble.  */
            else if (PyFloat_Check(next))
                tau->perturb[n] = PyFloat_AsDouble(next);

            /*  Invalid data type for one of the entries. Return with error.  */
            else
            {
                tau->error_occurred = rssringoccs_True;
                tau->error_message = rssringoccs_strdup(
                    "\rError Encountered: rss_ringoccs\n"
                    "\r\trssringoccs_Get_Py_Perturb\n\n"
                    "\rInput perturb has entries that are not real numbers.\n"
                    "\rAll five entries for the perturb list must be numbers.\n"
                );
                return;
            }
        }
    }

    /*  The input was not a list. Return with error.                          */
    else
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Get_Py_Perturb\n\n"
            "\rInput perturb is not a list.\n"
        );
        return;
    }
}
