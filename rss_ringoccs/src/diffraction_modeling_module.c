


#define VarToString(Var) (#Var)

#define _array_from_three(x1, x2, x3, y, dim, intype, outtype, f) ({           \
    long i;                                                                    \
    outtype *out = (outtype *)malloc(sizeof(outtype)*dim);                     \
    for (i=0; i<dim; ++i)                                                      \
        out[i] = (*f)(((intype *)x1)[i], x2, x3);                              \
    y = out;                                                                   \
})

#define _array_from_four(x1, x2, x3, x4, y, dim, intype, outtype, f) ({        \
    long i;                                                                    \
    outtype *out = (outtype *)malloc(sizeof(outtype)*dim);                     \
    for (i=0; i<dim; ++i)                                                      \
        out[i] = (*f)(((intype *)x1)[i], x2, x3, x4);                          \
    y = out;                                                                   \
})

/*  The Python wrappers for several of the functions is completely identical  *
 *  with the exception of which C function needs to be called, and what the   *
 *  function will be called at the Python level. We create a macro with two   *
 *  inputs:                                                                   *
 *      FuncName:                                                             *
 *          The name of the function at the python level.                     *
 *      CName:                                                                *
 *          The name of the function at the C level with the data type        *
 *          excluded. For example, if you wanted to wrap the Sinc_Double      *
 *          function defined in special_functions/, CName whould be Sinc.     *
 *  NOTE:                                                                     *
 *      Functions built using this macro MUST require four variables of the   *
 *      following types:                                                      *
 *          rho:                                                              *
 *              A numpy array or list of real numbers, a float, or an int.    *
 *          a:                                                                *
 *              A positive real number.                                       *
 *          b:                                                                *
 *              A positive real number GREATER than a.                        *
 *          F:                                                                *
 *              A positive real number.                                       *
 *      Ideally, only use this for the diffraction modeling functions which   *
 *      require four inputs. rho should be ring radius, a and b should be     *
 *      specific radii, and F is the fresnel scale. rho, a, b, and F should   *
 *      all have the same units, preferably kilometers (km).                  */
#define RSS_RINGOCCSDiffModelingFourVars(FuncName, CName)                      \
static PyObject *FuncName(PyObject *self, PyObject *args)                      \
{                                                                              \
    /*  We'll need output and capsule for safely creating the output array   */\
    /*  and ensuring we don't have a memory leak. rho is the input variable. */\
    PyObject *output, *capsule,  *rho;                                         \
                                                                               \
    /*  a, b, and F are the input real numbers. Currently the code only      */\
    /*  accepts floats and ints, and does not deal with fancy Python objects.*/\
    double a, b, F;                                                            \
                                                                               \
    /*  Variable for the size of the input array which we'll retrieve later. */\
    long dim;                                                                  \
                                                                               \
    /*  Try to extract the values passed from the user, return error if fail.*/\
    if (!PyArg_ParseTuple(args, "Oddd", &rho, &a, &b, &F))                     \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For all functions using this macro, a <= b. Check this.              */\
    if (a >= b)                                                                \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_ValueError,                                                  \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not less than outer radius (i.e. a >= b).\n",   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  We're modeling in polar coordinates, so radius values should not be  */\
    /*  negative. Check this and raise an error if a < 0.                    */\
    else if (a < 0.0)                                                          \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not positive. (i.e. a<0)\n",                    \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  The Fresnel scale should never be zero or negative. Check.           */\
    else if (F <= 0.0)                                                         \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rFresnel scale is not positive (i.e. F<=0).\n",                  \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If the user supplied a float or int, simply pass the inputs to the   */\
    /*  C function and return the value back to Python.                      */\
    if (PyLong_Check(rho) || PyFloat_Check(rho))                               \
    {                                                                          \
        /*  Convert a python float or int into a C double.                   */\
        double val = PyFloat_AsDouble(rho);                                    \
                                                                               \
        /*  The function returns a complex variable, so declare one.         */\
        complex double out;                                                    \
                                                                               \
        /*  Now get the computed value and convert the C complex double to a */\
        /*  Python complex value and return. The ## symbol tells the         */\
        /*  preprocessor to concatenate the words for the compiler. So if we */\
        /*  pass Gap_Diffraction as CName, the compiler will see             */\
        /*      out = Gap_Diffraction_Double(val, a, b, F);                  */\
        /*  which is exactly what we want.                                   */\
        out = CName##_Double(val, a, b, F);                                    \
                                                                               \
        /*  To build a Python complex value we need to pass the real and     */\
        /*  imaginary parts of out into PyComplex_FromDoubles. We can easily */\
        /*  access these values using creal and cimag, which are functions   */\
        /*  found in complex.h.                                              */\
        return PyComplex_FromDoubles(creal(out), cimag(out));                  \
    }                                                                          \
                                                                               \
    /*  As of v1.3 we allow lists to be passed to these routines.            */\
    else if (PyList_Check(rho))                                                \
    {                                                                          \
        /*  If the user passed a list, we'll need to return one. Create a    */\
        /*  variable for indexing over the input list.                       */\
        long i;                                                                \
                                                                               \
        /*  Create another variable for indexing over the list. This will be */\
        /*  the object corresponding to the value of the index i.            */\
        PyObject *ith_item;                                                    \
                                                                               \
        /*  And lastly, a double for storing the input value and a complex   */\
        /*  double for storing the output value.                             */\
        double val;                                                            \
        double complex out;                                                    \
                                                                               \
        /*  Get the number of elements in the list.                          */\
        dim = PyList_Size(rho);                                                \
                                                                               \
        /*  If the input list is empty, return with an error.                */\
        if (dim == 0)                                                          \
        {                                                                      \
            PyErr_Format(                                                      \
                PyExc_TypeError,                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\tdiffrec.%s\n\n"                                           \
                "\rInput list is empty.\n", VarToString(FuncName)              \
            );                                                                 \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  Create a new list to be returned to the user.                    */\
        output = PyList_New(dim);                                              \
                                                                               \
        /*  Loop over the entries of the list, examine values for error, and */\
        /*  try to perform computations if legal inputs are given.           */\
        for (i=0; i<dim; ++i)                                                  \
        {                                                                      \
            /*  Extract the python object in the ith slot of the list.       */\
            ith_item = PyList_GET_ITEM(rho, i);                                \
                                                                               \
            /*  The list should be homogeneous with real numbers only. Check */\
            /*  each item and raise an error otherwise.                      */\
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tdiffrec.%s\n\n"                              \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            /*  Convert the python float to a C double.                      */\
            val = PyFloat_AsDouble(ith_item);                                  \
                                                                               \
            /*  Compute the value, set the item in the list, and move on to  */\
            /*  to the next one. As explained before, ## concatenates words  */\
            /*  so CName##_Double compiles as CName_Double.                  */\
            out = CName##_Double(val, a, b, F);                                \
            PyList_SET_ITEM(output, i,                                         \
                            PyComplex_FromDoubles(creal(out), cimag(out)));    \
        }                                                                      \
                                                                               \
        /*  The list created without errors, so just return it.              */\
        return output;                                                         \
    }                                                                          \
                                                                               \
    /*  If the input was not a list or a number, it must be a numpy array.   */\
    /*  Check this and return error otherwise.                               */\
    else if (!PyArray_Check(rho))                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check the input array to make sure it is valid. It is important to   */\
    /*  check the number of dimensions BEFORE setting the dim variable. If   */\
    /*  the input array is empty, PyArray_DIMS(rho) will be a NULL pointer.  */\
    /*  Trying to access it with dim = PyArray_DIMS(rho)[0] will create a    */\
    /*  segmentation fault, crashing the Python interpreter.                 */\
    if (PyArray_NDIM((PyArrayObject *)rho) != 1)                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is not one-dimensional.\n",                   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Useful information about the data. The numpy API functions want a    */\
    /*  PyArrayObject pointer, so we need to cast rho as this type.          */\
    int typenum = PyArray_TYPE((PyArrayObject *)rho);                          \
    void *data  = PyArray_DATA((PyArrayObject *)rho);                          \
    dim         = PyArray_DIMS((PyArrayObject *)rho)[0];                       \
                                                                               \
    /*  This variable is for the data type of the output variable. It should */\
    /*  the complex valued equivalent of the input. So NPY_FLOAT should      */\
    /*  return NPY_CLOAT. All integer types will return NPY_CDOUBLE.         */\
    /*  If the input array was an int-like array (long, short, char, etc)    */\
    /*  then the output will be complex double by default/                   */\
    int outtype;                                                               \
                                                                               \
    /*  This void pointer will point to the output data we'll create later.  */\
    void *y;                                                                   \
                                                                               \
    /*  If the input array is empty, return with an error.                   */\
    if (dim == 0)                                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is empty.\n", VarToString(FuncName)           \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For float, double, and long double precision (which corresponds to   */\
    /*  numpy.float32, numpy.float64, and numpy.float128, respectively) we   */\
    /*  have functions which compute at these respective precisions. Pass    */\
    /*  the values into the _array_from_four macro which will loop over the  */\
    /*  data and create a new output pointer, which y will then point to.    */\
    if (typenum == NPY_FLOAT)                                                  \
    {                                                                          \
        _array_from_four(data, a, b, F, y, dim, float,                         \
                         complex float, CName##_Float);                        \
        outtype = NPY_CFLOAT;                                                  \
    }                                                                          \
    else if (typenum == NPY_DOUBLE)                                            \
    {                                                                          \
        _array_from_four(data, a, b, F, y, dim, double,                        \
                         complex double, CName##_Double);                      \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
    else if (typenum == NPY_LONGDOUBLE)                                        \
    {                                                                          \
        _array_from_four(data, a, b, F, y, dim, long double,                   \
                         complex long double, CName##_Long_Double);            \
        outtype = NPY_CLONGDOUBLE;                                             \
    }                                                                          \
                                                                               \
    /*  For all other types we try to convert to double and compute.         */\
    else                                                                       \
    {                                                                          \
                                                                               \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);          \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(newrho))                                                         \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tdiffrec.%s\n\n"                                   \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
        else                                                                   \
            _array_from_four(data, a, b, F, y, dim, double,                    \
                             complex double, CName##_Double);                  \
                                                                               \
        /*  Retrieve a pointer to the data inside of newrho and set outtype. */\
        data = PyArray_DATA((PyArrayObject *)newrho);                          \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
                                                                               \
    /*  Create a numpy array from the newly computed data and set a capsule  */\
    /*  for it. This ensures when the corresponding Python variable is       */\
    /*  deleted or removed the memory allocated to the respective C pointer  */\
    /*  is freed. Skipping this will result in memory leaks!                 */\
    output  = PyArray_SimpleNewFromData(1, &dim, outtype, y);                  \
    capsule = PyCapsule_New(y, NULL, capsule_cleanup);                         \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

/*  Same as RSS_RINGOCCSDiffModelingFourVars but for when three variables are *
 *  needed instead of four.                                                   *
 *  inputs:                                                                   *
 *      FuncName:                                                             *
 *          The name of the function at the python level.                     *
 *      CName:                                                                *
 *          The name of the function at the C level with the data type        *
 *          excluded. For example, if you wanted to wrap the Sinc_Double      *
 *          function defined in special_functions/, CName whould be Sinc.     *
 *  NOTE:                                                                     *
 *      Functions built using this macro MUST require three variables of the  *
 *      following types:                                                      *
 *          rho:                                                              *
 *              A numpy array or list of real numbers, a float, or an int.    *
 *          a:                                                                *
 *              A positive real number.                                       *
 *          F:                                                                *
 *              A positive real number.                                       *
 *      Ideally, only use this for the diffraction modeling functions which   *
 *      require three inputs. rho should be ring radius, a should be a        *
 *      specific radii, and F is the fresnel scale. rho, a, and F should all  *
 *      have the same units, preferably kilometers (km).                      */
#define RSS_RINGOCCSDiffModelingThreeVars(FuncName, CName)                     \
static PyObject *FuncName(PyObject *self, PyObject *args)                      \
{                                                                              \
    /*  We'll need output and capsule for safely creating the output array   */\
    /*  and ensuring we don't have a memory leak. rho is the input variable. */\
    PyObject *output, *capsule,  *rho;                                         \
                                                                               \
    /*  a and F are the input real numbers. Currently the code only          */\
    /*  accepts floats and ints, and does not deal with fancy Python objects.*/\
    double a, F;                                                               \
                                                                               \
    /*  Variable for the size of the input array which we'll retrieve later. */\
    long dim;                                                                  \
                                                                               \
    /*  Try to extract the values passed from the user, return error if fail.*/\
    if (!PyArg_ParseTuple(args, "Odd", &rho, &a, &F))                          \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  We're modeling in polar coordinates, so radius values should not be  */\
    /*  negative. Check this and raise an error if a < 0.                    */\
    if (a < 0.0)                                                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not positive. (i.e. a<0)\n",                    \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  The Fresnel scale should never be zero or negative. Check.           */\
    else if (F <= 0.0)                                                         \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rFresnel scale is not positive (i.e. F<=0).\n",                  \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If the user supplied a float or int, simply pass the inputs to the   */\
    /*  C function and return the value back to Python.                      */\
    if (PyLong_Check(rho) || PyFloat_Check(rho))                               \
    {                                                                          \
        /*  Convert a python float or int into a C double.                   */\
        double val = PyFloat_AsDouble(rho);                                    \
                                                                               \
        /*  The function returns a complex variable, so declare one.         */\
        complex double out;                                                    \
                                                                               \
        /*  Now get the computed value and convert the C complex double to a */\
        /*  Python complex value and return. The ## symbol tells the         */\
        /*  preprocessor to concatenate the words for the compiler. So if we */\
        /*  pass Gap_Diffraction as CName, the compiler will see             */\
        /*      out = Gap_Diffraction_Double(val, a, b, F);                  */\
        /*  which is exactly what we want.                                   */\
        out = CName##_Double(val, a, F);                                       \
                                                                               \
        /*  To build a Python complex value we need to pass the real and     */\
        /*  imaginary parts of out into PyComplex_FromDoubles. We can easily */\
        /*  access these values using creal and cimag, which are functions   */\
        /*  found in complex.h.                                              */\
        return PyComplex_FromDoubles(creal(out), cimag(out));                  \
    }                                                                          \
                                                                               \
    /*  As of v1.3 we allow lists to be passed to these routines.            */\
    else if (PyList_Check(rho))                                                \
    {                                                                          \
        /*  If the user passed a list, we'll need to return one. Create a    */\
        /*  variable for indexing over the input list.                       */\
        long i;                                                                \
                                                                               \
        /*  Create another variable for indexing over the list. This will be */\
        /*  the object corresponding to the value of the index i.            */\
        PyObject *ith_item;                                                    \
                                                                               \
        /*  And lastly, a double for storing the input value and a complex   */\
        /*  double for storing the output value.                             */\
        double val;                                                            \
        double complex out;                                                    \
                                                                               \
        /*  Get the number of elements in the list.                          */\
        dim = PyList_Size(rho);                                                \
                                                                               \
        /*  If the input list is empty, return with an error.                */\
        if (dim == 0)                                                          \
        {                                                                      \
            PyErr_Format(                                                      \
                PyExc_TypeError,                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\tdiffrec.ringlet_diffraction\n\n"                          \
                "\rInput list is empty.\n", VarToString(FuncName)              \
            );                                                                 \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  Create a new list to be returned to the user.                    */\
        output = PyList_New(dim);                                              \
                                                                               \
        /*  Loop over the entries of the list, examine values for error, and */\
        /*  try to perform computations if legal inputs are given.           */\
        for (i=0; i<dim; ++i)                                                  \
        {                                                                      \
            /*  Extract the python object in the ith slot of the list.       */\
            ith_item = PyList_GET_ITEM(rho, i);                                \
                                                                               \
            /*  The list should be homogeneous with real numbers only. Check */\
            /*  each item and raise an error otherwise.                      */\
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tdiffrec.%s\n\n"                              \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            /*  Convert the python float to a C double.                      */\
            val = PyFloat_AsDouble(ith_item);                                  \
                                                                               \
            /*  Compute the value, set the item in the list, and move on to  */\
            /*  to the next one. As explained before, ## concatenates words  */\
            /*  so CName##_Double compiles as CName_Double.                  */\
            out = CName##_Double(val, a, F);                                   \
            PyList_SET_ITEM(output, i,                                         \
                            PyComplex_FromDoubles(creal(out), cimag(out)));    \
        }                                                                      \
                                                                               \
        /*  The list created without errors, so just return it.              */\
        return output;                                                         \
    }                                                                          \
                                                                               \
    /*  If the input was not a list or a number, it must be a numpy array.   */\
    /*  Check this and return error otherwise.                               */\
    else if (!PyArray_Check(rho))                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check the input array to make sure it is valid. It is important to   */\
    /*  check the number of dimensions BEFORE setting the dim variable. If   */\
    /*  the input array is empty, PyArray_DIMS(rho) will be a NULL pointer.  */\
    /*  Trying to access it with dim = PyArray_DIMS(rho)[0] will create a    */\
    /*  segmentation fault, crashing the Python interpreter.                 */\
    if (PyArray_NDIM((PyArrayObject *)rho) != 1)                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is not one-dimensional.\n",                   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Useful information about the data. The numpy API functions want a    */\
    /*  PyArrayObject pointer, so we need to cast rho as this type.          */\
    int typenum = PyArray_TYPE((PyArrayObject *)rho);                          \
    void *data  = PyArray_DATA((PyArrayObject *)rho);                          \
    dim         = PyArray_DIMS((PyArrayObject *)rho)[0];                       \
                                                                               \
    /*  This variable is for the data type of the output variable. It should */\
    /*  the complex valued equivalent of the input. So NPY_FLOAT should      */\
    /*  return NPY_CLOAT. All integer types will return NPY_CDOUBLE.         */\
    /*  If the input array was an int-like array (long, short, char, etc)    */\
    /*  then the output will be complex double by default/                   */\
    int outtype;                                                               \
                                                                               \
    /*  This void pointer will point to the output data we'll create later.  */\
    void *y;                                                                   \
                                                                               \
    /*  If the input array is empty, return with an error.                   */\
    if (dim == 0)                                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is empty.\n", VarToString(FuncName)           \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For float, double, and long double precision (which corresponds to   */\
    /*  numpy.float32, numpy.float64, and numpy.float128, respectively) we   */\
    /*  have functions which compute at these respective precisions. Pass    */\
    /*  the values into the _array_from_three macro which will loop over the */\
    /*  data and create a new output pointer, which y will then point to.    */\
    if (typenum == NPY_FLOAT)                                                  \
    {                                                                          \
        _array_from_three(data, a, F, y, dim, float,                           \
                          complex float, CName##_Float);                       \
        outtype = NPY_CFLOAT;                                                  \
    }                                                                          \
    else if (typenum == NPY_DOUBLE)                                            \
    {                                                                          \
        _array_from_three(data, a, F, y, dim, double,                          \
                          complex double, CName##_Double);                     \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
    else if (typenum == NPY_LONGDOUBLE)                                        \
    {                                                                          \
        _array_from_three(data, a, F, y, dim, long double,                     \
                         long double, CName##_Long_Double);                    \
        outtype = NPY_CLONGDOUBLE;                                             \
    }                                                                          \
                                                                               \
    /*  For all other types we try to convert to double and compute.         */\
    else                                                                       \
    {                                                                          \
                                                                               \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);          \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(newrho))                                                         \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tdiffrec.%s\n\n"                                   \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
        else                                                                   \
            _array_from_three(data, a, F, y, dim, double,                      \
                              complex double, CName##_Double);                 \
                                                                               \
        /*  Retrieve a pointer to the data inside of newrho and set outtype. */\
        data = PyArray_DATA((PyArrayObject *)newrho);                          \
        outtype = NPY_CDOUBLE;                                                 \
    }                                                                          \
                                                                               \
    /*  Create a numpy array from the newly computed data and set a capsule  */\
    /*  for it. This ensures when the corresponding Python variable is       */\
    /*  deleted or removed the memory allocated to the respective C pointer  */\
    /*  is freed. Skipping this will result in memory leaks!                 */\
    output  = PyArray_SimpleNewFromData(1, &dim, outtype, y);                  \
    capsule = PyCapsule_New(y, NULL, capsule_cleanup);                         \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

/*  Same as RSS_RINGOCCSDiffModelingFourVars but when the output is real and  *
 *  not complex. It has the following inputs:                                 *
 *      FuncName:                                                             *
 *          The name of the function at the python level.                     *
 *      CName:                                                                *
 *          The name of the function at the C level with the data type        *
 *          excluded. For example, if you wanted to wrap the Sinc_Double      *
 *          function defined in special_functions/, CName whould be Sinc.     *
 *  NOTE:                                                                     *
 *      Functions built using this macro MUST require four variables of the   *
 *      following types:                                                      *
 *          rho:                                                              *
 *              A numpy array or list of real numbers, a float, or an int.    *
 *          a:                                                                *
 *              A positive real number.                                       *
 *          b:                                                                *
 *              A positive real number GREATER than a.                        *
 *          F:                                                                *
 *              A positive real number.                                       *
 *      Ideally, only use this for the diffraction modeling functions which   *
 *      require four inputs and return reals, such as the phase computing     *
 *      routines. rho should be ring radius, a and b should be specific radii,*
 *      and F is the fresnel scale. rho, a, b, and F should all have the same *
 *      units, preferably kilometers (km). Output is usually radians.         */
#define RSS_RINGOCCSDiffModelingFourVarsREAL(FuncName, CName)                  \
static PyObject *FuncName(PyObject *self, PyObject *args)                      \
{                                                                              \
    /*  We'll need output and capsule for safely creating the output array   */\
    /*  and ensuring we don't have a memory leak. rho is the input variable. */\
    PyObject *output, *capsule,  *rho;                                         \
                                                                               \
    /*  a, b, and F are the input real numbers. Currently the code only      */\
    /*  accepts floats and ints, and does not deal with fancy Python objects.*/\
    double a, b, F;                                                            \
                                                                               \
    /*  Variable for the size of the input array which we'll retrieve later. */\
    long dim;                                                                  \
                                                                               \
    /*  Try to extract the values passed from the user, return error if fail.*/\
    if (!PyArg_ParseTuple(args, "Oddd", &rho, &a, &b, &F))                     \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For all functions using this macro, a <= b. Check this.              */\
    if (a >= b)                                                                \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_ValueError,                                                  \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not less than outer radius (i.e. a >= b).\n",   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  We're modeling in polar coordinates, so radius values should not be  */\
    /*  negative. Check this and raise an error if a < 0.                    */\
    else if (a < 0.0)                                                          \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInner radius is not positive. (i.e. a<0)\n",                    \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  The Fresnel scale should never be zero or negative. Check.           */\
    else if (F <= 0.0)                                                         \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rFresnel scale is not positive (i.e. F<=0).\n",                  \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  If the user supplied a float or int, simply pass the inputs to the   */\
    /*  C function and return the value back to Python.                      */\
    if (PyLong_Check(rho) || PyFloat_Check(rho))                               \
    {                                                                          \
        /*  Convert a python float or int into a C double.                   */\
        double val = PyFloat_AsDouble(rho);                                    \
                                                                               \
        /*  Now get the computed value and convert the C complex double to a */\
        /*  Python complex value and return. The ## symbol tells the         */\
        /*  preprocessor to concatenate the words for the compiler. So if we */\
        /*  pass Gap_Diffraction as CName, the compiler will see             */\
        /*      out = Gap_Diffraction_Double(val, a, b, F);                  */\
        /*  which is exactly what we want.                                   */\
        val = CName##_Double(val, a, b, F);                                    \
                                                                               \
        /*  Now convert the C double back into a Python float and return.    */\
        return PyFloat_FromDouble(val);                                        \
    }                                                                          \
                                                                               \
    /*  As of v1.3 we allow lists to be passed to these routines.            */\
    else if (PyList_Check(rho))                                                \
    {                                                                          \
        /*  If the user passed a list, we'll need to return one. Create a    */\
        /*  variable for indexing over the input list.                       */\
        long i;                                                                \
                                                                               \
        /*  Create another variable for indexing over the list. This will be */\
        /*  the object corresponding to the value of the index i.            */\
        PyObject *ith_item;                                                    \
                                                                               \
        /*  And lastly, a double for storing the input value.                */\
        double val;                                                            \
                                                                               \
        /*  Get the number of elements in the list.                          */\
        dim = PyList_Size(rho);                                                \
                                                                               \
        /*  If the input list is empty, return with an error.                */\
        if (dim == 0)                                                          \
        {                                                                      \
            PyErr_Format(                                                      \
                PyExc_TypeError,                                               \
                "\n\rError Encountered: rss_ringoccs\n"                        \
                "\r\tdiffrec.%s\n\n"                                           \
                "\rInput list is empty.\n", VarToString(FuncName)              \
            );                                                                 \
            return NULL;                                                       \
        }                                                                      \
                                                                               \
        /*  Create a new list to be returned to the user.                    */\
        output = PyList_New(dim);                                              \
                                                                               \
        /*  Loop over the entries of the list, examine values for error, and */\
        /*  try to perform computations if legal inputs are given.           */\
        for (i=0; i<dim; ++i)                                                  \
        {                                                                      \
            /*  Extract the python object in the ith slot of the list.       */\
            ith_item = PyList_GET_ITEM(rho, i);                                \
                                                                               \
            /*  The list should be homogeneous with real numbers only. Check */\
            /*  each item and raise an error otherwise.                      */\
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))           \
            {                                                                  \
                PyErr_Format(PyExc_TypeError,                                  \
                             "\n\rError Encountered: rss_ringoccs\n"           \
                             "\r\tdiffrec.%s\n\n"                              \
                             "\rInput list must contain real numbers only.\n", \
                             VarToString(FuncName));                           \
                return NULL;                                                   \
            }                                                                  \
                                                                               \
            /*  Convert the python float to a C double.                      */\
            val = PyFloat_AsDouble(ith_item);                                  \
                                                                               \
            /*  Compute the value, set the item in the list, and move on to  */\
            /*  to the next one. As explained before, ## concatenates words  */\
            /*  so CName##_Double compiles as CName_Double.                  */\
            val = CName##_Double(val, a, b, F);                                \
            PyList_SET_ITEM(output, i, PyFloat_FromDouble(val));               \
        }                                                                      \
                                                                               \
        /*  The list created without errors, so just return it.              */\
        return output;                                                         \
    }                                                                          \
                                                                               \
    /*  If the input was not a list or a number, it must be a numpy array.   */\
    /*  Check this and return error otherwise.                               */\
    else if (!PyArray_Check(rho))                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rCould not parse inputs. Legal inputs are:\n"                    \
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"       \
            "\r\ta:     Positive constant (Float)\n"                           \
            "\r\tb:     Positive constant (Float) greater than a\n"            \
            "\r\tF      Positive constant (Float)\n\n"                         \
            "\rNotes:\n"                                                       \
            "\r\trho must be a non-empty one dimensional numpy array,\n"       \
            "\r\ta float, or an int.\n", VarToString(FuncName)                 \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Check the input array to make sure it is valid. It is important to   */\
    /*  check the number of dimensions BEFORE setting the dim variable. If   */\
    /*  the input array is empty, PyArray_DIMS(rho) will be a NULL pointer.  */\
    /*  Trying to access it with dim = PyArray_DIMS(rho)[0] will create a    */\
    /*  segmentation fault, crashing the Python interpreter.                 */\
    if (PyArray_NDIM((PyArrayObject *)rho) != 1)                               \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is not one-dimensional.\n",                   \
            VarToString(FuncName)                                              \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  Useful information about the data. The numpy API functions want a    */\
    /*  PyArrayObject pointer, so we need to cast rho as this type.          */\
    int typenum = PyArray_TYPE((PyArrayObject *)rho);                          \
    void *data  = PyArray_DATA((PyArrayObject *)rho);                          \
    dim         = PyArray_DIMS((PyArrayObject *)rho)[0];                       \
                                                                               \
    /*  This void pointer will point to the output data we'll create later.  */\
    void *y;                                                                   \
                                                                               \
    /*  If the input array is empty, return with an error.                   */\
    if (dim == 0)                                                              \
    {                                                                          \
        PyErr_Format(                                                          \
            PyExc_TypeError,                                                   \
            "\n\rError Encountered: rss_ringoccs\n"                            \
            "\r\tdiffrec.%s\n\n"                                               \
            "\rInput numpy array is empty.\n", VarToString(FuncName)           \
        );                                                                     \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    /*  For float, double, and long double precision (which corresponds to   */\
    /*  numpy.float32, numpy.float64, and numpy.float128, respectively) we   */\
    /*  have functions which compute at these respective precisions. Pass    */\
    /*  the values into the _array_from_four macro which will loop over the  */\
    /*  data and create a new output pointer, which y will then point to.    */\
    if (typenum == NPY_FLOAT)                                                  \
        _array_from_four(data, a, b, F, y, dim, float, float, CName##_Float);  \
    else if (typenum == NPY_DOUBLE)                                            \
        _array_from_four(data, a, b, F, y, dim, double,                        \
                         double, CName##_Double);                              \
    else if (typenum == NPY_LONGDOUBLE)                                        \
        _array_from_four(data, a, b, F, y, dim, long double,                   \
                         long double, CName##_Long_Double);                    \
                                                                               \
    /*  For all other types we try to convert to double and compute.         */\
    else                                                                       \
    {                                                                          \
                                                                               \
        /*  Try to convert the input numpy array to double and compute.      */\
        PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);          \
                                                                               \
        /*  If PyArray_FromObject failed, newrho should be NULL. Check this. */\
        if (!(newrho))                                                         \
        {                                                                      \
            PyErr_Format(PyExc_TypeError,                                      \
                        "\n\rError Encountered: rss_ringoccs\n"                \
                        "\r\tdiffrec.%s\n\n"                                   \
                        "\rInvalid data type for input array. Input should"    \
                        "\n\rbe a 1-dimensional array of real numbers.\n",     \
                        VarToString(FuncName));                                \
            return NULL;                                                       \
        }                                                                      \
        else                                                                   \
            _array_from_four(data, a, b, F, y, dim, double,                    \
                             double, CName##_Double);                          \
                                                                               \
        /*  Retrieve a pointer to the data inside of newrho.                 */\
        data = PyArray_DATA((PyArrayObject *)newrho);                          \
    }                                                                          \
                                                                               \
    /*  Create a numpy array from the newly computed data and set a capsule  */\
    /*  for it. This ensures when the corresponding Python variable is       */\
    /*  deleted or removed the memory allocated to the respective C pointer  */\
    /*  is freed. Skipping this will result in memory leaks!                 */\
    output  = PyArray_SimpleNewFromData(1, &dim, typenum, y);                  \
    capsule = PyCapsule_New(y, NULL, capsule_cleanup);                         \
                                                                               \
    /*  This frees the variable at the Python level once it's destroyed.     */\
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);                   \
                                                                               \
    /*  Return the results to Python.                                     */   \
    return Py_BuildValue("N", output);                                         \
}

RSS_RINGOCCSDiffModelingFourVars(gap_diffraction, Gap_Diffraction);
RSS_RINGOCCSDiffModelingFourVars(ringlet_diffraction, Ringlet_Diffraction);
RSS_RINGOCCSDiffModelingFourVarsREAL(ringlet_diffraction_phase,
                                     Ringlet_Diffraction_Phase);
RSS_RINGOCCSDiffModelingThreeVars(left_straightedge,
                                  Left_Straightedge_Diffraction);
RSS_RINGOCCSDiffModelingThreeVars(right_straightedge,
                                  Right_Straightedge_Diffraction);

static PyObject *square_wave_diffraction(PyObject *self, PyObject *args)
{
    /*  We'll need output and capsule for safely creating the output array    *
     *  and ensuring we don't have a memory leak. rho is the input variable.  */
    PyObject *output, *capsule, *rho;

    /*  Variables for the width of the waves and the Fresnel scale.           */
    double W, F;

    /*  The number of waves.                                                  */
    long N;

    /*  Variable for the size of the input array or list.                     */
    long dim;

    /*  Try to parse the user input, returning error if this fails.           */
    if (!PyArg_ParseTuple(args, "OddK", &rho, &W, &F, &N))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\tW:     Positive constant (Float)\n"
            "\r\tF:     Positive constant (Float)\n"
            "\r\tN:     Positive Integer (Int)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array."
        );
        return NULL;
    }

    /*  The width must be positive. Check this.                               */
    if (W <= 0.0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rWidth is not positive. (i.e. W<=0)\n"
        );
        return NULL;
    }

    /*  As always, the Fresnel scale is a positive value. Check.              */
    else if (F <= 0.0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rFresnel scale is not positive. (i.e. F<=0)\n"
        );
        return NULL;
    }

    /*  The number of waves to sum over must be positive as well.             */
    else if (N <= 0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rThe number of waves is not positive (i.e. N<=0)\n"
        );
        return NULL;
    }

    /*  If the user supplied a float or int, simply pass the inputs to the    *
     *  C function and return the value back to Python.                       */
    if (PyLong_Check(rho) || PyFloat_Check(rho))
    {
        /*  Convert a python float or int into a C double.                    */
        double val = PyFloat_AsDouble(rho);

        /*  The output is a complex double, so declare one.                   */
        complex double out;

        /*  Now get the computed value and convert the C complex double to a  *
         *  Python complex value and return.                                  */
        out = Square_Wave_Diffraction_Double(val, W, F, N);

        /*  To build a Python complex value we need to pass the real and      *
         *  imaginary parts of out into PyComplex_FromDoubles. We can easily  *
         *  access these values using creal and cimag, which are functions    *
         *  found in complex.h.                                               */
        return PyComplex_FromDoubles(creal(out), cimag(out));
    }

    /*  As of v1.3 we allow lists to be passed to these routines.             */
    else if (PyList_Check(rho))
    {
        /*  If the user passed a list, we'll need to return one. Create a     *
         *  variable for indexing over the input list.                        */
        long i;

        /*  Create another variable for indexing over the list. This will be  *
         *  the object corresponding to the value of the index i.             */
        PyObject *ith_item;

        /*  And lastly, a double for storing the input value and a complex    *
         *  double for storing the output.                                    */
        double val;
        complex double out;

        /*  Get the number of elements in the list.                           */
        dim = PyList_Size(rho);

        /*  If the input list is empty, return with an error.                 */
        if (dim == 0)
        {
            PyErr_Format(
                PyExc_TypeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.square_wave_diffraction\n\n"
                "\rInput list is empty.\n"
            );
            return NULL;
        }

        /*  Create a new list to be returned to the user.                     */
        output = PyList_New(dim);

        /*  Loop over the entries of the list, examine values for error, and  *
         *  try to perform computations if legal inputs are given.            */
        for (i=0; i<dim; ++i)
        {
            /*  Extract the python object in the ith slot of the list.        */
            ith_item = PyList_GET_ITEM(rho, i);

            /*  The list should be homogeneous with real numbers only. Check  *
             *  each item and raise an error otherwise.                       */
            if (!PyFloat_Check(ith_item) && !PyLong_Check(ith_item))
            {
                PyErr_Format(PyExc_TypeError,
                             "\n\rError Encountered: rss_ringoccs\n"
                             "\r\tdiffrec.%s\n\n"
                             "\rInput list must contain real numbers only.\n");
                return NULL;
            }

            /*  Convert the python float to a C double.                       */
            val = PyFloat_AsDouble(ith_item);

            /*  Compute the value, set the item in the list, and move on to   *
             *  to the next one. We'll need to convert the complex double     *
             *  to a Python complex value as well. Again, we'll use creal and *
             *  cimag, found in complex.h, to perform this task.              */
            out = Square_Wave_Diffraction_Double(val, W, F, N);

            PyList_SET_ITEM(output, i,
                            PyComplex_FromDoubles(creal(out), cimag(out)));
        }

        /*  The list created without errors, so just return it.               */
        return output;
    }

    /*  If the input is not a float, int, or list, it must be a numpy array.  *
     *  Check this and return with error if not.                              */
    else if (!PyArray_Check(rho))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rCould not parse inputs. Legal inputs are:\n"
            "\r\trho:   Numpy Array of positive real numbers (Floats)\n"
            "\r\ta:     Positive constant (Float)\n"
            "\r\tb:     Positive constant (Float) greater than a\n"
            "\r\tF      Positive constant (Float)\n\n"
            "\rNotes:\n"
            "\r\trho must be a non-empty one dimensional numpy array,\n"
            "\r\ta float, or an int.\n"
        );
        return NULL;
    }

    /*  Variables needed for the array data.                                  */
    void *data;
    complex double *T_hat;

    /*  Try to convert rho to double. Raise an error otherwise.               */
    PyObject *newrho = PyArray_FromObject(rho, NPY_DOUBLE, 1, 1);

    /*  If PyArray_FromObject failed, newrho is NULL. Check this.             */
    if (!(newrho))
    {
        PyErr_Format(PyExc_TypeError,
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\tdiffrec.square_wave_diffraction\n\n"
                    "\rInvalid data type for input array. Input should"
                    "\n\rbe a 1-dimensional array of real numbers.\n");
        return NULL;
    }

    /*  Check the inputs to make sure they're valid.                          */
    if (PyArray_NDIM((PyArrayObject *)newrho) != 1)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rInput numpy array is not one-dimensional.\n"
        );
        return NULL;
    }

    /*  Get the size of the input numpy array.                                */
    dim = PyArray_DIMS((PyArrayObject *)newrho)[0];

    /*  Check that the array isn't empty. Raise error otherwise.              */
    if (dim == 0)
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.square_wave_diffraction\n\n"
            "\rInput numpy array is empty.\n"
        );
        return NULL;
    }

    /*  Get a pointer to the actual data from the array. Allocate memory for  *
     *  the data of the output numpy array, which we'll call T_hat.           */
    data  = PyArray_DATA((PyArrayObject *)newrho);
    T_hat = (complex double *)malloc(dim*sizeof(complex double));

    /*  Variable for indexing over the array.                                 */
    long i;

    /*  Loop over the elements of the array and compute.                      */
    for (i=0; i<dim; ++i)
        T_hat[i] = Square_Wave_Diffraction_Double(((double *)data)[i], W, F, N);

    /*  Set the output and capsule, ensuring no memory leaks occur.           */
    output = PyArray_SimpleNewFromData(1, &dim, NPY_CDOUBLE, (void *)T_hat);
    capsule = PyCapsule_New((void *)T_hat, NULL, capsule_cleanup);

    /*  This frees the variable at the Python level once it's destroyed.      */
    PyArray_SetBaseObject((PyArrayObject *)output, capsule);

    /*  Return the results to Python.                                         */
    return Py_BuildValue("N", output);
}

static PyMethodDef diffrec_methods[] =
{
    {
        "square_wave_diffraction",
        square_wave_diffraction,
        METH_VARARGS,
        "Compute the normalized equivalent width of an array."
    },
    {
        "gap_diffraction",
        gap_diffraction,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.gap_diffraction\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of an annular gap in the plane.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, the inner radius of the annulus.\n\r\t\t"
        "b (float)\n\r\t\t\t"
        "A positive real number, the outter radius of the annulus.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.gap_diffraction(x, 45, 55, 0.05)"
    },
    {
        "ringlet_diffraction",
        ringlet_diffraction,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.gap_diffraction\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of a ringlet in the plane.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, the inner radius of the ring.\n\r\t\t"
        "b (float)\n\r\t\t\t"
        "A positive real number, the outter radius of the ring.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.ringlet_diffraction(x, a, b, F)"
    },
    {
        "ringlet_diffraction_phase",
        ringlet_diffraction_phase,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.ringlet_diffraction_phase\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the phase of a ringlet in the plane.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, the inner radius of the ring.\n\r\t\t"
        "b (float)\n\r\t\t\t"
        "A positive real number, the outter radius of the ring.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import special_functions\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> b = 55\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.ringlet_diffraction_phase(x, a, b, F)"
    },
    {
        "right_straightedge",
        right_straightedge,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.right_straightedge\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of a straightedge.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, starting point of the straightedge.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.right_straightedge(x, a, F)"
    },
    {
        "left_straightedge",
        left_straightedge,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.left_straightedge\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the diffraction pattern of a straightedge.\n\r\t"
        "Arguments:\n\r\t\t"
        "rho (numpy.ndarray):\n\r\t\t\t"
        "A numpy array of positive real numbers.\n\r\t\t"
        "a (float)\n\r\t\t\t"
        "A positive real number, starting point of the straightedge.\n\r\t\t"
        "F (float)\n\r\t\t\t"
        "A positive real number, the Fresnel scale (same units as rho).\n\r\t"
        "Outputs:\n\r\t\t"
        "T_hat (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of complex numbers equal to the diffraction pattern.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> x = numpy.arange(0,100,0.01)\n\r\t\t"
        ">>> a = 45\n\r\t\t"
        ">>> F = 0.05\n\r\t\t"
        ">>> y = diffrec.left_straightedge(x, a, F)"
    },
    {
        "fresnel_scale",
        fresnel_scale,
        METH_VARARGS,
        "\r\t"
        "Function:\n\r\t\t"
        "diffrec.fresnel_scale\n\r\t"
        "Purpose\n\r\t\t"
        "Compute the Fresnel scale.\n\r\t"
        "Arguments:\n\r\t\t"
        "Lambda (numpy.ndarray):\n\r\t\t\t"
        "Numpy array of positive real numbers, the wavelength.\n\r\t\t"
        "D (numpy.ndarray)\n\r\t\t\t"
        "Numpy array of real numbers, distance to ring intercept point.\n\r\t\t"
        "phi (numpy.ndarray)\n\r\t\t\t"
        "Numpy array of real numbers, the ring azimuth angle (radians).\n\r\t"
        "B (numpy.ndarray)\n\r\t\t\t"
        "Numpy array of real numbers, the ring opening angle (radians).\n\r\t"
        "Outputs:\n\r\t\t"
        "F (numpy.ndarray):\n\r\t\t\t"
        "The Fresnel scale.\n\r\t"
        "Example:\n\r\t\t"
        ">>> import numpy\n\r\t\t"
        ">>> import diffrec\n\r\t\t"
        ">>> lambda = 3.6e-5"
        ">>> phi = numpy.arange(0,3.14,0.01)\n\r\t\t"
        ">>> B = 0.3\n\r\t\t"
        ">>> D = 2.0e6\n\r\t\t"
        ">>> y = diffrec.fresnel_scale(lambda, D, phi, B)"
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "diffrec",
    NULL,
    -1,
    diffrec_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_diffraction_modeling(void)
{
    PyObject *m = PyModule_Create(&moduledef);
    if (!m) return NULL;

    import_array();

    return m;
}

