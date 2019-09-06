void Get_Float_Array(float *x, float *y, long dim, float (*f)(float))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Double_Array(double *x, double *y, long dim, double (*f)(double))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}

void Get_Long_Double_Array(long double *x, long double *y, long dim,
                           long double (*f)(long double))
{
    long i;
    for (i=0; i<dim; ++i){
        y[i] = (*f)(x[i]);
    }
}