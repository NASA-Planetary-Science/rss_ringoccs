
double Single_Slit_Fraunhofer_Diffraction(double x, double z, double a){
    return np.square(np.sinc(a*x/z))
}

double Double_Slit_Fraunhofer_Diffraction(double x, double z,
                                          double a, double d){
    f1 = np.square(np.sinc(a*x/z))
    f2 = np.square(np.sin(TWO_PI*d*x/z))
    f3 = 4.0*np.square(np.sin(np.pi*d*x/z))

    return f1*f2/f3
}