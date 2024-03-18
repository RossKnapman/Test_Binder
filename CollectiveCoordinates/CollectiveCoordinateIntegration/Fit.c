/*
C functions to calculate d/dt(R) and d/dt(eta), as this is faster than using Python functions (approximately halved
the time required to obtain the phase diagram in my testing)

Called by CCIntegration.sage

Must be compiled prior to running CCIntegration.sage with
gcc -fPIC -shared -o c_fit.o Fit.c
*/

#include <stdlib.h>
#include <math.h>


////////////////////////////
// Define Fits to Be Used //
////////////////////////////

double * inverse_fourth_order_fit_array(const double * arr, int N, double a, double b, double c, double d, double e, double f)
{
    double * result = (double *)malloc(sizeof(double) * N);
    for (int i=0; i<N; i++)
    {
        double R = arr[i];
        result[i] = a/(R*R*R*R) + b/(R*R*R) + c/(R*R) + d/R + f + e*R;
    }
    return result;
}

double inverse_fourth_order_fit(const double R, double a, double b, double c, double d, double e, double f)
{
    return (1/(R*R)) * (a/(R*R) + b/R + c) + d/R + f + e*R;
}

double linear_fit(const double R, double m, double c)
{
    return m*R + c;
}

double inverse_linear_fit(const double R, double a, double b, double c, double d)
{
    return a/(R-c) + b*(R-c) + d;
}

double inverse_quadratic_linear_fit(const double R, double a, double b, double c, double d)
{
    return a/(R*R) + b/R + c + d*R;
}

double quadratic_fit(const double R, double a, double b, double c)
{
    return a*R*R + b*R + c;
}


////////////////////////////////////////
// Function Returning dR/dt and dη/dt //
////////////////////////////////////////

double * time_derivatives(double t, double R, double eta, double alpha, double Ez, double Bz,
    double Gamma11ParamA, double Gamma11ParamB, double Gamma11ParamC, double Gamma11ParamD,
    double Gamma22ParamA, double Gamma22ParamB,
    double G12FitParamA, double G12FitParamB,
    double F_RexFitParamA, double F_RexFitParamB, double F_RexFitParamC, double F_RexFitParamD, double F_RexFitParamE, double F_RexFitParamF,
    double Xi1FitParamA, double Xi1FitParamB, double Xi1FitParamC, double Xi1FitParamD, double Xi1FitParamE, double Xi1FitParamF,
    double Xi2FitParamA, double Xi2FitParamB, double Xi2FitParamC, double Xi2FitParamD, double Xi2FitParamE, double Xi2FitParamF,
    double Xi3FitParamA, double Xi3FitParamB, double Xi3FitParamC, double Xi3FitParamD, double Xi3FitParamE, double Xi3FitParamF,
    double Xi4FitParamA, double Xi4FitParamB)
{

    // t is time
    // R is radius
    // eta is helicity
    // alpha is Gilbert damping constant
    // Ez is electric field
    // Bz is magnetic field
    // Gamma11<A,B,C,D> are the fit parameters for the dissipative tensor Γ_RR
    // Gamma22<A,B> are the fit parameters for the dissipative tensor Γ_ηη
    // F_RexFitParam<A,B,C,D,E,F> are the fit parameters for -dU_ex/DR
    // Xi1FitParams<A,B,C,D,E,F> are the fit parameters for the integral of cos(2θ) dθ / dR over ρ
    // X12FitParams<A,B,C,D,E,F> are the fit paarameters for the integral of ρ * d^2θ / dRdρ over ρ
    // Xi3FitParams<A,B,C,D,E,F> are the fit parameters for the integral of cos(θ)sin(θ) over ρ
    // Xi4FitParams<A,B> are the fit parameters for the integral of ρ * dθ / dρ over ρ

    // Obtain the quantities described above from their fits
    double Gamma11 = alpha * inverse_linear_fit(R, Gamma11ParamA, Gamma11ParamB, Gamma11ParamC, Gamma11ParamD);
    double Gamma22 = alpha * linear_fit(R, Gamma22ParamA, Gamma22ParamB);
    double G12 = linear_fit(R, G12FitParamA, G12FitParamB);
    double F_Rex = inverse_fourth_order_fit(R, F_RexFitParamA, F_RexFitParamB, F_RexFitParamC, F_RexFitParamD, F_RexFitParamE, F_RexFitParamF);
    double Xi1 = inverse_fourth_order_fit(R, Xi1FitParamA, Xi1FitParamB, Xi1FitParamC, Xi1FitParamD, Xi1FitParamE, Xi1FitParamF);
    double Xi2 = inverse_fourth_order_fit(R, Xi2FitParamA, Xi2FitParamB, Xi2FitParamC, Xi2FitParamD, Xi2FitParamE, Xi2FitParamF);
    double Xi3 = inverse_fourth_order_fit(R, Xi3FitParamA, Xi3FitParamB, Xi3FitParamC, Xi3FitParamD, Xi3FitParamE, Xi3FitParamF);
    double Xi4 = linear_fit(R, Xi4FitParamA, Xi4FitParamB);

    double F_R = F_Rex - Bz*G12 + Ez*cos(eta)*(Xi1 + Xi2);
    double F_eta = -Ez*sin(eta)*(Xi3 + Xi4);

    // The prefactor for both \dot{R} and \dot{\eta}
    double prefactor = 1 / (G12*G12 + Gamma11*Gamma22);

    double * return_array = (double *)malloc(sizeof(double) * 2);

    return_array[0] = prefactor * (Gamma22*F_R + G12*F_eta);
    return_array[1] = prefactor * (Gamma11*F_eta - G12*F_R);

    return return_array;

}
