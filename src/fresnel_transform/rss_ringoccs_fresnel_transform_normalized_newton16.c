#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>
#include <libtmpl/include/tmpl_vec2.h>
#include <libtmpl/include/tmpl_vec3.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

#define C00 (+1.422222222222222222222220E+01)
#define C01 (-4.977777777777777777777780E+00)
#define C02 (+1.810101010101010101010100E+00)
#define C03 (-5.656565656565656565656570E-01)
#define C04 (+1.392385392385392385392390E-01)
#define C05 (-2.486402486402486402486400E-02)
#define C06 (+2.841602841602841602841600E-03)
#define C07 (-1.554001554001554001554000E-04)

#define C10 (+4.551111111111111111111110E+02)
#define C11 (-7.964444444444444444444440E+01)
#define C12 (+1.930774410774410774410770E+01)
#define C13 (-4.525252525252525252525250E+00)
#define C14 (+8.911266511266511266511270E-01)
#define C15 (-1.326081326081326081326080E-01)
#define C16 (+1.299018441875584732727590E-02)
#define C17 (-6.216006216006216006216010E-04)

#define C20 (-1.920285089443184681279920E+03)
#define C21 (+1.627833114638447971781310E+03)
#define C22 (-6.562984614397947731281060E+02)
#define C23 (+2.121325509058842392175730E+02)
#define C24 (-5.301925728592395259061930E+01)
#define C25 (+9.545521286473667426048380E+00)
#define C26 (-1.096277746944413611080280E+00)
#define C27 (+6.014297519059423821328580E-02)

#define C30 (-6.144912286218190980095740E+04)
#define C31 (+2.604532983421516754850090E+04)
#define C32 (-7.000516922024477580033140E+03)
#define C33 (+1.697060407247073913740580E+03)
#define C34 (-3.393232466299132965799630E+02)
#define C35 (+5.090944686119289293892470E+01)
#define C36 (-5.011555414603033650652700E+00)
#define C37 (+2.405719007623769528531430E-01)

#define C40 (+9.152430529492455418381340E+04)
#define C41 (-9.990973153047227121301200E+04)
#define C42 (+5.554688246272246272246270E+04)
#define C43 (-1.979804405679368642331610E+04)
#define C44 (+5.165923389093759464129830E+03)
#define C45 (-9.515565560365560365560370E+02)
#define C46 (+1.107794587853847113106370E+02)
#define C47 (-6.130901964976039050113120E+00)

#define C50 (+2.928777769437585733882030E+06)
#define C51 (-1.598555704487556339408190E+06)
#define C52 (+5.925000796023729357062690E+05)
#define C53 (-1.583843524543494913865280E+05)
#define C54 (+3.306190969020006057043090E+04)
#define C55 (-5.074968298861632194965530E+03)
#define C56 (+5.064203830189015374200560E+02)
#define C57 (-2.452360785990415620045250E+01)

#define C60 (-2.067472141595140113658630E+06)
#define C61 (+2.529970186033313737017440E+06)
#define C62 (-1.665160880084656084656080E+06)
#define C63 (+6.973441364530668234371940E+05)
#define C64 (-1.967287208308837938467570E+05)
#define C65 (+3.780976423280423280423280E+04)
#define C66 (-4.515678262982559278855570E+03)
#define C67 (+2.540788210856359004507150E+02)

#define C70 (-6.615910853104448363707620E+07)
#define C71 (+4.047952297653301979227910E+07)
#define C72 (-1.776171605423633156966490E+07)
#define C73 (+5.578753091624534587497550E+06)
#define C74 (-1.259063813317656280619240E+06)
#define C75 (+2.016520759082892416225750E+05)
#define C76 (-2.064310063077741384619690E+04)
#define C77 (+1.016315284342543601802860E+03)

#define C80 (+2.453459187436857870720300E+07)
#define C81 (-3.191451913682147756221830E+07)
#define C82 (+2.312000968917107583774250E+07)
#define C83 (-1.086892688981383499902020E+07)
#define C84 (+3.407389116676464824612970E+06)
#define C85 (-6.993254258100277147896200E+05)
#define C86 (+8.705872617675876935136190E+04)
#define C87 (-5.034902592872539962487050E+03)

#define C90 (+7.851069399797945186304970E+08)
#define C91 (-5.106323061891436409954930E+08)
#define C92 (+2.466134366844914756025870E+08)
#define C93 (-8.695141511851067999216150E+07)
#define C94 (+2.180729034672937487752300E+07)
#define C95 (-3.729735604320147812211300E+06)
#define C96 (+3.979827482366115170347970E+05)
#define C97 (-2.013961037149015984994820E+04)

#define C100 (-1.560092373087517146776410E+08)
#define C101 (+2.103734402449131883205960E+08)
#define C102 (-1.616015015691697931697930E+08)
#define C103 (+8.210883624042043005005970E+07)
#define C104 (-2.812659180891097928134970E+07)
#define C105 (+6.280279259182299182299180E+06)
#define C106 (-8.312494380087292679885270E+05)
#define C107 (+5.019326503440044180784920E+04)

#define C110 (-4.992295593880054869684500E+09)
#define C111 (+3.365975043918611013129530E+09)
#define C112 (-1.723749350071144460477790E+09)
#define C113 (+6.568706899233634404004770E+08)
#define C114 (-1.800101875770302674006380E+08)
#define C115 (+3.349482271563892897226230E+07)
#define C116 (-3.799997430897048082233270E+06)
#define C117 (+2.007730601376017672313970E+05)

#define C120 (+4.998759507035704487556340E+08)
#define C121 (-6.894840699359592396629430E+08)
#define C122 (+5.500202466989129389129390E+08)
#define C123 (-2.945977389726371296741670E+08)
#define C124 (+1.078825598938257901220860E+08)
#define C125 (-2.603646138219706219706220E+07)
#define C126 (+3.736714365037541333837630E+06)
#define C127 (-2.410783461314542796024280E+05)

#define C130 (+1.599603042251425436018030E+10)
#define C131 (-1.103174511897534783460710E+10)
#define C132 (+5.866882631455071348404680E+09)
#define C133 (-2.356781911781097037393330E+09)
#define C134 (+6.904483833204850567813530E+08)
#define C135 (-1.388611273717176650509980E+08)
#define C136 (+1.708212281160018895468630E+07)
#define C137 (-9.643133845258171184097110E+05)

#define C140 (-6.303854353700198762632620E+08)
#define C141 (+8.825396095180278267685680E+08)
#define C142 (-7.220778623329318582651920E+08)
#define C143 (+4.011543679627399212584400E+08)
#define C144 (-1.542901415241307389455540E+08)
#define C145 (+3.967460782049076144314240E+07)
#define C146 (-6.171605660965229557822150E+06)
#define C147 (+4.408289757832306827015820E+05)

#define C150 (-2.017233393184063604042440E+10)
#define C151 (+1.412063375228844522829710E+10)
#define C152 (-7.702163864884606488162040E+09)
#define C153 (+3.209234943701919370067520E+09)
#define C154 (-9.874569057544367292515440E+08)
#define C155 (+2.115979083759507276967590E+08)
#define C156 (-2.821305445012676369290130E+07)
#define C157 (+1.763315903132922730806330E+06)

#define COEFF_POLY_EVAL(z, N) ( \
    C##N##0*z[0] +\
        C##N##1*z[1] +\
            C##N##2*z[2] +\
                C##N##3*z[3] +\
                    C##N##4*z[4] +\
                        C##N##5*z[5] +\
                            C##N##6*z[6] +\
                                C##N##7*z[7]\
)

#define ODD_POLY_EVAL(z) \
    z * (\
        coeffs[0] + z##sq * (\
            coeffs[2] + z##sq * (\
                coeffs[4] + z##sq * (\
                    coeffs[6] + z##sq * (\
                        coeffs[8] + z##sq * (\
                            coeffs[10] + z##sq * (\
                                coeffs[12] + z##sq * coeffs[14]\
                            )\
                        )\
                    )\
                )\
            )\
        )\
    )

#define EVEN_POLY_EVAL(z) \
    z##sq * (\
        coeffs[1] + z##sq * (\
            coeffs[3] + z##sq * (\
                coeffs[5] + z##sq * (\
                    coeffs[7] + z##sq * (\
                        coeffs[9] + z##sq * (\
                            coeffs[11] + z##sq * (\
                                coeffs[13] + z##sq * coeffs[15]\
                            )\
                        )\
                    )\
                )\
            )\
        )\
    )

void
rssringoccs_Fresnel_Transform_Normalized_Newton16(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const double * TMPL_RESTRICT const x_arr,
    const double * TMPL_RESTRICT const w_func,
    size_t nw_pts,
    size_t center
)
{
    size_t n, ind[16];

    double scale_factor;
    tmpl_ComplexDouble w_exp_minus_psi_left, w_exp_minus_psi_right;
    tmpl_ComplexDouble T_left, T_right, integrand;

    double coeffs[16], diff[8], mean[8], psi[16];

    const size_t shift = nw_pts >> 3;
    const size_t l_ind = center - nw_pts;
    const size_t r_ind = center + nw_pts;

    const double width_actual = 16.0 * tau->dx_km * TMPL_CAST(shift, double);
    const double rcpr_width_actual = 1.0 / width_actual;

    tmpl_ComplexDouble norm = tmpl_CDouble_One;
    tau->T_out[center] = tau->T_in[center];

    ind[0] = center - 8 * shift;
    ind[1] = center - 7 * shift;
    ind[2] = center - 6 * shift;
    ind[3] = center - 5 * shift;
    ind[4] = center - 4 * shift;
    ind[5] = center - 3 * shift;
    ind[6] = center - 2 * shift;
    ind[7] = center - shift;
    ind[8] = center + shift;
    ind[9] = center + 2 * shift;
    ind[10] = center + 3 * shift;
    ind[11] = center + 4 * shift;
    ind[12] = center + 5 * shift;
    ind[13] = center + 6 * shift;
    ind[14] = center + 7 * shift;
    ind[15] = center + 8 * shift;

    for (n = 0; n < 16; ++n)
    {
        const tmpl_TwoVectorDouble rho0 = tmpl_2DDouble_Polard(
            tau->rho_km_vals[ind[n]], tau->phi_deg_vals[ind[n]]
        );

        const tmpl_TwoVectorDouble rho = tmpl_2DDouble_Polard(
            tau->rho_km_vals[center], tau->phi_deg_vals[ind[n]]
        );

        const tmpl_ThreeVectorDouble R = tmpl_3DDouble_Rect(
            tau->rx_km_vals[ind[n]],
            tau->ry_km_vals[ind[n]],
            tau->rz_km_vals[ind[n]]
        );

        psi[n] = tmpl_Double_Stationary_Cyl_Fresnel_Psi(
            tau->k_vals[center], &rho, &rho0, &R, tau->EPS, tau->toler
        );
    }

    diff[0] = psi[8] - psi[7];
    diff[1] = psi[9] - psi[6];
    diff[2] = psi[10] - psi[5];
    diff[3] = psi[11] - psi[4];
    diff[4] = psi[12] - psi[3];
    diff[5] = psi[13] - psi[2];
    diff[6] = psi[14] - psi[1];
    diff[7] = psi[15] - psi[0];

    mean[0] = (psi[8] + psi[7]) * 0.5;
    mean[1] = (psi[9] + psi[6]) * 0.5;
    mean[2] = (psi[10] + psi[5]) * 0.5;
    mean[3] = (psi[11] + psi[4]) * 0.5;
    mean[4] = (psi[12] + psi[3]) * 0.5;
    mean[5] = (psi[13] + psi[2]) * 0.5;
    mean[6] = (psi[14] + psi[1]) * 0.5;
    mean[7] = (psi[15] + psi[0]) * 0.5;

    coeffs[0] = COEFF_POLY_EVAL(diff, 0);
    coeffs[1] = COEFF_POLY_EVAL(mean, 1);
    coeffs[2] = COEFF_POLY_EVAL(diff, 2);
    coeffs[3] = COEFF_POLY_EVAL(mean, 3);
    coeffs[4] = COEFF_POLY_EVAL(diff, 4);
    coeffs[5] = COEFF_POLY_EVAL(mean, 5);
    coeffs[6] = COEFF_POLY_EVAL(diff, 6);
    coeffs[7] = COEFF_POLY_EVAL(mean, 7);
    coeffs[8] = COEFF_POLY_EVAL(diff, 8);
    coeffs[9] = COEFF_POLY_EVAL(mean, 9);
    coeffs[10] = COEFF_POLY_EVAL(diff, 10);
    coeffs[11] = COEFF_POLY_EVAL(mean, 11);
    coeffs[12] = COEFF_POLY_EVAL(diff, 12);
    coeffs[13] = COEFF_POLY_EVAL(mean, 13);
    coeffs[14] = COEFF_POLY_EVAL(diff, 14);
    coeffs[15] = COEFF_POLY_EVAL(mean, 15);

    for (n = 0; n < nw_pts; ++n)
    {
        const double x = x_arr[n] * rcpr_width_actual;
        const double xsq = x * x;

        const double psi_odd = ODD_POLY_EVAL(x);
        const double psi_even = EVEN_POLY_EVAL(x);
        const double psi_left = psi_even + psi_odd;
        const double psi_right = psi_even - psi_odd;

        w_exp_minus_psi_left = tmpl_CDouble_Polar(w_func[n], -psi_left);
        w_exp_minus_psi_right = tmpl_CDouble_Polar(w_func[n], -psi_right);

        T_left = tau->T_in[l_ind + n];
        T_right = tau->T_in[r_ind - n];

        tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_left);
        tmpl_CDouble_AddTo(&norm, &w_exp_minus_psi_right);

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_left, T_left);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);

        integrand = tmpl_CDouble_Multiply(w_exp_minus_psi_right, T_right);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
    }

    scale_factor = tmpl_Double_Rcpr_Sqrt_Two / tmpl_CDouble_Abs(norm);
    integrand = tmpl_CDouble_Rect(scale_factor, scale_factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(integrand, tau->T_out[center]);
}
