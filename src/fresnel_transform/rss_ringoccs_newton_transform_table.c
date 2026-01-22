#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

const rssringoccs_FresnelNewtonTransform
rssringoccs_newton_transform_table[7] = {
    rssringoccs_Fresnel_Transform_Newton_Riemann,
    rssringoccs_Fresnel_Transform_Newton_Filon01,
    rssringoccs_Fresnel_Transform_Newton_Filon11,
    rssringoccs_Fresnel_Transform_Newton_Filon02,
    rssringoccs_Fresnel_Transform_Newton_Filon12,
    rssringoccs_Fresnel_Transform_Elliptical_Newton,
    rssringoccs_Fresnel_Transform_Perturbed_Newton
};
