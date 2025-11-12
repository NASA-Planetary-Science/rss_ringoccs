#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

const rssringoccs_FresnelNewtonTransform
rssringoccs_newton_transform_table[5] = {
    rssringoccs_Fresnel_Transform_Newton_Riemann,
    rssringoccs_Fresnel_Transform_Newton_Linear_Filon,
    NULL,
    rssringoccs_Fresnel_Transform_Elliptical_Newton,
    rssringoccs_Fresnel_Transform_Perturbed_Newton
};
