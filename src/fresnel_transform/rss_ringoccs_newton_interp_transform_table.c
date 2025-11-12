#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

const rssringoccs_FresnelTransform
rssringoccs_newton_interp_transform_table[6] = {
    rssringoccs_Fresnel_Transform_Newton4,
    rssringoccs_Fresnel_Transform_Newton8,
    rssringoccs_Fresnel_Transform_Newton16,
    NULL,
    NULL,
    NULL
};
