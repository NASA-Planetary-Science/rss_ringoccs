#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

const rssringoccs_FresnelNewtonTransform
rssringoccs_newton_transform_table[2][2][3] = {

    /*  Inverse transforms, integration is done with respect to rho0.         */
    {
        /*  Inverse transforms without normalizing by the window width.       */
        {
            rssringoccs_Fresnel_Transform_Newton,
            rssringoccs_Fresnel_Transform_Elliptical_Newton,
            rssringoccs_Fresnel_Transform_Perturbed_Newton
        },

        /*  Inverse transforms with normalizing by the window width.          */
        {
            rssringoccs_Fresnel_Transform_Normalized_Newton,
            rssringoccs_Fresnel_Transform_Normalized_Elliptical_Newton,
            rssringoccs_Fresnel_Transform_Normalized_Perturbed_Newton
        }
    },

    /*  Forward transforms, integration is done with respect to rho.          */
    {
        /*  Forward transform without normalizing by the window width.        */
        {
            NULL,
            NULL,
            NULL
        },

        /*  Forward transform with normalizing by the window width.           */
        {
            rssringoccs_Fresnel_Transform_Normalized_Forward_Newton,
            NULL,
            NULL
        }
    }
};
