#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <stddef.h>

const rssringoccs_FresnelTransform
rssringoccs_newton_interp_transform_table[2][2][6] = {

    /*  Inverse transforms, integration is done with respect to rho0.         */
    {
        /*  Inverse transforms without normalizing by the window width.       */
        {
            rssringoccs_Fresnel_Transform_Newton4,
            NULL,
            NULL,
            NULL,
            NULL,
            NULL
        },

        /*  Inverse transforms with normalizing by the window width.          */
        {
            rssringoccs_Fresnel_Transform_Normalized_Newton4,
            NULL,
            rssringoccs_Fresnel_Transform_Normalized_Newton8,
            NULL,
            NULL,
            NULL
        }
    },

    /*  Forward transforms, integration is done with respect to rho.          */
    {
        /*  Forward transform without normalizing by the window width.        */
        {
            NULL,
            NULL,
            NULL,
            NULL,
            NULL,
            NULL
        },

        /*  Forward transform with normalizing by the window width.           */
        {
            NULL,
            NULL,
            NULL,
            NULL,
            NULL,
            NULL
        }
    }
};
