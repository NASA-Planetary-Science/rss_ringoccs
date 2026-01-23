#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <libtmpl/include/helper/tmpl_array_size.h>
#include <stddef.h>

static const char * const rssringoccs_rev_number_vals[84] = {
    "007RI", "007RI", "008RI", "008RI",
    "009RI", "010RI", "010RI", "011RI",
    "012RI", "012RI", "013RI", "014RI",
    "028RI", "028RI", "044RI", "046RI",
    "053RI", "054RI", "056RI", "057RI",
    "058RI", "060RI", "063RI", "064RI",
    "067RI", "079RI", "081RI", "082RI",
    "084RI", "089RI", "123RI", "123RI",
    "125RI", "125RI", "133RI", "133RI",
    "137RI", "137RI", "167RI", "168RI",
    "169RI", "170RI", "174RI", "179RI",
    "180RI", "189RI", "190RI", "191RI",
    "193RI", "194RI", "196RI", "197RI",
    "236RI", "237RI", "237RI", "238RI",
    "247RI", "248RI", "250RI", "251RI",
    "253RI", "255RI", "256RI", "257RI",
    "266RI", "268RI", "270RI", "273RI",
    "273RI", "273RI", "PERIO", "274RI",
    "274RI", "275RI", "275RI", "275RI",
    "276RI", "278RI", "278RI", "280RI",
    "280RI", "282RI", "284RI", "284RI"
};

static const unsigned int rssringoccs_year_vals[84] = {
    2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005,
    2005, 2006, 2006, 2007, 2007, 2007, 2007, 2008, 2008, 2008, 2008,
    2008, 2008, 2008, 2008, 2008, 2008, 2008, 2008, 2009, 2009, 2010,
    2010, 2010, 2010, 2010, 2010, 2012, 2012, 2012, 2012, 2012, 2013,
    2013, 2013, 2013, 2013, 2013, 2013, 2013, 2013, 2016, 2016, 2016,
    2016, 2016, 2016, 2016, 2016, 2016, 2017, 2017, 2017, 2017, 2017,
    2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017,
    2017, 2017, 2017, 2017, 2017, 2017, 2017
};

static const unsigned int rssringoccs_doy_vals[84] = {
    123, 123, 141, 141, 159, 177, 177, 196, 214, 214, 232, 248, 258,
    259, 130, 162, 337, 353,  15,  27,  39,  62,  92, 102, 130, 217,
    232, 239, 254, 291, 359, 360,  26,  27, 169, 170, 245, 245, 156,
    180, 204, 225, 315,  18,  31, 130, 140, 151, 175, 187, 220, 244,
    158, 182, 182, 206, 307, 317, 333, 340, 354,   3,  10,  17,  81,
     96, 110, 129, 129, 129, 135, 135, 136, 142, 142, 142, 148, 161,
    161, 174, 174, 187, 200, 200
};

const char *
rssringoccs_Date_to_Rev(unsigned int year, unsigned int doy)
{
    size_t number_of_revs, n;

    /*  There are no Cassini data sets before 2005.                           */
    if (year < 2005)
        return NULL;

    /*  The Cassini grand finale was in 2017, no data sets after that year.   */
    if (year > 2017)
        return NULL;

    /*  Leap year has 366 days. Any day-of-year beyond 366 is invalid.        */
    if (doy > 366)
        return NULL;

    /*  Total number of revs that include ring occultations.                  */
    number_of_revs = TMPL_ARRAY_SIZE(rssringoccs_rev_number_vals);

    /*  Loop through and see if the given day / year matches a rev number.    */
    for (n = 0; n < number_of_revs; ++n)
    {
        if (year != rssringoccs_year_vals[n])
            continue;

        if (doy != rssringoccs_doy_vals[n])
            continue;

        return rssringoccs_rev_number_vals[n];
    }

    return NULL;
}
