#include <rss_ringoccs/include/rss_ringoccs_history.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

int main(void)
{
	rssringoccs_GeoCSV geo;
	rssringoccs_GeoCSV_Init(&geo);
	rssringoccs_GeoCSV_Write_History(&geo, "my_bogus_file");
	rssringoccs_History_Print(geo.history);
	rssringoccs_GeoCSV_Destroy_Members(&geo);
	return 0;
}
