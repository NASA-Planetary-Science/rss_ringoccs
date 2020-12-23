#include <stdlib.h>
#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>

void rssringoccs_Set_Tau_Error_Message(const char *mess,
                                       rssringoccs_TAUObj *tau)
{
    long length = strlen(mess);
    tau->error_message = malloc(sizeof(*tau->error_message) * length);
    strcpy(tau->error_message, mess);
    tau->error_occurred = rssringoccs_True;
}
