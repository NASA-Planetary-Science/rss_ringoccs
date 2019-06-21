import sys

sys.path.append("../")

import rss_ringoccs
import spiceypy


geo007 = "/Volumes/rmaguire002/Research/TC2017/data/Archived_Cassini_RSS_RingOccs_2018/data/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO.TAB"
cal007 = "/Volumes/rmaguire002/Research/TC2017/data/Archived_Cassini_RSS_RingOccs_2018/data/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_CAL.TAB"
dlp007 = "/Volumes/rmaguire002/Research/TC2017/data/Archived_Cassini_RSS_RingOccs_2018/data/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_DLP_500M.TAB"
tau007 = "/Volumes/rmaguire002/Research/TC2017/data/Archived_Cassini_RSS_RingOccs_2018/data/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_TAU_01KM.TAB"

data007 = rss_ringoccs.tools.CSV_tools.ExtractCSVData(geo007, cal007, dlp007, tau=tau007, verbose=False)
rec007 = rss_ringoccs.diffrec.diffraction_correction.DiffractionCorrection(data007, 1.0, psitype="fresnel")
oet_spm = rec007.t_oet_spm_vals
set_spm = rec007.t_set_spm_vals
oet_et = rss_ringoccs.tools.spm_to_et(oet_spm, 123, 2005, kernels='../tables/e2e_kernels.ker')
set_et = rss_ringoccs.tools.spm_to_et(set_spm, 123, 2005, kernels='../tables/e2e_kernels.ker')

Sun2Cassini0, ltimeC0 = spiceypy.spkezp(-82, set_et[0], 'J2000', 'XCN', 0)
Sun2Cassini1, ltimeC1 = spiceypy.spkezp(-82, set_et[-1], 'J2000', 'XCN', 0)
Sun2Earth0, ltimeE0 = spiceypy.spkezp(399, oet_et[0], 'J2000', 'XCN', 0)
Sun2Earth1, ltimeE1 = spiceypy.spkezp(399, oet_et[-1], 'J2000', 'XCN', 0)

e2c0 = Sun2Cassini0 - Sun2Earth0
e2c1 = Sun2Cassini1 - Sun2Earth1

ltime0 = ltimeC0 - ltimeE0
ltime1 = ltimeC1 - ltimeE1

print(spiceypy.vnorm(e2c1-e2c0))

e2c0, ltime0 = spiceypy.spkezp(-82, oet_et[0], 'J2000', 'XLT+S', 399)
e2c1, ltime1 = spiceypy.spkezp(-82, oet_et[-1], 'J2000', 'XLT+S', 399)

print(spiceypy.vnorm(e2c1-e2c0))
