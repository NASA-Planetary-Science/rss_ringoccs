

import rss_ringoccs; import numpy; import matplotlib.pyplot as plt
geou = "/Users/ryan/Downloads/rss_uringoccs/output/E/G/VGR2_X43_E_URING_G_GEO_20210101_0001.TAB"
dlpu = "/Users/ryan/Downloads/rss_uringoccs/output/E/G/VGR2_X43_E_URING_G_DLP_005M_20210101_0001.TAB"
taug = "/Users/ryan/Downloads/rss_uringoccs/URINGS_GRESH1989_050M/RU1P2XGE.TAB"
datau = rss_ringoccs.tools.CSV_tools.GetUranusData(geou, dlpu)
rho_dg, tau_dg = numpy.loadtxt(taug, usecols=(0,1), delimiter=',').T
rng  = [47627, 47630]
peri = -0.3851087541499648

rece = rss_ringoccs.diffrec.DiffractionCorrection(datau, 0.05, psitype="ellipse", ecc=0.0, peri=0.0, rng=rng, res_factor=1.0)
plt.subplot(611)
plt.plot(rece.rho_km_vals, rece.tau_vals, label="ellipse ecc=0")
plt.plot(rho_dg, tau_dg, label="Gresh")
plt.xlim(47627, 47630)
plt.ylim(8, -1)
plt.legend()

rece = rss_ringoccs.diffrec.DiffractionCorrection(datau, 0.05, psitype="ellipse", ecc=1.09e-4, peri=peri, rng=rng, res_factor=1.0)
plt.subplot(612)
plt.plot(rece.rho_km_vals, rece.tau_vals, label="ellipse ecc=1.04e-4")
plt.plot(rho_dg, tau_dg, label="Gresh")
plt.xlim(47627, 47630)
plt.ylim(8, -1)
plt.legend()

rece = rss_ringoccs.diffrec.DiffractionCorrection(datau, 0.05, psitype="ellipse", ecc=1.09e-3, peri=peri, rng=rng, res_factor=1.0)
plt.subplot(613)
plt.plot(rece.rho_km_vals, rece.tau_vals, label="ellipse ecc=1.04e-3")
plt.plot(rho_dg, tau_dg, label="Gresh")
plt.xlim(47627, 47630)
plt.ylim(8, -1)
plt.legend()

rece = rss_ringoccs.diffrec.DiffractionCorrection(datau, 0.05, psitype="ellipse", ecc=1.09e-5, peri=peri, rng=rng, res_factor=1.0)
plt.subplot(614)
plt.plot(rece.rho_km_vals, rece.tau_vals, label="ellipse ecc=1.04e-5")
plt.plot(rho_dg, tau_dg, label="Gresh")
plt.xlim(47627, 47630)
plt.ylim(8, -1)
plt.legend()

rece = rss_ringoccs.diffrec.DiffractionCorrection(datau, 0.05, psitype="newtond", rng=rng, res_factor=1.0)
plt.subplot(615)
plt.plot(rece.rho_km_vals, rece.tau_vals, label="newtond")
plt.plot(rho_dg, tau_dg, label="Gresh")
plt.xlim(47627, 47630)
plt.ylim(8, -1)
plt.legend()

rece = rss_ringoccs.diffrec.DiffractionCorrection(datau, 0.05, psitype="newton", rng=rng, res_factor=1.0)
plt.subplot(616)
plt.plot(rece.rho_km_vals, rece.tau_vals, label="newton")
plt.plot(rho_dg, tau_dg, label="Gresh")
plt.xlim(47627, 47630)
plt.ylim(8, -1)
plt.legend()
