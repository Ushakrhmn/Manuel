import numpy as np
import matplotlib.pyplot as plt

# Load data from file (assumes whitespace-separated values)
data_new = np.loadtxt("Prob_T2K_NH_dcp-90_matter_newphysics.dat")

# Column indices: 0-based
energy = data_new[:, 0]  # 1st column: Energy in GeV
p_mue_new = data_new[:, 1]  # 3rd column: P(nu_mu → nu_e)
p_muebar_new = data_new[:, 2]  # 3rd column: P(nu_mu_bar → nu_e_bar)
data_SI = np.loadtxt("Prob_T2K_NH_dcp-90_matter_SI.dat")

# Column indices: 0-based
energy_SI = data_SI[:, 0]  # 1st column: Energy in GeV
p_mue_SI = data_SI[:, 1]  # 3rd column: P(nu_mu → nu_e)
p_muebar_SI = data_SI[:, 2]  # 3rd column: P(nu_mu_bar → nu_e_bar)

# Plot
plt.figure(figsize=(8, 5))

plt.plot(energy_SI, p_mue_SI, '-', color='blue', label=r"$P(\nu_\mu \rightarrow \nu_e)$, SI")
plt.plot(energy, p_mue_new, '--', color='blue', label=r"$P(\nu_\mu \rightarrow \nu_e)$, new")
plt.plot(energy, p_muebar_SI, '-', color='red', label=r"$P(\bar{\nu}_\mu \rightarrow \bar{\nu}_e)$, SI")
plt.plot(energy, p_muebar_new, '--', color='red', label=r"$P(\bar{\nu}_\mu \rightarrow \bar{\nu}_e)$, new")


plt.ylim(0, 0.2)
plt.xlabel("Energy [GeV]")
plt.ylabel(r"Probability")
plt.title(r"T2K, NH, $\delta_{CP}=-90^\circ$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

