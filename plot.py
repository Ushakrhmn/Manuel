import numpy as np
import matplotlib.pyplot as plt

# Load data from file (assumes whitespace-separated values)
data = np.loadtxt("Prob_NOvA_NH_MINOS.dat")

# Column indices: 0-based
energy = data[:, 0]  # 1st column: Energy in GeV
p_mumu = data[:, 2]  # 3rd column: P(nu_mu → nu_mu)
data_matter = np.loadtxt("Prob_NOvA_NH_MINOS.dat")

# Column indices: 0-based
energy_matter = data_matter[:, 0]  # 1st column: Energy in GeV
p_mumu_matter = data_matter[:, 2]  # 3rd column: P(nu_mu → nu_mu)
# Plot
plt.figure(figsize=(8, 5))
plt.plot(energy, p_mumu, label=r"$P(\nu_\mu \rightarrow \nu_\mu), vacuum$")
plt.plot(energy_matter, p_mumu_matter, label=r"$P(\nu_\mu \rightarrow \nu_\mu), matter")
plt.xlabel("Energy [GeV]")
plt.ylabel("Survival Probability")
plt.title(r"New Physics probability for MINOS")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

