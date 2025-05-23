import ctypes
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import faulthandler
faulthandler.enable()
# Load the shared library
lib = ctypes.CDLL("./libchi2.so")

# Define function signatures
#lib.init_globes.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.compute_chi2.argtypes = [ctypes.c_double] * 7
lib.compute_chi2.restype = ctypes.c_double
lib.compute_chi2_ih.argtypes = [ctypes.c_double] * 7
lib.compute_chi2_ih.restype = ctypes.c_double
lib.free_globes.argtypes = []

# Initialize GLoBES with both NOvA appearance and disappearance experiments
#lib.init_globes(b"2024_nova_app.glb", b"2024_nova_disapp.glb")

def chi2_objective_nh(params):
    return lib.compute_chi2(*params)

def chi2_objective_ih(params):
    return lib.compute_chi2_ih(*params)

# Initial guesses for NH and IH
x0_nh = [
    np.arcsin(np.sqrt(0.02237)),   # theta13
    np.arcsin(np.sqrt(0.52)),      # theta23
    0.0,                           # delta_cp
    2.528e-3,                      # dm31
    5.0,                           # mu
    10.0,                          # N
    0.2                            # m0
]

x0_ih = [
    np.arcsin(np.sqrt(0.02259)),   # theta13
    np.arcsin(np.sqrt(0.52)),      # theta23
    0.0,                           # delta_cp
    -2.510e-3,                     # dm32 (IH)
    5.0,                           # mu
    10.0,                          # N
    0.2                            # m0
]

# Bounds for all parameters
bounds_nh = [
    (np.arcsin(np.sqrt(0.02237 - 3*0.00066)), np.arcsin(np.sqrt(0.02237 + 3*0.00066))),  # theta13
    (np.arcsin(np.sqrt(0.43)), np.arcsin(np.sqrt(0.61))),                                # theta23
    (-np.pi, np.pi),                                                                      # delta_cp
    (2.528e-3 - 3*0.000031, 2.528e-3 + 3*0.000031),                                        # dm31
    (1.0, 1000.0),                                                                         # mu
    (2.0, 1000.0),                                                                         # N
    (0.0, 0.4)                                                                             # m0
]

bounds_ih = [
    (np.arcsin(np.sqrt(0.02259 - 3*0.00065)), np.arcsin(np.sqrt(0.02259 + 3*0.00065))),  # theta13
    (np.arcsin(np.sqrt(0.43)), np.arcsin(np.sqrt(0.61))),                                # theta23
    (-np.pi, np.pi),                                                                      # delta_cp
    (-2.510e-3 - 3*0.000031, -2.510e-3 + 3*0.000031),                                      # dm32
    (1.0, 1000.0),                                                                         # mu
    (2.0, 1000.0),                                                                         # N
    (0.0, 0.4)                                                                             # m0
]

# Optimize to get minimum chi2 for NH and IH
res_nh = minimize(chi2_objective_nh, x0_nh, bounds=bounds_nh, method='L-BFGS-B')
res_ih = minimize(chi2_objective_ih, x0_ih, bounds=bounds_ih, method='L-BFGS-B')

min_chi2_global = min(res_nh.fun, res_ih.fun)

# Scan over delta_cp and theta23 with marginalization (NH)
theta23_vals = np.linspace(np.arcsin(np.sqrt(0.43)), np.arcsin(np.sqrt(0.61)), 30)
dcp_vals = np.linspace(-np.pi, np.pi, 30)
contour_nh = np.zeros((len(theta23_vals), len(dcp_vals)))

for i, t23 in enumerate(theta23_vals):
    for j, dcp in enumerate(dcp_vals):
        def marginalized(params):
            return chi2_objective_nh([params[0], t23, dcp, params[1], params[2], params[3], params[4]])

        x0 = [x0_nh[0], x0_nh[3], x0_nh[4], x0_nh[5], x0_nh[6]]
        bounds = [bounds_nh[0], bounds_nh[3], bounds_nh[4], bounds_nh[5], bounds_nh[6]]

        result = minimize(marginalized, x0, bounds=bounds, method='L-BFGS-B')
        contour_nh[i, j] = result.fun - min_chi2_global

# Same for IH
contour_ih = np.zeros((len(theta23_vals), len(dcp_vals)))
for i, t23 in enumerate(theta23_vals):
    for j, dcp in enumerate(dcp_vals):
        def marginalized(params):
            return chi2_objective_ih([params[0], t23, dcp, params[1], params[2], params[3], params[4]])

        x0 = [x0_ih[0], x0_ih[3], x0_ih[4], x0_ih[5], x0_ih[6]]
        bounds = [bounds_ih[0], bounds_ih[3], bounds_ih[4], bounds_ih[5], bounds_ih[6]]

        result = minimize(marginalized, x0, bounds=bounds, method='L-BFGS-B')
        contour_ih[i, j] = result.fun - min_chi2_global

# Convert theta23 values to sin^2(theta23)
sin2_theta23 = np.sin(theta23_vals)**2

# Plotting
X, Y = np.meshgrid(dcp_vals * 180/np.pi, sin2_theta23)
fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

cs1 = ax[0].contour(X, Y, contour_nh, levels=[2.30, 11.83], colors=['blue', 'blue'], linestyles=['-', '--'])
ax[0].clabel(cs1, fmt={2.30: '1σ', 11.83: '3σ'})
ax[0].set_title("NH")

cs2 = ax[1].contour(X, Y, contour_ih, levels=[2.30, 11.83], colors=['red', 'red'], linestyles=['-', '--'])
ax[1].clabel(cs2, fmt={2.30: '1σ', 11.83: '3σ'})
ax[1].set_title("IH")

for a in ax:
    a.set_xlabel("δₑₚ [degrees]")
    a.set_ylabel("sin²θ₂₃")
    a.grid(True)

plt.tight_layout()
plt.show()

# Cleanup
lib.free_globes()
