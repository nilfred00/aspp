import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

# Load data
exp_data = np.load("I_q_IPA_exp.npy")
q_exp = exp_data[:,0]
I_exp = exp_data[:,1]

model_data = np.load("I_q_IPA_model.npy")
q_model = model_data[:,0]
I_model = model_data[:,1]

# Interpolate experimental data
interp_model = interp1d(q_exp, I_exp)

# Take q-values of the model in the experimental q-range for the scaling
i_scale = (q_model<=max(q_exp)) & (q_model>=min(q_exp)) & (q_model>0)
q_scale = q_model[i_scale]

# Interpolated intensity at the scaling q-values for the scaling
I_interp = interp_model(q_scale)

# Model data to scale to the interpolated intensity
I_scale = I_model[i_scale]

# Function to minimize
def fun(scale):
    return np.sum((I_interp - scale * I_scale)**2)

# Find the best scaling factor
res = minimize_scalar(fun)
scale = res.x
print(scale)

I_model_scaled = scale*I_model

# Plot result
plt.figure()
plt.plot(q_exp, I_exp, '*', label='Experimental data')
plt.plot(q_scale, I_interp, '*', label="Interpolation data")
plt.plot(q_model, I_model_scaled, label="Scaled model data")
plt.legend()
plt.xlabel(r'$q$ (Å$^{-1}$)')
plt.ylabel("Intensity (a.u.)")
plt.show()