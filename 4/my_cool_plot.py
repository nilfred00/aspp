import numpy as np
import matplotlib.pyplot as plt

# Functions
def gaussian(x, A, x0, sigma):
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2))

def lorentzian(x, A, x0, gamma):
    return A * (gamma**2 / ((x - x0)**2 + gamma**2))

def pseudovoigt(x, A, x0, sigma, gamma, eta):
    G = gaussian(x, A, x0, sigma)
    L = lorentzian(x, A, x0, gamma)
    return eta * L + (1 - eta) * G, G, L

# x range 
x = np.linspace(-10, 10, 500)

# Functions
params = [1.0, 0.0, 1.5, 1.0, 0.4]
pV, G, L = pseudovoigt(x, *params)

# Generate noisy data
noise_level = 0.05
y_noisy = pV + np.random.normal(0, noise_level, size=x.shape)

# Plotting
fig, ax = plt.subplots(facecolor="#dda084", figsize=(8,5))
ax.set_facecolor("#e0ebd1")

plt.scatter(x, y_noisy, s=10, alpha=0.5, label='Data')
plt.plot(x, pV, label='pseudo-Voigt', linewidth=2.5)
plt.plot(x, G, '--', label='Gaussian', linewidth=2)
plt.plot(x, L, '--', label='Lorentzian', linewidth=2)

plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
ax.set_title('My cool plot', color="#151815")
plt.xlim(-10,10)
plt.legend()
plt.tight_layout()

fig.savefig('my_cool_plot.pdf', format='pdf', bbox_inches='tight')

plt.show()


