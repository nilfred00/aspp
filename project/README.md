# Project description
This project will take the code from https://github.com/KlaraHackenstrass/pi-pi_stacking/tree/main/analysis_pi-pi-stacking, which is a code for calculating pi-pi stacking of aromatic rings from molecular dynamics trajectories of lignin. The idea is to improve it by transforming it into a package and making it more available for general use as well as adding a jupyter-notebook example of how to use the package.

# Project outcome
Added a class for the aromatic rings, changed the code to a package and added numpy docstrings and ValueError:s to the functions. Also incorporated MDAnalysis to load and treat the trajectories rather than going in to the coordinate file and extracting all information. This makes the code applicable to other systems than the ones tested in the article. The example.ipynb notebook goes over an example, which should be enough for a new user to use the code.    

