import numpy as np

data = np.genfromtxt("new_anni_events.dat")
ellipsoid = np.genfromtxt("WorkingEllipsoidConfigurations.dat")

np.save("new_anni_events.npy", data)
np.save("WorkingEllipsoidConfigurations.npy", ellipsoid)
