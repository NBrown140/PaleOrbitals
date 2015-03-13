
import math as m
import numpy as np
import matplotlib.pyplot as plt

import insol_fcts
import parameters_fcts


def plot_lat_insol(S0,month,e,obliq,prec):
	"""
	orbital parameters in degrees 
	
	
	OUT:
		(matplotlib figure, insol)................................. tuple
			matplotlib figure...................................... matplotlib figure
			insol: insolation array................................ 1-D array
	"""
	
	lat = np.linspace(-90,90,181)
	insol = np.empty(len(lat)) * np.nan
	
	for i in range(len(lat)):
		insol[i] = insol_fcts.wmcal(S0,month,e,m.radians(obliq),m.radians(prec),m.radians(lat[i]))
	
	plt.figure(1)
	plt.xlim(-90, 90)
	plt.ylim(0, 600)
	plt.xlabel('Latitude (deg)')
	plt.ylabel('Insolation (W m-2)')
	plt.title('Average Monthly Insolation')
	plt.grid(True)
	
	plot = plt.plot(lat,insol)
	plt.show()
	
	return plot,insol
	



