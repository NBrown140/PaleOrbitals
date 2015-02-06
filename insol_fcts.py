
"""
Functions for the computation of insolation

Code is inspired from Laskar 2004 fortran code, which can be found here:
http://www.imcce.fr/Equipes/ASD/insola/earth/La2004/insolsub.f

Modified by Nicolas Brown (02/2015)

Modifications:

"""




def wmcal(month,e,eps,pibarh,phi,w)
	"""
	Monthly insolation for a given latitude 
	
	Dependencies:
		trueavg()
		
	
	IN:
		month: 
		e,eps,pibarh: orbital elements for date tps
			e: eccentricity
			eps: obliquity
			pibarh: longitude of perihelion from equinox at given date
		phi: latitude of point on Earth
		
	OUT:
		w: monthly insolation at given latitude
	"""
	
	
	pibar = pibarh + pi