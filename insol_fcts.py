
"""
Functions for the computation of insolation

Code is inspired from Laskar 2004 fortran code, which can be found here:
http://www.imcce.fr/Equipes/ASD/insola/earth/La2004/insolsub.f

Modified by Nicolas Brown (02/2015)

Modifications:

"""
import math as m


def true_avg(hl,e,pibar):
	"""
	Calculates average longitude(hlm) from true longitude(hl), eccentricity(e)
	and longitude of perihelion(pibar).
	
	IN:
		hl: true longitude
		e: eccentricity
		pibar: longitude of perihelion
		
	OUT:
		hlm: average longitude
	"""
	eV = hl-pibar
	return hl-2*e*m.sin(eV)+(3*e**2/4+e**4/8)*m.sin(2*eV)-(e**3/3+e**5/8)*m.sin(3*eV)+5*e**4/32*m.sin(4*eV)-3*e**5/40*m.sin(5*eV)


def avg_true(hlm,e,pibar):
	"""
	Calculates true longitude(hl) from average longitude(hlm), eccentricity(e)
	and longitude of perihelion(pibar).
	
	IN:
		hlm: average longitude
		e: eccentricity
		pibar: longitude of perihelion
		
	OUT:
		hl: true longitude
	"""
	eM = hlm-pibar
	return hlm+(2*e-e**3/4 + 5*e**5/96)*m.sin(eM)+(5*e**2/4 - 11*e**4/24)*m.sin(2*eM)+(13*e**3/12 - 43*e**5/64)*m.sin(3*eM)+ 103*e**4/96*m.sin(4*eM) + 1097*e**5/960*m.sin(5*eM)


def wmcal(month,e,eps,pibarh,phi,w):
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