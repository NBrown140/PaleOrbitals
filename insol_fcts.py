
"""
Functions for the computation of insolation

Code is inspired from Laskar 2004 fortran code, which can be found here:
http://www.imcce.fr/Equipes/ASD/insola/earth/La2004/insolsub.f

Kept same function(subroutine in fortran) names except for vraimoy and moyvrai,
which became true_avg and avg_true, respectively.

Modified by Nicolas Brown (02/2015)

Modifications:

"""
import math as m
import scipy.integrate.romberg



def true_avg(hl,e,pibar):
	"""
	Calculates average longitude(hlm) from true longitude(hl), eccentricity(e)
	and longitude of perihelion(pibar).
	
	IN:
		hl: true longitude....................................... float or int
		e: eccentricity.......................................... float or int
		pibar: longitude of perihelion........................... float or int
		
	OUT:
		hlm: average longitude................................... float
	"""
	
	# Convert input to floats if not already done
	hl,e,pibar = float(hl),float(e),float(pibar)
	
	eV = hl-pibar
	return hl-2*e*m.sin(eV)+(3*e**2/4+e**4/8)*m.sin(2*eV)-(e**3/3+e**5/8)*m.sin(3*eV) \
	+5*e**4/32*m.sin(4*eV)-3*e**5/40*m.sin(5*eV)


def avg_true(hlm,e,pibar):
	"""
	Calculates true longitude(hl) from average longitude(hlm), eccentricity(e)
	and longitude of perihelion(pibar).
	
	IN:
		hlm: average longitude................................... float or int
		e: eccentricity.......................................... float or int
		pibar: longitude of perihelion........................... float or int
		
	OUT:
		hl: true longitude....................................... float
	"""
	
	# Convert input to floats if not already done
	hlm,e,pibar = float(hlm),float(e),float(pibar)
	
	eM = hlm-pibar
	return hlm+(2*e-e**3/4 + 5*e**5/96)*m.sin(eM)+(5*e**2/4 - 11*e**4/24)*m.sin(2*eM) \
	+(13*e**3/12 - 43*e**5/64)*m.sin(3*eM)+ 103*e**4/96*m.sin(4*eM) \
	+1097*e**5/960*m.sin(5*eM)



def cwj(S0,wd,e,pibar,eps,phi):
	"""
	Daily insolation for a given latitude
	
	IN:
		S0: solar constant....................................... float or int
		wd: true longitude of sun (in radians) from
		    equinox at date...................................... float or int
		e: eccentricity.......................................... float or int
		pibar: longitude of perihelion from equinox
		       at date + pi...................................... float or int
		eps: obliquity........................................... float or int
		phi: latitude on Earth................................... float or int
		
	OUT:
		w: daily insolation...................................... float
	"""
	
	# Convert input to floats if not already done
	S0,wd,e,pibar,eps,phi = float(S0),float(wd),float(e),float(pibar),float(eps),float(phi)
	
	# True anomaly
	v = wd - pibar
	
	# Earth-Sun distance
	rho = (1-e**2) / (1+e*m.cos(v))
	
	# Sun declination
	delta = m.asin(m.sin(eps)*m.sin(wd))
	
	# Latitude can either have: (1) sunrise and sunset, (2) no sunset or
	# (3) no sunrise  
	aux = m.pi/2. - abs(delta)
	a1 = m.pi/2. - delta
	a2 = m.pi/2. + delta
	
	## (1) Sunrise and sunset
	if (-1*aux < phi) and (phi < aux):
		# Hour angle of sunrise and sunset
		ho = m.acos(-1*m.tan(phi) * m.tan(delta))
		# Insolation
		return (ho*m.sin(phi)*m.sin(delta) + m.cos(phi)*m.cos(delta)*m.sin(ho)) \
		* S0/(m.pi*rho**2)
		
	## (2) No sunset
	elif (phi > a1) or (phi < -1*a2):
		return S0 * m.sin(phi)*m.sin(delta) / rho**2
		
	## (3) No sunrise
	elif (phi < -1*a1) or (phi > a2):
		return 0.
	
	print 'Should never get here! if-statements should contain all possibilities'



def wjour(S0,date,e,eps,pibarh,phi):
	"""
	Daily insolation for a given latitude
	
	IN:
		S0: solar constant....................................... float or int
		date: true longitude of sun (in degrees) from
		      true equinox at date............................... float or int
		e: eccentricity.......................................... float or int
		pibarh: longitude of perihelion from equinox
		        at date.......................................... float or int
		eps: obliquity........................................... float or int
		phi: latitude on Earth................................... float or int
	OUT:
		w: daily insolation...................................... float
	"""
	
	# Convert input to floats if not already done
	S0,date,e,pibarh,eps,phi = float(S0),float(date),float(e),float(pibarh),float(eps),float(phi)

	# Convert pibarh to longitude of perihelion from equinox at date + pi
	pibar = pibarh + m.pi
	# Convert date to radians
	hl = m.radians(date)
	
	return cwj(S0,hl,e,pibar,eps,phi)



def wjcal(S0,datecal,e,eps,pibarh,phi):
	"""
	Daily insolation for a given latitude
	
	IN:
		S0: solar constant....................................... float or int
		datecal: avg longitude of sun (in degrees) from
		         equinox at date................................. float or int
		e: eccentricity.......................................... float or int
		pibarh: longitude of perihelion from equinox
		        at date.......................................... float or int
		eps: obliquity........................................... float or int
		phi: latitude on Earth................................... float or int
	OUT:
		w: daily insolation...................................... float
	"""
	
	# Convert input to floats if not already done
	S0,datecal,e,pibarh,eps,phi = float(S0),float(datecal),float(e),float(pibarh),float(eps),float(phi)

	# Convert pibarh to longitude of perihelion from equinox at date + pi
	pibar = pibarh + m.pi	
	# Convert datecal to radians, obtain true longitude
	hlm = m.radians(datecal)
	# Get true longitude on March 21st
	hlm0 = true_avg(0.,e,pibar)
	# Avg longitude at date
	hlm = hlm0 + hlm
	# True longitude at date
	wd = avg_true(hlm,e,pibar,)
	
	return cwj(S0,wd,e,pibar,eps,phi)



def wam(S0,e):
	""""
	Annual global average insolation	
	
	IN:
		e: eccentricity.......................................... float or int
		
	OUT:
		w: annual global average insolation...................... float
	"""
	
	# Convert input to floats if not already done
	S0,e = float(S0),float(e)
	
	return S0/(4.*m.sqrt(1.-e**2))
	
	

def wmcal(month,e,eps,pibarh,phi):
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
	
	# Convert pibarh to longitude of perihelion from equinox at date + pi
	pibar = pibarh + m.pi	
	# Average longitude on March 21st
	hlm0 = true_avg(0.,2,pibar)
	# Average longitude at beginning of month
	hlm1 = hlm0 + (month-4)*m.pi*30./180
	# Average longitude at end of month
	hlm2 = hlm1 + 30.*m.pi/180
	
	# Define function F to integrate
	def F(hlm):
		hl = avg_true(hlm,e,pibar)
		return cwj(hl,e,pibar,eps,phi)
		
	# Calculate insolation with Romberg Integration
	w = scipy.integrate.romberg(F,hlm1,hlm2)
	w = w/30./m.pi*180.
	return w
	
	