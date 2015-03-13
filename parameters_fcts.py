
"""

"""

import numpy as np

laskar04 = np.load('datasets/laskar04.npy')
laskar04_age   = laskar04[:,0]
laskar04_ecc = laskar04[:,1]
laskar04_obl = laskar04[:,2]
laskar04_pre = laskar04[:,3]

def get_orbs(age,dataset):
	"""
	Returns orbital parameters at given age (linearly interpolates if necessary)
	
	IN:
		age: time from present in kyrs, past >0....................... float or int
		dataset: dataset from which orbital parameters
				 will be obtained. Options: 'laskar04'................ string
	OUT:
		(e,obliq,prec): 3 orbital parameters.......................... tuple
			e: eccentricity........................................... float
			obliq: obliquity.......................................... float
			prec: precession.......................................... float
	"""
	if dataset=='laskar04':
		assert(age>0 and age<5000)
		
		
		
		
		return e,obliq,prec