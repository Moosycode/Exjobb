#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 16:45:53 2023

@author: frost
"""
#------------------------------------------------------------------------------
import os

#------------------------------------------------------------------------------
import numpy as np                       ## for handling of data
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
import json
#------------------------------------------------------------------------------
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.family'] = 'monospace'
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['k','r','g','b']) 
plt.rcParams['figure.autolayout'] = True
#------------------------------------------------------------------------------
#import ImportData as ID
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
class ASCIIfile():
#------------------------------------------------------------------------------
    def __init__( self, filename, lines_in_header=0, lines_in_footer=0, delim=None ,encoding=None ):
        self.header = list()
        self.footer = list()
        self.data   = list()
        self.name  = os.path.splitext(os.path.basename(filename))[0]
        self.filename = filename
        self.lines_in_header = lines_in_header
        self.lines_in_footer = lines_in_footer
    # read the file contents into a list
        datalist = list()
        if encoding: file = open(filename, 'r', encoding=encoding)
        else: file = open(filename, 'r')
        for line in file:
            if delim: temp = line.strip('\n').split(delim)
            else: temp = line.strip('\n').split()
            datalist.append(temp)
        file.close()
    # seperate the data
        if lines_in_header > 0: self.header = datalist[0:lines_in_header]
        if lines_in_footer > 0: self.footer = datalist[len(datalist)-(lines_in_footer+1):len(datalist)-1]
        for row in datalist[lines_in_header:len(datalist)-(lines_in_footer)]:
            try:
                self.data.append([float(s) for s in row])
            except:
                print(f'could not convert line {row} to data')
        self.data = np.array(self.data)


#------------------------------------------------------------------------------
# --- element dictionary -----------------------------------
#------------------------------------------------------------------------------
element = {
	'N':{'A':14,'Z':7},
	'Si':{'A':28,'Z':14},
	'Cl':{'A':35,'Z':17},
	'Ti':{'A':48,'Z':22},
	'Cu':{'A':63,'Z':29},
	'Br':{'A':79,'Z':35},
	'Mo':{'A':96,'Z':42},
	'Ag':{'A':107,'Z':47},
	'I':{'A':127,'Z':53},
	'Ta':{'A':181,'Z':73},
	'Au':{'A':197,'Z':79}
	}
#------------------------------------------------------------------------------


def correct_loss( x, c, loss ):
	corr = loss*x                                 # fraction lost for each depth
#	abs_corr = c/corr                             # absolute loss for each depth
#	c_corr = c+ c*loss#abs_corr
	c_corr = np.zeros(len(c))
	for i in range(len(c)):
		if c[i] > 0 and c[i]+corr[i]>0:	c_corr[i] = c[i]+corr[i]
		else: c_corr[i] = 0
	return c_corr
	# -----------------------------


def tscs( Ap, Zp, Ar, Zr, E ):
	cs = Ap * Zp**2 * Zr**2 / ( Ar * E )
	return cs


def tscs_2( Ap, Zp, Ar1, Zr1, Ar2, Zr2, E ):
	cs1 = Ap * Zp**2 * Zr1**2 / ( Ar1 * E )
	cs2 = Ap * Zp**2 * Zr2**2 / ( Ar2 * E )
	return (cs1+cs2)/2



#------------------------------------------------------------------------------
#--- generates the correction factor from system paramiter and empitcal coeffs
#------------------------------------------------------------------------------

def get_loss( Ap, Zp, Ar, Zr, E, co_A, co_B, verbose=False ):

	if verbose: print('\n# .... taking the log of the primary-ion mass number ....')
	log_Ap = np.log(Ap)
	if verbose: print(f'  log(Ap) = {log_Ap}')
	
	
	if verbose: print('\n# put into the empirical formula ....')
	log_factor = co_A*log_Ap + co_B
	if verbose: print(f'  ln[f($\sigma$,$\phi$)] = {log_factor}')
#	alt = co_a*np.log(Ap)-co_b
#	print(f'alt_fac = {alt}')
	
	if verbose: print('\n# ... take the exponent to get the slope of loss in counts with stopping cross-section ....')
	factor = np.exp(log_factor)
#	factor = np.exp(co_a*np.log(Ap)-co_b)
	if verbose: print(f'  f($\sigma$,$\phi$) = {factor}')
	
	if verbose: print('\n# ... calculate the stopping cross-section .... ')
	cross_sec = tscs( Ap, Zp, Ar, Zr, E )
	if verbose: print(f'  $\sigma$ = {cross_sec}')
	
	if verbose: print('\n# ... multiply ....')
	phi = factor * cross_sec
	print(f'  $\phi$= {phi}')

	return phi

#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#--- generates the correction factor from system paramiter and empitcal coeffs
#------------------------------------------------------------------------------

def get_loss_2( Ap, Zp, Ar1, Zr1, Ar2, Zr2, E, co_A, co_B, verbose=False ):

	if verbose: print('\n# .... taking the log of the primary-ion mass number ....')
	log_Ap = np.log(Ap)
	if verbose: print(f'  log(Ap) = {log_Ap}')
	
	
	if verbose: print('\n# put into the empirical formula ....')
	log_factor = co_A*log_Ap + co_B
	if verbose: print(f'  ln[f($\sigma$,$\phi$)] = {log_factor}')
#	alt = co_a*np.log(Ap)-co_b
#	print(f'alt_fac = {alt}')
	
	if verbose: print('\n# ... take the exponent to get the slope of loss in counts with stopping cross-section ....')
	factor = np.exp(log_factor)
#	factor = np.exp(co_a*np.log(Ap)-co_b)
	if verbose: print(f'  f($\sigma$,$\phi$) = {factor}')
	
	if verbose: print('\n# ... calculate the stopping cross-section .... ')
	cross_sec = tscs_2( Ap, Zp, Ar1, Zr1, Ar2, Zr2, E )
	if verbose: print(f'  $\sigma$ = {cross_sec}')
	
	if verbose: print('\n# ... multiply ....')
	phi = factor * cross_sec
	print(f'  $\phi$= {phi}')

	return phi

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
def read_profiles( folder_path, sample_keys, profile_keys, beam_keys, rebin=False, renorm=False, verbose=False, plot=False ):
	
	dp_list=[]
	if verbose: print('\n starting !!!!!!!!!')
	
	# --- determin if we are already in a request folder and 
	#	otherwise build a list of request folders --------------------
	if '.potku' in folder_path:
		print(' ... folder path is a Potku folder ...')
		files = [folder_path]
	else:
		print( '... searching folder path for Potku folders ...')
		files = [x for x in os.listdir(in_folder) if '.potku' in x]
		for i in range(len(files)):
			files[i]=folder_path+files[i]
			
	# ----------------------------------------------------------------
	
	# --- cycle through the list of request folders ------------------
	for file in files: 
		# --- read the 'measurement' file to get ion and energy -----
		mdat = json.load(open(file+'/Default/Default.measurement'))
		ion = mdat['beam']['ion']
		Sion = ''.join(filter(lambda x: x.isalpha(), ion))   # chemical symbol for ion
		energy = mdat['beam']['energy']
		print(f'found "{file}" with beam: {ion} at {energy} MeV')
		# -----------------------------------------------------------
		for x in filter(lambda bk: bk in ion+str(energy), beam_keys):           # check the ion-data for matching beam_keys
			print(' ... matched beam key ...')
			if verbose: print('\n'+file)
			# ---------------------------------------------------------
			subfiles = os.listdir(file+'/')
			for subfile in subfiles:
				for x in filter(lambda sk: sk in subfile, sample_keys):           # search sub-file names for the sample_keys
					if verbose: print('... '+subfile)
					measfiles = os.listdir(file+'/'+subfile+'/')
					for measfile in measfiles:
						if 'Measurement_' in measfile:
							dp_path = file+'/'+subfile+'/'+measfile+'/Depth_profiles/'
							dps = [x for x in os.listdir(dp_path) if 'depth.' in x and '.total' not in x]
							for dp in filter(lambda dp: dp in dps, profile_keys):           # search depth-profile names for the profile_keys
								if verbose: print('...... '+dp)
								
								# --- read the profile data ----
								pdata = ASCIIfile( dp_path+dp, lines_in_header=0, delim=None ).data 	# read in the profile data
								x = pdata[:,0]
								c = pdata[:,3]
								
								# --- assign label for the ion ----
								#ion_label = f'$^{element[Sion]["A"]}${Sion}'
								ion_zorder = element[Sion]["A"]
								#print(ion_label)
								
								# --- assign label for the profile -----
								Sdp = dp[6:]
								dp_zorder = element[Sdp]["A"]
								#print(dp_label)
								
								# --- take the (centered) derivative -----
								dc = np.zeros(len(c))
								for i in range(1,len(c)-1):
									dc[i]=(c[i+1]-c[i-1])/2.0;
								
								# --- rebin the data (double bin width) -----
								if rebin:
									c = c[1::2]+c[:-1:2]
									x = x[:-1:2]
								
								# --- renormailse the data to specific bin -----
								if renorm:
									#surf = c[sbz+1]					# massive hack, need defining properly
									#surf = c[0]					# massive hack, need defining properly
									surf = max(c)
									c /= surf
									
								dp_list.append([x, c, dc, energy, Sion, Sdp, ion_zorder, dp_zorder])

	if plot:
		plt.figure(figsize=[4,3])
		for dpt in dp_list: 
			plt.step(dpt[1],dpt[2], label=ion_label)
		plt.xlim([-1000,   5000])
		plt.ylim([    0.0, 1.2])
		plt.xlabel('depth [10$^{15}$ atoms/cm$^2$]')
		plt.ylabel('atomic fraction')
		plt.grid(linestyle='--')
		plt.legend()

	
	return dp_list




def simple_func( x, m ):
	y = m*x +1
	return y

#------------------------------------------------------------------------------
def lin_func( x, m, c ):
	y = m*x + c
	return y
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def quad_func( x, a, b, c ):
	y = a*x**2 + b*x + c
	return y
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def fit_profile_slop( start, stop, x, c, func ):
	
	x_ = []
	c_ = []
	for i in range(len(x)):
		if c[i] > 0:
			x_.append(x[i])
			c_.append(c[i])
			if x_[-1] < start:
				x_start = len(x_)-1
			if x_[-1] < stop:
				x_stop = len(x_)-1
	
	print(x_[x_start:x_stop])
	print(c_[x_start:x_stop])
	
	coeff, Cov = curve_fit( func, x_[x_start:x_stop], c_[x_start:x_stop] )
	sdev = []
	for i in range(len(Cov)): sdev.append(Cov[i,i]**(1/2))
	for co,sd in zip(coeff,sdev): print(f'   {np.round(co,5)}\t +/- {np.round(sd,5)}')
	
	return(coeff)
#------------------------------------------------------------------------------





#####################################################################################


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

in_folder = '/Users/niwi9751/potku/requests/'

SAVE_FIGURE = False
PLOT_FIGURE = False

start = 100
stop = 3000
#func = lin_func

#func = simple_func


save_path = 'pyplots2/'


print('--- reading in deapth profiles ....')
print('\n')





#####################################################################################




#%%#################################################################################
###   Load in the data and fit the slope of the fluece with depth for each   #####
####################################################################################

#####################################################################################
# LOOK AT EVERYTHING !!!!!!!!!!!

sample_key  = ['Sample_02','Sample_07','Sample_08']
profile_key = ['depth.Ta','depth.Mo','depth.Cu']        # use 'depth.' for all elements
ion_key     = ['24','32','36']              # this can be energy (e.g. 36), element (e.g. I) or both (e.g. I36)
rmass = 180

dpt_list = read_profiles(in_folder, sample_key, profile_key, ion_key, rebin=False, renorm=False, verbose=False, plot=False )

### x, c, dc, energy, ion_label, dp_label, ion_zorder, dp_zorder ###

#%%

title = 'Ta measured with different primary-ions at 24 MeV'
if SAVE_FIGURE: 
	plt.savefig(save_path+title, dpi=500)
elif PLOT_FIGURE:
	plt.title(title)

rZ_list = []
rmass_list = []
mass_list = []
eng_list = []
slope_list = []

plt.figure()
for fit_data in dpt_list:
	
	x = fit_data[0]
	c = fit_data[1]
	energy = fit_data[3]
	ion_label = fit_data[4]
	dp_label = fit_data[5]
	#print(fit_data)
	
	print('\n')
#	print(f'eng = {fit_data[4]}, mass = {fit_data[5]}, rmass = {rmass}')
#	print(f'mass/(rmass*eng) = {fit_data[5]/(rmass*fit_data[4])}')
	coeff = fit_profile_slop( start, stop, x, c, simple_func )
	plt.plot( x, c, 'o', linestyle=None )
	plt.plot( x, simple_func(x,coeff[0]) )
	
	rZ_list.append(element[dp_label]["Z"])
	rmass_list.append(element[dp_label]["A"])
	mass_list.append(element[ion_label]["A"])
	eng_list.append(energy)
	slope_list.append(coeff[0])
	
plt.xlim([-1000,   5000])
plt.ylim([    0.0, 1.5])
plt.xlabel('depth [10$^{15}$ atoms/cm$^2$]')
plt.ylabel('atomic fraction')
plt.grid(linestyle='--')
plt.title(title)

# --- convert everything to numpy arrays -----
rZ_list = np.array(rZ_list)
rmass_list = np.array(rmass_list)
mass_list = np.array(mass_list)
eng_list = np.array(eng_list)
slope_list = np.array(slope_list)
# --------------------------------------------

#fit_data = dpt_list[0]
#fit_profile_slop( start, stop, fit_data[1], fit_data[2], func )

#####

plt.figure()
for dpt in dpt_list: plt.plot(0,0,label=dpt[3]) 
plt.legend()
title = 'X measured with different primary-ions at 24 MeV - legend'
if SAVE_FIGURE: 
	plt.savefig(save_path+title, dpi=500)
elif PLOT_FIGURE:
	plt.title(title)
	
#####################################################################################

print('---------------------------------')


#%%#################################################################################
###   Plot the gradients of the fits vs the total-stopping cross-section   #####
####################################################################################


Au_list = []
I_list = []
Ag_list = []
Br_list = []
Cl_list = []

plt.figure()
for rz,rmass,mass,eng,slope in zip(rZ_list,rmass_list,mass_list,eng_list,slope_list):
	
	ls=None
	m=''
	c='c'
	
	if rmass==181: 
		m='o'
		
	if rmass==96: 
		m='s'

	if rmass==63: 
		m= '*'

		
	if mass==197: 
		c='k'
		z=79
		Au_list.append([rz, rmass, z, mass, eng, slope])
	if mass==127: 
		c='r'
		z=53
		I_list.append([rz, rmass, z, mass, eng, slope])
	if mass==107: 
		c= 'g'
		z=47
		Ag_list.append([rz, rmass, z, mass, eng, slope])
	if mass==79: 
		c= 'b'
		z=35
		Br_list.append([rz, rmass, z, mass, eng, slope])
	if mass==35: 
		c= 'm'
		z=17
		Cl_list.append([rz, rmass, z, mass, eng, slope])
		
	if eng== 24:
		em='v'
	if eng== 32:
		em='<'
	if eng==36:
		em='>'
	if eng==48:
		em='^'
		
	
	x =   mass * z**2 * rz**2 / ( rmass * eng ) 
	#x =  (1 / ((rmass+mass)**2)) * z**2 * rz**2 / ( eng ) 
	
	plt.plot(x,-slope, linestyle=None, marker=m, markersize=10, c=c)
	plt.plot(x,-slope, linestyle=None, marker=em, markersize=5, c='c')

plt.xlabel('[ A$_p$$\cdot$Z$_p$$^2$$\cdot$Z$_r$$^2$ / (A$_r$$\cdot$E ) ]')
plt.ylabel('[ -1 / 10$^{15}$ atoms / cm$^2$ ]')
plt.grid(linestyle='--')

plt.legend(['shapes sperate target masses','colours seperate ion masses'])


x = eng * (rmass*mass) / (rmass+mass)**2 

y = mass * z**2 * rz**2 / (rmass * eng)





#%%#################################################################################
###   Perfrom a linear fit to the distributions for each primary ion used     #####
###      (very messy! - needs cleaning up)
####################################################################################

plt.figure()
x_ = np.linspace(-0.1e6,2e6,100)

def func2(x, m):
	y = x*m
	return y



rz_l = []
rmass_l  = []
z_l  = []
mass_l = []
eng_l = []
slope_l = []

for i in Au_list:
	
	rz = i[0]
	rmass  = i[1]
	z  = i[2]
	mass = i[3]
	eng = i[4]
	slope = i[5]
	
	#x =  (1 / ((rmass+mass)**2)) * z**2 * rz**2 / ( eng ) 
	x =   mass * z**2 * rz**2 / ( rmass * eng ) 
	y = -slope
	
	print(i)

	plt.plot(x,y,'ko')
	
	rz_l.append(i[0])
	rmass_l.append(i[1])
	z_l.append(i[2])
	mass_l.append(i[3])
	eng_l.append(i[4])
	slope_l.append(-i[5])

rz_l = np.array(rz_l)
rmass_l  = np.array(rmass_l)
z_l  = np.array(z_l)
mass_l = np.array(mass_l)
eng_l = np.array(eng_l)
slope_l = np.array(slope_l)


x_l =   mass_l * z_l**2 * rz_l**2 / ( rmass_l * eng_l )

coeff, Cov = curve_fit( func2, x_l, slope_l )

sdev = []
for i in range(len(Cov)): sdev.append(Cov[i,i]**(1/2))
for co,sd in zip(coeff,sdev): print(f'   {np.round(co,5)}\t +/- {np.round(sd,5)}')



plt.plot( x_, func2(x_,coeff[0]), 'k', linestyle='--' )

Au_factor = coeff[0]

#######

rz_l = []
rmass_l  = []
z_l  = []
mass_l = []
eng_l = []
slope_l = []

for i in I_list:
	
	rz = i[0]
	rmass  = i[1]
	z  = i[2]
	mass = i[3]
	eng = i[4]
	slope = i[5]
	
	#x =  (1 / ((rmass+mass)**2)) * z**2 * rz**2 / ( eng ) 
	x =   mass * z**2 * rz**2 / ( rmass * eng ) 
	y = -slope
	
	print(i)

	plt.plot(x,y,'ro')
	
	rz_l.append(i[0])
	rmass_l.append(i[1])
	z_l.append(i[2])
	mass_l.append(i[3])
	eng_l.append(i[4])
	slope_l.append(-i[5])

rz_l = np.array(rz_l)
rmass_l  = np.array(rmass_l)
z_l  = np.array(z_l)
mass_l = np.array(mass_l)
eng_l = np.array(eng_l)
slope_l = np.array(slope_l)


x_l =   mass_l * z_l**2 * rz_l**2 / ( rmass_l * eng_l )

coeff, Cov = curve_fit( func2, x_l, slope_l )

sdev = []
for i in range(len(Cov)): sdev.append(Cov[i,i]**(1/2))
for co,sd in zip(coeff,sdev): print(f'   {np.round(co,5)}\t +/- {np.round(sd,5)}')



plt.plot( x_, func2(x_,coeff[0]), 'r', linestyle='--' )

I_factor = coeff[0]

#######

rz_l = []
rmass_l  = []
z_l  = []
mass_l = []
eng_l = []
slope_l = []

for i in Ag_list:
	
	rz = i[0]
	rmass  = i[1]
	z  = i[2]
	mass = i[3]
	eng = i[4]
	slope = i[5]
	
	#x =  (1 / ((rmass+mass)**2)) * z**2 * rz**2 / ( eng ) 
	x =   mass * z**2 * rz**2 / ( rmass * eng ) 
	y = -slope
	
	print(i)

	plt.plot(x,y,'go')
	
	rz_l.append(i[0])
	rmass_l.append(i[1])
	z_l.append(i[2])
	mass_l.append(i[3])
	eng_l.append(i[4])
	slope_l.append(-i[5])

rz_l = np.array(rz_l)
rmass_l  = np.array(rmass_l)
z_l  = np.array(z_l)
mass_l = np.array(mass_l)
eng_l = np.array(eng_l)
slope_l = np.array(slope_l)


x_l =   mass_l * z_l**2 * rz_l**2 / ( rmass_l * eng_l )

coeff, Cov = curve_fit( func2, x_l, slope_l )

sdev = []
for i in range(len(Cov)): sdev.append(Cov[i,i]**(1/2))
for co,sd in zip(coeff,sdev): print(f'   {np.round(co,5)}\t +/- {np.round(sd,5)}')



plt.plot( x_, func2(x_,coeff[0]), 'g', linestyle='--' )

Ag_factor = coeff[0]

#######

rz_l = []
rmass_l  = []
z_l  = []
mass_l = []
eng_l = []
slope_l = []

for i in Br_list:
	
	rz = i[0]
	rmass  = i[1]
	z  = i[2]
	mass = i[3]
	eng = i[4]
	slope = i[5]
	
	#x =  (1 / ((rmass+mass)**2)) * z**2 * rz**2 / ( eng ) 
	x =   mass * z**2 * rz**2 / ( rmass * eng ) 
	y = -slope
	
	print(i)

	plt.plot(x,y,'bo')
	
	rz_l.append(i[0])
	rmass_l.append(i[1])
	z_l.append(i[2])
	mass_l.append(i[3])
	eng_l.append(i[4])
	slope_l.append(-i[5])

rz_l = np.array(rz_l)
rmass_l  = np.array(rmass_l)
z_l  = np.array(z_l)
mass_l = np.array(mass_l)
eng_l = np.array(eng_l)
slope_l = np.array(slope_l)


x_l =   mass_l * z_l**2 * rz_l**2 / ( rmass_l * eng_l )

coeff, Cov = curve_fit( func2, x_l, slope_l )

sdev = []
for i in range(len(Cov)): sdev.append(Cov[i,i]**(1/2))
for co,sd in zip(coeff,sdev): print(f'   {np.round(co,5)}\t +/- {np.round(sd,5)}')



plt.plot( x_, func2(x_,coeff[0]), 'b', linestyle='--' )

Br_factor = coeff[0]

#######

rz_l = []
rmass_l  = []
z_l  = []
mass_l = []
eng_l = []
slope_l = []

for i in Cl_list:
	
	rz = i[0]
	rmass  = i[1]
	z  = i[2]
	mass = i[3]
	eng = i[4]
	slope = i[5]
	
	#x =  (1 / ((rmass+mass)**2)) * z**2 * rz**2 / ( eng ) 
	x =   mass * z**2 * rz**2 / ( rmass * eng ) 
	y = -slope
	
	print(i)

	plt.plot(x,y,'mo')
	
	rz_l.append(i[0])
	rmass_l.append(i[1])
	z_l.append(i[2])
	mass_l.append(i[3])
	eng_l.append(i[4])
	slope_l.append(-i[5])

rz_l = np.array(rz_l)
rmass_l  = np.array(rmass_l)
z_l  = np.array(z_l)
mass_l = np.array(mass_l)
eng_l = np.array(eng_l)
slope_l = np.array(slope_l)


x_l =   mass_l * z_l**2 * rz_l**2 / ( rmass_l * eng_l )

coeff, Cov = curve_fit( func2, x_l, slope_l )

sdev = []
for i in range(len(Cov)): sdev.append(Cov[i,i]**(1/2))
for co,sd in zip(coeff,sdev): print(f'   {np.round(co,5)}\t +/- {np.round(sd,5)}')



plt.plot( x_, func2(x_,coeff[0]), 'm', linestyle='--' )

Cl_factor = coeff[0]







plt.xlabel(' A$_p$$\cdot$Z$_p$$^2$$\cdot$Z$_r$$^2$ / (A$_r$$\cdot$E ) ')
plt.ylabel(' -1 / 10$^{15}$ atoms / cm$^2$ ')


plt.xlim([-0.1e6,2.1e6])
plt.ylim([-0.0002,0.003])

plt.grid(linestyle='--')





#%%#################################################################################
###   Plot the gradients of the fits vs primary-ion mass on ln-ln and fit    #####
####################################################################################


log_masses  = np.log(np.array([197,127,107,79,35]))
log_factors = np.log(np.array([Au_factor,I_factor,Ag_factor,Br_factor,Cl_factor,]))


print(log_factors)
print(log_masses)


plt.figure()
for l_mass,l_fact in zip(log_masses,log_factors):
	plt.plot(l_mass,l_fact, '-o')



coeff, Cov = curve_fit( lin_func, log_masses, log_factors )

sdev = []
for i in range(len(Cov)): sdev.append(Cov[i,i]**(1/2))
print('\n')
print('###################################################')
for co,sd in zip(coeff,sdev): print(f'   {np.round(co,5)}\t +/- {np.round(sd,5)}')
print('###################################################')
print('\n')

x_ = np.linspace(3,6,3)

plt.plot( x_, lin_func(x_,coeff[0],coeff[1]), 'b', linestyle='--' )
plt.xlabel(' ln[A$_p$] ')
plt.ylabel(' ln[f($\sigma$,$\phi$)]')
plt.grid(linestyle='--')




#####################################################################
#####################################################################
# THE CORRECTION COEFFICENTS ARE: co[0] and co[1]
#####################################################################
#####################################################################



"""
#%%#################################################################################
###       #####
####################################################################################

# 1st order correction for any ion-energy-target combination is performed by:

# .... getting the sytem configuration for ....


# .... the beam: 44 MeV 127-I on Xe ....
Ap = 127
Zp = 53
E  = 44 

# .... then target element: Xe ....
Ar_Xe = 131
Zr_Xe = 54

# .... then target element: Kr ....
Ar_Kr = 84
Zr_Kr = 36

# .... then target element: U ....
Ar_U = 238
Zr_U = 92

# .... then target element: O ....
Ar_O = 16
Zr_O = 8


co_a = co[0]#-3.22203
co_b = co[1]#-4.59366



"""






#%%#################################################################################
###   redload the original data and apply correction    #####
####################################################################################


#sample_key  = ['Sample_02']
#profile_key = ['depth.Cu']        # use 'depth.' for all elements
#ion_key     = ['36']              # this can be energy (e.g. 36), element (e.g. I) or both (e.g. I36)
#E = 36






sample_key  = ['Sample_02','Sample_07','Sample_08']
profile_key = ['depth.Ta','depth.Mo','depth.Cu']        # use 'depth.' for all elements
ion_key     = ['24','32','36','48']              # this can be energy (e.g. 36), element (e.g. I) or both (e.g. I36)


for profile in profile_key:   # load one sample at a time ...

	dpt_list = read_profiles( in_folder, sample_key, [profile], ion_key, rebin=True, renorm=True, verbose=True, plot=False )
	
	# x, c, dc, energy, Sion, Sdp, ion_zorder, dp_zorder
	
	for d in dpt_list:
		
		# --- get the beam paramiters --------------------------------------
		beam_elem = d[4]    # Sdp
		Ap = element[beam_elem]['A']
		Zp = element[beam_elem]['Z']
		E  = d[3]
		print(f'beam element is: {beam_elem}, A = {Ap}, Z = {Zp} @ {E} MeV')
		# ------------------------------------------------------------------
		
		# --- get the target paramiters ------------------------------------
		target_elem = d[5]    # Sdp
		Ar = element[target_elem]['A']
		Zr = element[target_elem]['Z']
		print(f'target element is: {target_elem}, A = {Ar}, Z = {Zr} ')
		# ------------------------------------------------------------------
		
		# ------------------------------------------------------------------
		x = d[0]
		c = d[1]
		
		loss = get_loss( Ap, Zp, Ar, Zr, E, coeff[0], coeff[1] )
		print(f'\n---------------\nloss = {loss}\n---------------------\n')
		dc = np.zeros(len(c))
		for i in range(1,len(c)-1):
			dc[i-1] = c[i-1]-c[i]
			
		corr = loss*x
		c_corr = np.zeros(len(c))
		for i in range(len(c)):
			if c[i] > 0 and c[i]+corr[i]>0:	c_corr[i] = c[i]+corr[i]
			else: c_corr[i] = 0
		# ------------------------------------------------------------------
		
		# ------------------------------------------------------------------
		div_data = abs(1-c)
		div_corr = abs(1-c_corr)
		div_div  = (div_data/div_corr)
		# ------------------------------------------------------------------
		
		
		
		
		
		
		
		"""
		plt.figure()
		
		plt.plot( x, c, 'o', linestyle=None, label=target_elem )
		plt.plot(x,corr,'--',label='corr_fact')
		plt.plot(x,c_corr,'--',label='corrected')
		
		plt.xlim([-1000,   5000])
		plt.ylim([    -0.2, 1.5])
		plt.xlabel('depth [10$^{15}$ atoms/cm$^2$]')
		plt.ylabel('atomic fraction')
		plt.grid(linestyle='--')
		plt.legend()
		plt.title(f'{element[beam_elem]["A"]}-{beam_elem} @ {E} MeV')
		"""
		

		
		fig, ax = plt.subplots(ncols=2, figsize=(9, 3), layout='constrained')

		ax[0].plot( x, c, 'o', linestyle=None, label=target_elem )
		ax[0].plot(x,c_corr,'--',label='corrected')
		ax[0].set_xlim([-1000,   5000])
		ax[0].set_ylim([    -0.1, 1.4])
		ax[0].set_xlabel('depth [10$^{15}$ atoms/cm$^2$]')
		ax[0].set_ylabel('atomic fraction')
		ax[0].grid(linestyle='--')
		ax[0].legend()
		
		ax[1].plot(x,div_div, '-', label='div_div')
		ax[1].plot(x,div_data,'--',label='div_data')
		ax[1].plot(x,div_corr,'--',label='div_corr')
		ax[1].set_xlim([-1000,   5000])
#		ax[1].set_ylim([    -0.75, 0.75])
		ax[1].set_xlabel('depth [10$^{15}$ atoms/cm$^2$]')
		ax[1].set_ylabel('atomic fraction')
		ax[1].grid(linestyle='--')
		ax[1].legend()
		
		plt.suptitle(f'{element[beam_elem]["A"]}-{beam_elem} @ {E} MeV')
		
		
		
		### TODO: plot the differences for 1 to show improvement in the profiles
		
		### TODO apply to another system ... prefereable duel element !!!!
		
		


#%%##########################################################################################



sample_key  = ['Sample_05']
profile_key = ['depth.Ti','depth.N']        # use 'depth.' for all elements
ion_key     = ['24','32','36']              # this can be energy (e.g. 36), element (e.g. I) or both (e.g. I36)

folders = [x for x in os.listdir(in_folder) if '.potku' in x]

for folder in folders:

	dpt_list = read_profiles( folder, sample_key, profile_key, ion_key, rebin=False, renorm=False, verbose=True, plot=False )
	
	
	labels = []
	Ar = []
	Zr = []
	c  = []

	for d in dpt_list:
		
		# --- get the beam paramiters --------------------------------------
		beam_elem = d[4]    # Sdp
		Ap = element[beam_elem]['A']
		Zp = element[beam_elem]['Z']
		E  = d[3]
		print(f'beam element is: {beam_elem}, A = {Ap}, Z = {Zp} @ {E} MeV')
		# ------------------------------------------------------------------
		
		# --- get the target paramiters ------------------------------------
		target_elem = d[5]    # Sdp
		Ar = element[target_elem]['A']
		Zr = element[target_elem]['Z']
		print(f'target element is: {target_elem}, A = {Ar}, Z = {Zr} ')
		# ------------------------------------------------------------------
		
		# ------------------------------------------------------------------
		x = d[0]
		c = d[1]
		# ------------------------------------------------------------------

		
		if target_elem == 'Ti':
			Ar_Ti = Ar
			Zr_Ti = Zr
			lab_Ti = target_elem
			c_Ti = c
		
		if target_elem == 'N':
			Ar_N = Ar
			Zr_N = Zr
			lab_N = target_elem
			c_N = c
		# ------------------------------------------------------------------
		
	
	
	loss = get_loss_2( Ap, Zp, Ar_Ti, Zr_Ti, Ar_N, Zr_N, E, coeff[0], coeff[1] )
	print(f'\n---------------\nloss = {loss}\n---------------------\n')
	
	
	
	c_corr_Ti = correct_loss(x, c_Ti, loss)
	c_corr_N  = correct_loss(x, c_N, loss)
	
	
	
	
	fig, ax = plt.subplots(ncols=2, figsize=(9, 3), layout='constrained')
	
	ax[0].plot( x, c_Ti,      'o',  linestyle=None, label=lab_Ti )
	ax[0].plot( x, c_corr_Ti, '--',                 label='Ti corrected' )
	ax[0].plot( x, c_N,       'o',  linestyle=None, label=lab_N )
	ax[0].plot( x, c_corr_N,  '--',                 label='N corrected' )
		
	ax[0].set_xlim([-1000,   5000])
	ax[0].set_ylim([    -0.2, 1.5])
	ax[0].set_xlabel('depth [10$^{15}$ atoms/cm$^2$]')
	ax[0].set_ylabel('atomic fraction')
	ax[0].grid(linestyle='--')
	ax[0].legend()
	

	ax[1].plot( x, c_Ti/c_N,      'o',  linestyle=None, label='Ti:N ratio' )
	ax[1].plot( x, c_corr_Ti/c_corr_N, '--',                 label='corrected ratio' )
		
	ax[1].set_xlim([-1000,   5000])
	ax[1].set_ylim([    0.75, 1.25])
	ax[1].set_xlabel('depth [10$^{15}$ atoms/cm$^2$]')
	ax[1].set_ylabel('atomic fraction')
	ax[1].grid(linestyle='--')
	ax[1].legend()

	plt.suptitle(f'{element[beam_elem]["A"]}-{beam_elem} @ {E} MeV')





#%%####################################################################################

#folders = os.listdir(in_folder)
#for folder in folders:
#	print(folder)





















