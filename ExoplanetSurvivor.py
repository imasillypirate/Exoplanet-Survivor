import os 
import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
import random as rnd
import Tkinter as Tk
import tkMessageBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backend_bases import key_press_handler

def intFromList(ifl_x,ifl_y):
#trapezoid integration of lists 
#inputs: ifl_x = independent variable (list)
#        ifl_y = dependent variable (same sized list)
#Output: ifl_out = trapezoid integration of y dx over list
	ifl_out = 0.0
	for ifl_i in range(len(ifl_x)-1):
		ifl_out = ifl_out + 1.0*(ifl_x[ifl_i+1]-ifl_x[ifl_i])*(ifl_y[ifl_i+1] + ifl_y[ifl_i])/2
		
	return ifl_out

def AppMag(AM_Dist, AM_L):
	#Determine the apparent magnitude of a star, given luminosity and distance.
	#Inputs: AM_Dist = Distance to the star in LY
	#           AM_L = Luminosity of the star in Stellar luminosities
	#Output:    AM_m = apparent magnitude of the star
	AM_DSun = 1.58*10**(-5)#Distance to sun in LY
	AM_m = -26.73 - 2.5*np.log10(AM_L*(AM_DSun/AM_Dist)**2)
	return AM_m

def Temp2Stats(T2S_Temp):
	#convert temp (k) into rough values for stats of a main sequence star of that temp
	#Input: T2S_Temp = Star temperature in k (scalar)
	#Output: T2S_Stats = [T2S_SpecClass, T2S_Mass, T2S_Radius, T2S_Lumin, T2S_Lifetime]
	#					T2S_SpecClass = Spectral classification OBAFGKM
	#                   T2S_Mass = Mass of star in solar masses
	#                   T2S_Radius = Radius of star in solar radii
	#                   T2S_Lumin = Luminosity in stellar luminosity relation
	#                   T2S_Lifetime = Lifetime of star in Millions of years
	#  Results generally only good to the order of magnitude
	T2S_Stats = [' ',0.0,0.0,0.0,0.0]
	
	if (T2S_Temp > 48000) or (T2S_Temp<2510):#check that temp range is acceptable
		print 'Error: Temp must be < 48000 K and > 2510 K. \n Given: '+str(T2S_Temp)
		os._exit(0)

	#-------------Find T2S_SpecClass
	T2S_TRange = [30000,10000,7500,6000,5000,3500,0] #min Temp for each class
	T2S_class = ['O', 'B', 'A', 'F', 'G', 'K', 'M'] #spectral classes
	for T2S_i in range(len(T2S_TRange)):
		if (T2S_Temp >= T2S_TRange[T2S_i]):
			T2S_Stats[0] = T2S_class[T2S_i]
			break	
	
	#-------------Find T2S_Lumin
	#Some main sequence temperatures and luminosities to interpolate
	T2S_T = [2510 ,  2900, 3470, 3850,4350,4590,5080,5250,5570,5770,5800,6030,6200,6440,6890,8200,9940,9600,10500,14000,22000,30000,33000,35800, 44500, 48000]
	T2S_L = [0.003,0.0034,0.036,0.077,0.15,0.19,0.37,0.42,0.66,0.79,   1, 1.5, 2.1, 3.2, 4.3,  14,25.4,  42,   95,  500,5700,52000,97000,170000,790000,990000] 
	T2S_Stats[3] = np.interp(T2S_Temp,T2S_T,T2S_L)
	
	#-------------Find T2S_Radius
	#L = 4 pi (sigma) R^2 T^4 #Stefan-Boltzmann Law
	# So, R = sqrt(L/(4 pi (sigma) T^4))
	T2S_sigma = ((5.6070367*10**(-8))/LSun)*(RSun**2)#Stefan-Boltzmann constant (L_sun /( R_sun^2 k^-4)
	T2S_Stats[2] = np.sqrt(T2S_Stats[3]/(4*np.pi*T2S_sigma*T2S_Temp**4))
	
	#-------------Find T2S_Mass
	#L = M^a
	T2S_a = 3.5 #common value for MS stars
	T2S_Stats[1] = 1.0*T2S_Stats[3]**(1/T2S_a)
	
	#-------------Find T2S_Lifetime
	# Lifetime = 10^10*M^(-2.5) years (give answer in mill years)
	T2S_Stats[4] = (10.0**4)*T2S_Stats[1]**(-2.5)
	
	return T2S_Stats
	
def KeplersThird(KT_a, KT_M, KT_m):
	#Use kepler's 3rd law to determine the orbital period of a planet given the masses and orbital dist
	# T^2 = (1/(M+m))a^3. If a in AU, M and m in M_solar, and T in earth years
	#Inputs: KT_a = orbital distance (in AU) we are assuming a circular orbital
	#        KT_M = star mass (in solar masses)
	#        KT_m = planet mass (in solar masses)
	#Output: KT_T = orbital period (in earth years)
	KT_T = np.sqrt((KT_a**3)/(KT_M + KT_m))
	return KT_T
	
def mkRadVelPlots(mR_Obt,mR_am, mR_cV, mR_Vmax, mR_P, mR_off,mR_name):#add m dependent noise, randomize time values, add orbital offset
#Make and display simulated radial velocity plot
#Inputs: mR_Obt  = list of observation times (days)
#        mR_am   = star's apparent magnitude (for noise)
#        mR_cV   = constant velocity of star system away (m/s)
#        mR_Vmax = Maximum radial velocity of star (m/s)
#        mR_P    = Orbital period (earth years)
#        mR_off  = Orbital offset (radians further in orbit than aligned with star at start)
#Output: mR_Obt  = radial velocity observation times
#        mR_V2   = radial velocity values
#Depends on: addRndNoise
	mR_numPoints = 80 #number of dtat points to generate
	mR_V = [(mR_cV + mR_Vmax*np.sin(-2*np.pi*mR_t/(mR_P*365.25) - mR_off) )for mR_t in mR_Obt]
	mR_noiseCoeff =  0.1*(np.max(mR_V) - np.min(mR_V))#+0.00001*(10**mR_am) #scales noise based on magnitude
	mR_V2 = addRndNoise(mR_V,1,mR_noiseCoeff)

	return mR_Obt,mR_V2
	
def mkUnevenSpc(mus_min,mus_max,mus_n):	
#generate a list of values ordered (low to high), but randomly distributed and unique in a given range
#Inputs: mus_min = minimum value
#        mus_max = maximum value
#        mus_n   = number of values
#Output: mus_out = mus_n length list containing unique, sorted, random values in range
	mus_szChng = 10 #random choices will be made from a list with this factor more elements than mus_n
	mus_big = [mus_min+((mus_max-mus_min)/(mus_n*mus_szChng))*mus_k for mus_k in range(mus_n*mus_szChng)] #choices pulled from here
	mus_out = []
	for mus_i in range(mus_n):
		mus_j = rnd.choice(range(mus_n*mus_szChng-mus_i))
		mus_out.append(mus_big[mus_j])
		del mus_big[mus_j]	
	mus_out.sort()
	del mus_big,mus_i,mus_szChng,mus_j,mus_k
	return mus_out

def addRndNoise(arn_in,arn_choice,arn_scl):
#adds random noise to a list
#inputs: arn_in  = list to have noise added to items
#        arn_choice = type of dist. 0= uniform, 1=normal
#        arn_scl = Size of noise. Either +- value (for uniform dist) or std (for normal).
#outputs: arn_out = arn_in with added noise 
	if (arn_choice == 0):
		arn_out = [rnd.uniform(arn_k-arn_scl,arn_k+arn_scl) for arn_k in arn_in]
	elif (arn_choice == 1):
		arn_out = [np.random.normal(arn_k,arn_scl) for arn_k in arn_in]
	else:
		print 'ERROR:: arn_choice must be 0 or 1. Entered '+str(arn_choice)
		os._exit(0)
		
	return arn_out

def calcShadeVals(csv_P,csv_off,csv_stR,csv_pR,csv_stCM,csv_pCM):
#For each observation time, calculate the ratio of the area of the star's disk that is shaded by the planet
#Inputs: csv_P   = orbital period (years)
#        csv_off = orbital offset (radians further in orbit than aligned with star at start)
#        csv_stR = star radius (km)
#        csv_pR = planet radius (km)
#        csv_stCM = Star's orbital radius (km)
#        csv_pCM = planet's orbital radius (km)
#Output: csv_shadA = list. For each observation (area shaded by planet)/(area of star disk)
#Depends on amountCover
	csv_shadA = [] #area shadded by planet (km^2) for each observation
	#print csv_pR, csv_stR
	csv_n = 200 #number of transit observations
	
	csv_a = 3*csv_stR #Offset used to determine range to explore
	csv_minT = (-np.arcsin((csv_stR+csv_a)/(csv_stCM+csv_pCM))-csv_off+2*np.pi)*(csv_P*365.25/(2*np.pi))#time before transit
	csv_maxT = (np.arcsin((csv_stR+csv_a)/(csv_stCM+csv_pCM))-csv_off+2*np.pi)*(csv_P*365.25/(2*np.pi))#time after transit
	csv_obs = mkUnevenSpc(csv_minT,csv_maxT,csv_n)#
	
	for csv_i in range(len(csv_obs)):
		csv_stTheta = -2*np.pi*csv_obs[csv_i]/(csv_P*365.25) - csv_off#star angular position from pl aligned (starts to the left) (rad)
		csv_pTheta = csv_stTheta + np.pi#planet angular position from pl aligned (starts to the right) (rad)
		csv_sp = [csv_stCM*np.sin(csv_stTheta), csv_stCM*np.cos(csv_stTheta)]#position of star center [on trasit axis, toward viewer]
		csv_pp = [csv_pCM*np.sin(csv_pTheta), csv_pCM*np.cos(csv_pTheta)]#position of planet center [on trasit axis, toward viewer]
		#print (csv_pp[0]- csv_pR),csv_pp[0],(csv_pp[0]+ csv_pR),'|',(csv_sp[0]- csv_stR),csv_sp[0],(csv_sp[0]+ csv_stR)
		if (csv_pp[1] <= 0 ) and ((csv_pp[0]+ csv_pR)>(csv_sp[0]-csv_stR)) and ((csv_pp[0]- csv_pR)<(csv_sp[0]+csv_stR)): #planet in front of star
			csv_tempShade = amountCover(csv_stR, csv_sp[0], csv_pR, csv_pp[0], csv_n)/(np.pi*(csv_stR**2))
			csv_shadA.append(csv_tempShade)
		elif (((csv_pp[0]+ csv_pR)>(csv_sp[0]+csv_stR)) and ((csv_pp[0]- csv_pR)<(csv_sp[0]-csv_stR))): #planet bigger than star?
			csv_shadA.append(1)
		else: #star not blocked at all
			csv_shadA.append(0)
			
	return csv_shadA,csv_obs

def amountCover(ac_stR, ac_stx, ac_pR, ac_px, ac_n):
#determine area covered by planet
#Inputs: ac_stR = Radius of star (km)
#        ac_stx = Position of center of star (km)
#        ac_pR  = Radius of planet (km)
#        ac_px  = position of center of planet along transit axis (km)
#        ac_n   = number of points to integrate over
#Output: ac_out = area of overlap
#Depends on intFromList
	if ((ac_px - ac_pR) >= (ac_stx-ac_stR)) and ((ac_px + ac_pR) <= (ac_stx+ac_stR)):#planet entirely in front of star
		ac_out = np.pi*(ac_pR**2)
	else:#partial overlap
		ac_yval = []
		if ((ac_px - ac_pR) <= (ac_stx-ac_stR)): #overlap on lh side
			ac_minx = ac_stx - ac_stR
			ac_maxx = ac_px + ac_pR
			ac_xval = np.linspace(ac_minx,ac_maxx,ac_n)
			ac_int = (ac_stR**2 - ac_pR**2 + ac_px**2 - ac_stx**2)/(2*(ac_px-ac_stx))# x val for junction between disks
			for ac_i in range(ac_n):
				if (ac_xval[ac_i] <= ac_int):#star defines upper boundary
					ac_tsr = (ac_stR**2 -(ac_xval[ac_i]-ac_stx)**2)
					if (ac_tsr>0):
						ac_yval.append(np.sqrt(ac_tsr))
					else:
						ac_yval.append(0)
				else:                        #planet defines upper boundary
					ac_tsr = (ac_pR**2 -(ac_xval[ac_i]-ac_px)**2)
					if (ac_tsr>0):
						ac_yval.append(np.sqrt(ac_pR**2 -(ac_xval[ac_i]-ac_px)**2))
					else:
						ac_yval.append(0)
		elif ((ac_px + ac_pR) >= (ac_stx+ac_stR)): #overlap on rh side
			ac_minx = ac_px - ac_pR
			ac_maxx = ac_stx + ac_stR
			ac_xval = np.linspace(ac_minx,ac_maxx,ac_n)
			ac_int = (ac_stR**2 - ac_pR**2 + ac_px**2 - ac_stx**2)/(2*(ac_px-ac_stx))# x val for junction between disks
			for ac_i in range(ac_n):
				if (ac_xval[ac_i] <= ac_int):#planet defines upper boundary
					ac_tsr = (ac_pR**2 -(ac_xval[ac_i]-ac_px)**2)
					if (ac_tsr>0):
						ac_yval.append(np.sqrt(ac_tsr))
					else:
						ac_yval.append(0)
				else:                        #Star defines upper boundary
					ac_tsr = (ac_stR**2 -(ac_xval[ac_i]-ac_stx)**2)
					if (ac_tsr>0):
						ac_yval.append(np.sqrt(ac_tsr))
					else:
						ac_yval.append(0)
		else:
			print 'ERROR:: something went wrong in amountCover. Not lh, rh or full cover.'
			os._exit(0)
							
		ac_out = 2*intFromList(ac_xval,ac_yval) #integrate output list, then multiply by 2 because only looking at top half
	return ac_out
		
	
def mkTransPlot(mtp_obs,mtp_sv,mtp_L,mtp_d,mtp_name):
#Make plot of transit
#inputs: mtp_obs = observation times(days)
#        mtp_sv = ratio of (shaded area)/(area of star disk) for each observation
#        mtp_L = Luminosity of star (Solar Luminosity)
#        mtp_d = distance to star (Ly)
#outputs: mtp_obs = transit observation values
#         mtp_m2 = transit magnitude values
	mtp_m = [AppMag(mtp_d, mtp_L*(1-mtp_k)) for mtp_k in mtp_sv]
	mtp_noise = 0.05*(np.max(mtp_m) - np.min(mtp_m))#scale for acceptable noise
	mtp_m2 = addRndNoise(mtp_m,1,mtp_noise) 
	return mtp_obs,mtp_m2
	
def determineHabitableZone(dhz_L):
#Determine the rough range of radii from star that contain the habitable zone
#Inputs: dhz_L = star's luminosity (L_sun)
#Output: dhz_out = radii bounding habitable zone [inner, outer] (AU)
	dhz_out = [0,0]
	dhz_out[0] = np.sqrt(dhz_L/1.1)
	dhz_out[1] = np.sqrt(dhz_L/0.53)
	return dhz_out

def exoPlanetSimulation(es_stT,es_stD,es_stV,es_stO,es_pR,es_pOR,es_pM,es_n,es_name):
#Simulate exoplanet and star based on stats, then make star radial velovity plot and transit light curve
#Inputs: es_stT = star temperature (k)
#        es_stD = distance to star (Ly)
#        es_stV = velocity away of star system (m/s)
#        es_stO = initial orbital offset (radians past aligned)
#        es_pR  = planet radius (earth radii)
#        es_pOR = planet orbital radius (AU)
#        es_pM  = planet mass (earth mass)
#        es_n   = number of observations to be made
#        es_name= planet name (for plots)
#output: es_stStatsFull = stats for star: [temperature(k),spec class,Mass(M_sun),Radius (R_sun),Luminosity(L_sun), Lifetime(Myr), apparent magnitude]
#        es_pStats      = stats for planet: [name, mass (M_earth),radius(R_earth),density(kg/m^3),orbital radius (AU),orbital period (yr)]
#        es_sysStats    = stats for system: [distance(Ly), inner hab zone(AU), outer hab zone(AU), velocity away (m/s), orbital offset (rad)]
#Depends on: Temp2Stats,AppMag,KeplersThird,mkUnevenSpc,mkRadVelPlots,calcShadeVals,mkTransPlot,determineHabitableZone

	global radObs,tranObs,radVals,tranVals

	es_stStatsFull = [es_stT,' ',0.0,0.0,0.0,0.0,0.0]
	es_pStats = [es_name,es_pM,es_pR,0.0,es_pOR,0.0]
	es_sysStats = [es_stD,0.0,0.0,es_stV,es_stO]

	es_plMassSolar = es_pM*(MEarth/MSun)#mass of planet in solar masses
	es_stStats = Temp2Stats(es_stT) #Get star stats from temp [SpecClass, Mass, Radius, Lumin, Lifetime]
	es_stm = AppMag(es_stD, es_stStats[3]) #Calc apparent magnitude
	es_plOrPer = KeplersThird(es_pOR,es_stStats[1],es_plMassSolar)#Orbital Period in years
	es_stCM = (es_plMassSolar*es_pOR)/(es_plMassSolar+es_stStats[1])#Star distance from center of mass
	es_plCM = es_pOR - es_stCM #planet distance from center of mass
	es_stVmax = (2* np.pi* es_stCM * Au_conv)/(es_plOrPer*SinYr)#maximum tangential v of star (m/s)
	es_obs_time = mkUnevenSpc(0,1.5*es_plOrPer*365.25,es_n)#Observation times for radial velocity(days)
	
	es_rOb,es_rV = mkRadVelPlots(es_obs_time,es_stm, es_stV, es_stVmax, es_plOrPer,es_stO,es_name)#Simulate radial velocity curves
	radObs.append(es_rOb)#save for later
	radVals.append(es_rV)
	
	es_shadeVal,es_obs_timeTrans = calcShadeVals(es_plOrPer,es_stO,(es_stStats[2]*RSun)/1000.,es_pR*REarth/1000.,es_stCM*Au_conv/1000.,es_plCM*Au_conv/1000.)
	es_tOb,es_tV = mkTransPlot(es_obs_timeTrans,es_shadeVal,es_stStats[3],es_stD,es_name)#make the transit plot
	tranObs.append(es_tOb)#save for later
	tranVals.append(es_tV)
	
	es_hz = determineHabitableZone(es_stStats[3])#rough range of habitable zone for star 

	es_stStatsFull[1:6] = es_stStats #record stats for output
	es_stStatsFull[6] = es_stm
	es_pStats[3] = (es_pM*MEarth)/((4./3)*np.pi*(es_pR*REarth)**3)
	es_pStats[5] = es_plOrPer
	es_sysStats[1:3] = es_hz
	
	return 	es_stStatsFull,es_pStats,es_sysStats
	
def pickSystem(ps_t,ps_h):
#pick stats for the planetary system
#inputs: ps_t = terrestrial (1=terrestrial, 0=jovian)
#        ps_h = habitable (0=habitable zone, 1=not in habitable zone)
#output: ps_out = [star temp, distance, system velocity, offset, planet R, planet M, orbit R]

	ps_out = [0,0,0,0,0,0,0]
	ps_out[1] = rnd.uniform(4,10)#random distance
	ps_out[2] = np.random.normal(0,10000)#random system velocity
	ps_out[3] = rnd.uniform(0,np.pi)#random orbital offset
	
	if (ps_t == 1): #Must be terrestrial 
		ps_out[5] = rnd.uniform(0.1,2)#mass on range near terrestrial planet (earth mass)
		ps_den = 5500*rnd.uniform(0.7,1.2)#random terrestrialish Density (kg/m^3)
		ps_out[4] = (((3.0*MEarth*ps_out[5])/(4.0*np.pi*ps_den))**(1./3))/REarth
		
	else: #Is jovian
		ps_out[5] = rnd.uniform(14,400)#mass on range near jovian planets (earth mass)
		ps_den = 5500*rnd.uniform(0.1,0.3)#jovianish Density (kg/m^3)
		ps_out[4] = (((3.0*MEarth*ps_out[5])/(4.0*np.pi*ps_den))**(1./3))/REarth
		
	if (ps_h == 1): #Must be in habitable zone
		ps_out[0] = rnd.uniform(3300,8000)#random temp value of reasonably behaved star
		ps_st = Temp2Stats(ps_out[0])#calc star stats
		ps_hz = determineHabitableZone(ps_st[3])
		ps_out[6] = rnd.uniform(ps_hz[0],ps_hz[1])#determine orbital radius in habitable zone
		
	else: #Not necessarily habitable... but might be...
		ps_temp = np.random.normal(3000,12000)#random temp value for star
		ps_out[0] = np.sqrt((ps_temp - 2700)**2)+2700 #made to fit range, low mass more probable
		ps_st = Temp2Stats(ps_out[0])#calc star stats
		ps_out[6] = rnd.uniform(0.2,60.0)#determine orbital radius 
	return ps_out

def PlanMission():
#open page to plan mission
	def quitPM():
		#close window and quit
		top3.destroy()
		raise SystemExit(0)
	def returnPM():
		#close planning window and return to display plots
		top3.destroy()
		DisplayPlots()
	def OnLaunch():
		#what to do when the launch button has been pressed
		ol_Sys = SysChoice.get()#retrieve planet choice
		
		if ('None' in ol_Sys): #chose to go nowhere
			ol_ResultMessage = '\nYou chose to go nowhere because you are boring.\n Once all the other humans have left for places unkown, you sit alone on the poisoned rock they left behind'
			ol_stayOpt = [' and live out your days as a grumpy, old hermit.',' and choke on the exhaust from their departing ships.', ' and wish you had made something of your life.', ' and cry while watching old episodes of Friends.','. At first you enjoy the peace and quiet, but one day you choke on your pie and, since there is nobody around to help, you suffocate and die. \nYou can rest assured knowing that you at least provided a good meal for the hungry scavangers.']
			ol_ResultMessage = ol_ResultMessage + rnd.choice(ol_stayOpt)
			verdictMessage.set('\n\nYOU FAIL TO SURVIVE ON A NEW PLANET... or do anything of interest.')

		else: #To infinity and BEYOND!
			buttonindex = SysNames.index(ol_Sys)
			ol_ResultMessage = '\nYou choose to go to system '+SysNames[buttonindex]+'. \n'
			systemStatsMessage = 'Planet stats:\n mass = '+str(round(plS[buttonindex][1],2))+' earth masses\n radius = '+str(round(plS[buttonindex][2],2))+' earth radii\n density = '+str(round(plS[buttonindex][3],2))+' kg/m^3 \n           = '+str(round((plS[buttonindex][3]/5500),2))+' earth densities \norbital period = '+str(round(plS[buttonindex][5],2))+' years \n orbital radius = '+str(round(plS[buttonindex][4],2))+' AU \n habitable zone: '+str(round(sysS[buttonindex][1],2))+' to '+str(round(sysS[buttonindex][2],2))+' AU'
			SysStatMes.set(systemStatsMessage)
			if (plS[buttonindex][3] < 0.4*5500): #Jovian planet
				ol_ResultMessage = ol_ResultMessage + 'When you arrive, you find a beautiful Jovian planet.\n Unfortunately you were not looking for a jovian planet.\n'
				ol_JovOpt = ['You try to land anyway and the ship burns up in the atmosphere.','Lacking the fuel to return to earth, you sit in your ship in orbit, wondering where you went wrong. When people start starving they choose to eat you first because this was all your fault.','You try to settle one of its moons instead, but there is no atmosphere and you quickly suffocate.']
				ol_ResultMessage = ol_ResultMessage + rnd.choice(ol_JovOpt)
				verdictMessage.set('\n\nYOU FAIL TO SURVIVE ON A NEW PLANET.')
				
			else:
				ol_ResultMessage = ol_ResultMessage + 'When you arrive, you find a terrestrial planet.\n'
				if (plS[buttonindex][4]<sysS[buttonindex][1]):#too close for habitable
					ol_ResultMessage = ol_ResultMessage + 'Unfortunately it is too close to the star to be habitable and your skin melts off.'
					verdictMessage.set('\n\nYOU FAIL TO SURVIVE ON A NEW PLANET.')
				
				elif (plS[buttonindex][4]>sysS[buttonindex][2]):#too far to be habitable
					ol_ResultMessage = ol_ResultMessage + 'Unfortunately it is too far from the star to be habitable and you starve to death while clutching the frozen corpse of a loved one.'
					verdictMessage.set('\n\nYOU FAIL TO SURVIVE ON A NEW PLANET.')
					
				else:
					ol_ResultMessage = ol_ResultMessage + 'It is in the habitable zone of its star.'
					ol_bt = buildChoi.get()#retrieve building type
					if plS[buttonindex][4]  < (((sysS[buttonindex][2]-sysS[buttonindex][1])/3)+sysS[buttonindex][1]):#inner third of habitable zone
						ol_ResultMessage = ol_ResultMessage + ' Being in the inner third of the habitable zone, you are in for some high temperatures.'
						if ('hot weather' in ol_bt):#brought warm weather buildings
							ol_ResultMessage = ol_ResultMessage + ' Good thing the buildings you brought the supplies to construct were designed for such a climate! \nThough there are certain times of day when you can not go outside, and you have to be careful to protect your crops from the harsh sun, you came prepared to do this and, thus, you flourish in your new home!'
							verdictMessage.set('\n\nYOU SUCCESSFULLY SURVIVE ON A NEW PLANET!')
							
						else: #brought a different building	
							ol_ResultMessage = ol_ResultMessage + ' Unfortunately the supplies you brought do not permit you to build the appropriate structures to provide adequate protection. You struggle for a while to survive, but your numbers dwindle over time as the crops fail and people succumb to heat stroke.\n Years later, a more prepared group finds its way to your planet and discovers your bones scattered among the dust of a wildly inappropriate building.\n "What could they have possibly been thinking?" these new explorers ask themselves as they seal it up and begin their own colony.'
							verdictMessage.set('\n\nYOU FAIL TO SURVIVE ON A NEW PLANET.')
							
					elif plS[buttonindex][4]  > (sysS[buttonindex][2] - ((sysS[buttonindex][2]-sysS[buttonindex][1])/3)):#outer third of habitable zone
						ol_ResultMessage = ol_ResultMessage + ' Being in the outer third of the habitable zone, you are in for some low temperatures.'
						if ('cold weather' in ol_bt):#brought cold weather buildings
							ol_ResultMessage = ol_ResultMessage + ' Good thing the buildings you brought the supplies to construct were designed for such a climate! \nThough there are certain times when you can not go outside, and you have to be careful to protect your crops from the freezing cold, you came prepared to do this and, thus, you flourish in your new home!'
							verdictMessage.set('\n\nYOU SUCCESSFULLY SURVIVE ON A NEW PLANET!')
							
						else: #brought a different building	
							ol_ResultMessage = ol_ResultMessage + ' Unfortunately the supplies you brought do not permit you to build the appropriate structures to provide adequate protection. You struggle for a while to survive, but your numbers dwindle over time as the crops fail and people freeze to death.\n Years later, a more prepared group finds its way to your planet and discovers your frozen bodies huddled among the rooms of a wildly inappropriate building.\n "What could they have possibly been thinking?" these new explorers ask themselves as they seal it up and begin their own colony.'
							verdictMessage.set('\n\nYOU FAIL TO SURVIVE ON A NEW PLANET.')
							
					else:#middle third of the habitable zone
						ol_ResultMessage = ol_ResultMessage + ' Being in the middle third of the habitable zone, you are in for pleasant weather. It is hard work, but your crops flourish in this beautiful new environment and so do you.'
						verdictMessage.set('\n\nYOU SUCCESSFULLY SURVIVE ON A NEW PLANET!')
		
		ol_paramSimp = [0,0] #basic planet survival parameters0=no, 1= yes [terrestrial?, habitable zone?]
		
		ResultMessage.set(ol_ResultMessage )


	top3 = Tk.Tk()

	
	backbutton = Tk.Button(top3, text = " Return to choices ",activeforeground='white',activebackground='gray', command = returnPM).grid(row = 1, column = 1)
	quitbutton = Tk.Button(top3, text = " Quit ",activeforeground='white',activebackground='gray', command = quitPM).grid(row = 1, column = 2)
	
	
	welcomeVar = Tk.StringVar()
	SysChoice = Tk.StringVar()
	ResultMessage = Tk.StringVar()
	verdictMessage = Tk.StringVar()
	SysStatMes = Tk.StringVar()
	buildChoi = Tk.StringVar()
	
	welcome = Tk.Message(top3,textvariable=welcomeVar,width=800,justify='center')
	welcomeVar.set('Now that you have decided on a destination, it is time to plan your mission.\n When you are done making selections, press the Launch button.')
	welcome.grid(row=2,columnspan=3)
	
	SysChoice.set(SysNames[0]) # set the default option

	pm_choices = [k for k in SysNames]#choices for pull-down menu
	pm_choices.append('None')#add "none" for the lazies
 
	planetChoice = Tk.OptionMenu(top3, SysChoice,*pm_choices)
	Tk.Label(top3, text='Choose your destination: ').grid(row = 3, column = 1)
	planetChoice.grid(row = 3, column =3)
	
	buildChoices = ['normal','cold weather (advised for outer third of habitable zone)','hot weather (advised for inner third of habitable zone)']
	buildChoi.set(buildChoices[0])
	BuildingChoice = Tk.OptionMenu(top3, buildChoi,*buildChoices)
	Tk.Label(top3, text='What types of buildings will you want to build there? ').grid(row = 5, column = 1)
	BuildingChoice.grid(row = 5, column =3)
	
	Launchbutton = Tk.Button(top3, text = " Launch! ",activeforeground='red',activebackground='gray', command = OnLaunch).grid(row = 6, column = 2)
	
	Result = Tk.Message(top3,textvariable=ResultMessage,width=300,justify='left')
	Result.grid(row = 7, columnspan = 7)
	
	Verdict = Tk.Message(top3,textvariable=verdictMessage,width=300,justify='center', foreground = 'red')
	Verdict.grid(row = 8, columnspan = 7)
	
	ActualStats = Tk.Message(top3,textvariable=SysStatMes,width=300,justify='left')
	ActualStats.grid(row = 7, column = 9)
	
	top3.mainloop()

def intermediary():
#do some window closing between top and top2
		top.destroy()#close window
		DisplayPlots()
	
def DisplayPlots():

	def on_key_event(event):
		if event.key in SysNames:
			buttonindex = SysNames.index(event.key)
			print 'You chose to examine ',SysNames[buttonindex]
			
			oke_Fig = plt.figure(2+buttonindex)#Radial velocity plots
			plt.subplot(2,1,1)
			plt.plot(radObs[buttonindex],radVals[buttonindex],'k.')
			plt.xlabel('time (days)')
			plt.ylabel('velocity away (m/s)')
			plt.title('Radial Velocity: '+plS[buttonindex][0])
			plt.tight_layout()
		
			plt.subplot(2,1,2)
			#mAveTemp = np.average(tranVals[buttonindex][0:10])
			#TempIntensity = [10**(0.4*(mAveTemp-k)) for k in tranVals[buttonindex]]#convert magnitude to relative intensity.
			plt.plot(tranObs[buttonindex],tranVals[buttonindex],'k.')
			plt.xlabel('time (days)')
			plt.ylabel('apparent magnitude')
			plt.title('Transit: '+plS[buttonindex][0])
			plt.ylim([np.max(tranVals[buttonindex]),np.min(tranVals[buttonindex])])
			plt.tight_layout()
			
			oke_Fig.show()
			
		elif event.key == 'escape':#quit
			top2.quit()
			top2.destroy()
			quit()
		elif event.key == 'enter':#done, next page
			top2.quit()
			top2.destroy()
			PlanMission()
		elif event.key == 'x':#directions
			DirectMess = 'You can Determine things in the following way:\n\n The orbital period is the amount of time it takes for the radial velocity plot to go through 1 cycle.\n\n'
			DirectMess = DirectMess + 'The orbital velocity will be the amplitude of the radial velocity plot (half the difference between max and min)\n\n'
			DirectMess = DirectMess + "The orbital distance can be determined using Kepler's 3rd law: P^2 = (1/M)a^3 where a is the semimajor axis and M is the mass of the star.\n\n"
			DirectMess = DirectMess + 'The planet mass can be determined using m = (M*(M_e/M_s)*v*P*ns)/(2 pi ((M*P^2)^(1/3))*nm) where M is the Mass of the star in solar masses,M_e is the mass of the earth in kg,M_s is the mass of the sun in kg, v is the orbital velocity in m/s, P is the orbital period in years, ns is the number of seconds in a year, and nm is the number of meters in an AU \n'
			DirectMess = DirectMess + 'M_s = 1.989*10^30 kg, M_e = 5.972*10^24 kg, ns = 31469940 s, nm = 1.496*10^11 m \n\n'
			DirectMess = DirectMess + 'The planet radius is determined using the following (\delta m) = (r_planet/r_star)^2 where (\delta m) is the change in magnitude during the transit.\n\n'
			DirectMess = DirectMess + 'Density is the mass divided by the volume (4/3 pi r^3). Hint: do this in kg/m^3 first, then divide by the density of the earth (5500 kg/m^3). \n\n'
			DirectMess = DirectMess + 'Densities less than about 0.4 times that of the earth are Jovian planets \n\n'
			DirectMess = DirectMess + 'The rough approximation of the habitable zone of a star is between sqrt(L/1.1) and sqrt(L/0.53) AU, if L is the luminosity of the star in solar luminosities.\n\n'
			
			
			tkMessageBox.showinfo('Directions', DirectMess)
		else:
			print event.key,' is not a valid option.'
		
	#top.quit()
	top2 = Tk.Tk()
	top2.attributes("-fullscreen", True)
	
	fig = plt.figure(1,figsize=(10, 8))#figures to show
	plt.text(0,0,"For some reason you have to click on this window before anything will work.\n To quit press ESCAPE. To examine the plots for any planet more closely, push the letter on your keyboard associated with its name.\n When you are done examining the options, hit ENTER. You will make your choice on the next screen.")
	plt.axis('off')
	
	for dp_i in range(nSys):#for each system
	
		a = fig.add_subplot(nSys+1,4,4*dp_i+1)
		a.text(0,0.5,plS[dp_i][0],fontsize = 20)
		plt.axis('off')
		
		b = fig.add_subplot(nSys+1,4,4*dp_i+2)#Radial velocity plots
		plt.plot(radObs[dp_i],radVals[dp_i],'k.')
		plt.xlabel('time (days)')
		plt.ylabel('velocity away (m/s)')
		plt.title('Radial Velocity: '+plS[dp_i][0])
		plt.tight_layout()
		
		c = fig.add_subplot(nSys+1,4,4*dp_i+3) #Transit plots
		#mAveTemp = np.average(tranVals[dp_i][0:10])
		#TempIntensity = [10**(0.4*(mAveTemp-k)) for k in tranVals[dp_i]]#convert magnitude to relative intensity.
		plt.plot(tranObs[dp_i],tranVals[dp_i],'k.')
		plt.xlabel('time (days)')
		plt.ylabel('apparent magnitude')
		plt.title('Transit: '+plS[dp_i][0])
		plt.ylim([np.max(tranVals[dp_i]),np.min(tranVals[dp_i])])
		
		dp_txtInfo = 'Star stats:\nclass: '+stS[dp_i][1]+'\nM = '+str(round(stS[dp_i][2],2))+' $M_{\odot}$\nR = '+str(round(stS[dp_i][3],2))+ '$R_{\odot}$\nL = '+str(round(stS[dp_i][4],2))+ '$L_{\odot}$\n'
		d = fig.add_subplot(nSys+1,4,4*dp_i+4)
		d.text(0,0,dp_txtInfo)
		plt.axis('off')
	
	canvas = FigureCanvasTkAgg(fig, master=top2)
	canvas.show()
	canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

	canvas.mpl_connect('key_press_event', on_key_event)
	
	top2.mainloop()
	
	
#################################################################################
###   MAIN   ###   MAIN   ###   MAIN   ###   MAIN   ###   MAIN   ###   MAIN   ###
#################################################################################
#major assumptions:
#Circular orbits. Viewed perfectly edge-on. No limb darkening. Only 1 planet per system.
	
global LSun,RSun,REarth,MSun,MEarth,Au_conv
global sysS,plS,stS,nSys, SysNames
global radObs,tranObs,radVals,tranVals
	
radObs = []
tranObs = []
radVals = [] 
tranVals = []

#constants
LSun = 3.828*10**(26) #Luminosity of sun in W
RSun = 6.957*10**(8) #Radius of sun in m
REarth = 6.3781*10**(6) #radius of earth in m
MSun = 1.989*10**(30) #Mass of sun (kg)
MEarth = 5.972*10**(24) #Mass of Earth (kg)
Au_conv = 149597870700.0 #1 AU = 149597870700 m
SinYr = 3600*24*365.25 #How many seconds in a sidereal year

numObs = 100#Number of observations made for system
nSys = 4#Number of systems to simulate

SysNames = ['a','b','c','d']
SysProp = [[1,1],[1,0],[0,1],[0,0]] #one of each type of system
plOrder = [0,1,2,3] #indices for SysProp
rnd.shuffle(plOrder) #Randomize indices

top = Tk.Tk()
var1 = Tk.StringVar()
var2 = Tk.StringVar()
var3 = Tk.StringVar()
ProgVar = Tk.StringVar()
Logo = Tk.PhotoImage(file = 'FreighterWBg.gif')

sysS = []#list of stats for each system
plS = []#List of stats for each planet
stS = []#List of stats for each star

logoIm = Tk.Label(image = Logo)
logoIm.grid(row=0,column=1)
label = Tk.Message(top,textvariable=var1,width=300,justify='center')
var1.set('Welcome to Exoplanet Survivor!\n Code written by Richard D Mellinger Jr')
label.grid(row=1,columnspan=3)
label2 = Tk.Message(top,textvariable=var2,width=300,justify='left')
var2.set("\n You will be shown the radial velocity curves and transit light curves for several simulated exoplanet systems. \n At least one of these will be a terrestrial planet in its star's habitable zone. \n Use the plots to calculate basic properties of the planets and decide which planet to move to.\n\n Once your destination has been chosen, you will plan your mission.\n\n Then we'll see if you survive!")
label2.grid(row=2,columnspan=3)

ProgVar.set('  ')
Tk.Label(top,textvariable=ProgVar).grid(row = 4, column = 1)

for i in range(nSys):#for each system
	sysLs = pickSystem(SysProp[plOrder[i]][0],SysProp[plOrder[i]][1])#[star temp, distance, system velocity, offset, planet R, planet M, orbit R]
	stStat,plStat,sysStat = exoPlanetSimulation(sysLs[0],sysLs[1],sysLs[2],sysLs[3],sysLs[4],sysLs[6],sysLs[5],numObs,SysNames[i])
	ProgVar.set('Generating systems... '+str(round((100.*(i+1)/nSys),1))+' % complete')#Update Progress
	sysS.append(sysStat)
	plS.append(plStat)
	stS.append(stStat)
	
ProgVar.set(' ')
Tk.Button(top, text = " Let's go! ",activeforeground='white',activebackground='gray', command = intermediary).grid(row = 5, column = 1)

top.mainloop()
