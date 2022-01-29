## Date: 01-13-22 
## THis is a modification to the earlier code. In the earlier code average current from the APD was input in toe the thermal calculation. From this average current the average optical intensity to generate that current was obtained and from that the PEAK optical irradiance which created the damage was calculated. 

## One issue with that approach is that every thing is calculated at an average basis. So when the time resolution was changed the damage did not show up.

#this code takes the peak power per pulse as an input and from that calculates the thermal equation.

#adding photocurrent generation via a pulse 
#optical pulses


### date: 01-24-2022

# modeling the effect of a series resistance on the thermal effect of the device

## experimental results shows that passives up to 10K does not change the optical damage. Also, limiting the current from the SMU to 0.1 mA does not have a big effect.


# R will limit the maximum dark current in the circuit and induce a fluctualation of the vbias due to the drop across the R and eventual time it takes to get re-armed.
## The primary effect of the above will be that the bias across the APD is not constant and will affect the heat generated in the device.


# Discharge across the interal resistance of the APD and charge using the R15

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
font = {'family': 'courier','weight':'normal','size':12}
matplotlib.rc("font",**font)

def heatflow(peakpower,PRF,R15):  #peak power in kW/cm2 #H20E bond h = 0.29, _jphmax is given in A/cm2, 5mA@200-um2 area is 15A/cm2
	print(PRF,peakpower)

	#R15 = 10E6

	Vbr = 41.7

	T0 = 25
	h = 0.29
	resolution = 3.7E-9
	t = np.arange(0,1E-2,resolution)
	T = np.zeros(np.size(t))

	J0 = 100*1E-6 		# A/cm2
	V = 40.0 			# V
	D = 0.05			# 1/C
	Ac = np.pi*50*50*1E-8/4	# cm2
	h = h/50E-4  			# W/cm2.K

	M = 10 #gain

	#	
	Idarklimit = V/R15 #A
	#print("Due to the series resistance of {} the dark current Id is limited under {} A".format(R15,Idarklimit))
	Capd = 128 * 1E-15 #F
	Tau = R15*Capd #s
	print("tau: {}".format(Tau))
	#print(" The time constant of the arming is {} s".format(Tau))
	JdarkLimit = Idarklimit/(np.pi*200*200*0.25*1E-8) #cm2 #Limit of the current drawn at 3.7 ns 
	print("Max dark current density due to R15  current limit {} A/cm2".format(JdarkLimit))
	##


	C = 1E-4 #define									# J/g.K, use the volumen and density to get the g)
	C = 5.5*np.pi*200*200*3*1E-12/4
	Aj = np.pi*200*200*1E-8/4
	Tambient = 20.0

	Tpulse = resolution #4E-9
	
	QE = 0.8
	Responsivity = QE * 1240/1550 
	averagepower = peakpower * 1000 * Tpulse * PRF # energy per second
	averagecurrentgenerated = averagepower * Responsivity * M 
	Jphmax = Responsivity * M * peakpower * 1000 #since peakpower is in Kw

	#if (J0 >= JdarkLimit):
	#	print(">>> Dark current is limited to {}".format(JdarkLimit))
	#	J0 = JdarkLimit

	TransientSlots = int(Tau/resolution)

	def Gain(Vbr,Vin):
		if R15 >= 1E6: #there might be gain issue at high series resistance
			return 1;
		else:
			return((1-Vin/Vbr)**(-1*0.73))

	def func(Temp,Time,fire,newVoltage):
		if(fire == 1):
			#return((V*Aj*J0*np.exp(D*(Temp-Tambient)) + (V*Aj*Jphmax) - (h*Ac*(Temp-Tambient)))/C)
			M = Gain(Vbr,newVoltage)
			#print("Gain: {}".format(M))
			Jphmax = Responsivity * M * peakpower * 1000 #since peakpower is in Kw, calculate gain for every voltage at the firing point
			##the dakcurrent has an upper limit as well
			if (J0*np.exp(D*(Temp-Tambient)) <= JdarkLimit):
				return((newVoltage*Aj*J0*np.exp(D*(Temp-Tambient)) + (newVoltage*Aj*Jphmax) - (h*Ac*(Temp-Tambient)))/C)
			else:
				#print("Hits the dark limit")
				return((newVoltage*Aj*JdarkLimit + (newVoltage*Aj*Jphmax) - (h*Ac*(Temp-Tambient)))/C)
		
		else:
			#return((V*Aj*J0*np.exp(D*(Temp-Tambient)) - (h*Ac*(Temp-Tambient)))/C)
			if (J0*np.exp(D*(Temp-Tambient)) <= JdarkLimit):
				return((newVoltage*Aj*J0*np.exp(D*(Temp-Tambient)) - (h*Ac*(Temp-Tambient)))/C)
			else:
				#print("Hits the dark limit")
				return((newVoltage*Aj*JdarkLimit - (h*Ac*(Temp-Tambient)))/C)


	T[0] = T0

	fire_step = 0
	fire_interval = int((1/PRF)/resolution)
	print("fire interval=%d"%fire_interval)

	
	transientstep = 0


	#Discharge Tau
	#dependent on the dynamic resistance the discharge will vary
	#dynamics resistance dependes on the series resistance from emperical data
	
	Rdynamics = 1590 #ohms without passives

	if (R15 >= 3000) and (R15 < 1E4):
		Rdynamics = 7110
	elif (R15 >= 1E4) and (R15 < 1E6):
		Rdynamics = 16300
	elif (R15 >= 1E6):
		Rdynamics = 87900

	dischargeTau = Rdynamics * 128 * 1E-15 

	if (dischargeTau >= Tau):
		TransientSlots = int(dischargeTau/resolution)
	else:
		TransientSlots = int(Tau/resolution)
		
	print("Timeslots where the diode is off or rechargin {}".format(TransientSlots))	

	##get APD voltage
	V_APD = []

	def ApdChargeVoltage(timestep):
		return (V*(1-np.exp(-(timestep*resolution)/Tau)))
	def ApdDischargeVoltage(timestep):
		return (V*np.exp(-(timestep*resolution)/dischargeTau))

	voltage = V	

	for i in range(1,np.size(t)):
		
		if (fire_step == fire_interval):
			if(transientstep < TransientSlots) : #if the timestep is less than Tau from the point the laser fired discharge the APD
				voltage = ApdDischargeVoltage(transientstep)
			else:
				voltage = ApdChargeVoltage(transientstep)

			T[i] = T[i-1] + (func(T[i-1],t[i-1],1,voltage)*(t[i] - t[i-1]))
			transientstep =0
			fire_step = 0
			V_APD.append(voltage)
	
		else:
			fire_step = fire_step + 1
			if(transientstep < TransientSlots) : #if the timestep is less than Tau from the point the laser fired discharge the APD
				voltage = ApdDischargeVoltage(transientstep)
			else:
				voltage = ApdChargeVoltage(transientstep)

			T[i] = T[i-1] + (func(T[i-1],t[i-1],0,voltage)*(t[i] - t[i-1]))
			transientstep = transientstep + 1
			V_APD.append(voltage)

	#plt.semilogx(t[1:],V_APD)
	#return

	current = np.pi*200*200*0.25*1E-8*Jphmax*1E3 
	#l = "PRF=%iKHz, Pph = %1.3f kW/cm2"%(PRF/1E3,peakpower)
	l = "PRF=%iKHz,Pph=%1.1fkW/cm2,R=%1.1gOhm"%(PRF/1E3,peakpower,R15)
	#l = "PRF=%iKHz, I = %1.1f mA, Pph = %1.3f kW/cm2"%(PRF/1E3,current,peakpower)
	plt.semilogx(t,T,label=l)
	#plt.plot(t,T,label=l)
	return 

def __main__():


	#heatflow(60,PRF=1E4, R15=1E4) #A/cm2, PRF
	#heatflow(0.3,PRF=970*1000, R15=0.1) #A/cm2, PRF
	
	for i in np.arange(60,90,10):
		heatflow(i,PRF=1E4,R15=1E6)

	#for i in np.arange(10,30,5):
	#	heatflow(i,PRF=970*1000,R15=1E6)
		
	#for i in np.arange(1,50,10):
	#	heatflow(i,970*1000,R15=10E6)

	#for i in np.arange(0.1,20,10):
	#	heatflow(i,3571*1000)

	#for i in  range(3,6,1):
	#	heatflow(15,10**i)



	plt.ylim(0,500)
	
	plt.xlabel("Time (s)")
	plt.ylabel("Temp (C)")
	plt.legend()
	plt.grid(b="True",which="major",color="k",linestyle='-')
	plt.minorticks_on()
	plt.grid(b="True",which="minor",color="k",linestyle='--',alpha=0.2)
	plt.show()


__main__()
