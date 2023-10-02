import chargeexchange as CX

#Constants    
mass_ion=1.0
mass_neutral=1.0

#Conversion of units
SI_to_Nektar_n=3.33333e-19

#Initial Conditions
weight=1.65E+16
initial_neutral_velocity=0.2
initial_fluid_momentum=0.8
initial_fluid_density=3.3e18

#Timestep info
dt_nektar=0.005 
dt_SI=1.66667e-07

#Rate info
CXrate = 5E-14

IC=CX.InitialConditions(initial_fluid_density,initial_fluid_momentum,weight,initial_neutral_velocity)
newPart = CX.Particle(mass_neutral,weight,initial_neutral_velocity)
newIonFluid = CX.IonFluid(mass_ion,initial_fluid_density,initial_fluid_momentum)

CXrunner = CX.runner(IC,newIonFluid,newPart,250,dt_nektar,dt_SI,CXrate)
CXrunner.runCX()
