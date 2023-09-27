#This python script executes charge exchange
import matplotlib.pyplot as plt #Used for plotting 
import numpy as np #Used for normal distribution
import unittest #Used for unit testing of code
R_CE = 5E-14 #Units m^3 s^-1
w = 1.65E+16 #No units
dt = 1.66667e-07 #Units Seconds
timestep=0 #No units
max_neutrals_density=3.3e18;
max_timesteps=800 #No units
fluid_velocity_sd=0.0

mass_neutral=1.0 #Nektar units
mass_ion=1.0 #Nektar units

neutrals_density=np.zeros(max_timesteps+1) #Nektar Units
ion_density=np.zeros(max_timesteps+1) #Nektar Units
fluid_velocity=np.zeros(max_timesteps+1) #Nektar Units
fluid_momentum=np.zeros(max_timesteps+1) #Nektar Units

neutrals_velocity=np.full(max_timesteps+1,0.5) #Nektar Units


n_to_nektar=3.33333e-19 #Conversion factor from SI to nektar units
dt_nektar=0.005 #Timestep in Nektar units
neutrals_momentum=np.zeros(max_timesteps+1) #Nektar units
time=np.zeros(max_timesteps+1) #Time in Nektar Units
total_momentum_injected=np.zeros(max_timesteps+1) #Nektar Units

particle_injection_timestep_cutoff=100 #Timestep at which code stops injecting particles

class InitialConditions:
    def __init__(self,fluid_density,fluid_momentum,neutrals_density,neutrals_velocity):
        self.fluid_density=fluid_density
        self.fluid_momentum=fluid_momentum
        self.neutrals_density=neutrals_density
        self.neutrals_velocity=neutrals_velocity
        
#Add mass of a single neutral        
class ParticleDat:
    def __init__(self, max_timesteps, initial_neutrals_density, initial_neutrals_velocity):
        self.max_timesteps=max_timesteps #No units
        self.neutrals_density=np.zeros(max_timesteps+1)
        self.neutrals_velocity=np.zeros(max_timesteps+1)
        self.neutrals_density[0]=initial_neutrals_density #Units SI
        self.neutrals_velocity[0]=initial_neutrals_velocity #Nektar Units

#Add mass of a single ion
class FluidDat:
    def __init__(self,max_timesteps, initial_ion_density, initial_fluid_momentum):
        self.max_timesteps=max_timesteps #No units
        self.ion_density=np.zeros(max_timesteps+1)
        self.fluid_momentum=np.zeros(max_timesteps+1)
        self.ion_density[0]=initial_ion_density #Units SI
        self.fluid_momentum[0]=initial_fluid_momentum #Nektar Units

class Diagnostics:
    def __init__(self,fluid_temp,particle_temp):
        self.fluid_temp=fluid_temp
        self.particle_temp=particle_temp

#Add function which defines takes ParticleDat and FluidDat
#and perform charge exchange between them
        
IC = InitialConditions(0.0,0.0,0.0,0.5)
Neutrals = ParticleDat(max_timesteps,IC.neutrals_density,IC.neutrals_velocity)    
Fluid = FluidDat(max_timesteps,IC.fluid_density,IC.fluid_momentum)

while timestep <= max_timesteps:
    if timestep <= particle_injection_timestep_cutoff:
        neutrals_density[timestep]=timestep*max_neutrals_density/float(particle_injection_timestep_cutoff)
        ion_density[timestep]=max_neutrals_density
    else:
        neutrals_density[timestep]=neutrals_density[timestep-1]
        ion_density[timestep]=ion_density[timestep-1]
        
    if timestep<= particle_injection_timestep_cutoff and timestep!=0:   
        injected_number_particles=(max_neutrals_density/float(particle_injection_timestep_cutoff))
        dp_injected=mass_neutral*injected_number_particles*neutrals_velocity[timestep]*n_to_nektar
        total_momentum_injected[timestep]=total_momentum_injected[timestep-1] + dp_injected           
    if timestep>particle_injection_timestep_cutoff:
        total_momentum_injected[timestep]=total_momentum_injected[timestep-1]
        
    time[timestep]=timestep*dt  
    #Units dw
    #=[m^3 s^-1]*[unitless]*[m^-3]*[s]
    #=[Unitless] 
    # weight value changes as we increase the neutrals_density.
    w=dt_nektar*neutrals_density[timestep]
    
    dw = R_CE*w*ion_density[timestep]*dt
    if dw>w:
        dw=w
    
    #Units neutrals_momentum (nektar_units)
    #=[nektar_mass]*[nektar_velocity]*[SI_number_density]*[conversion_factor_SI_to_nektar_for_density]
    neutrals_momentum[timestep]=mass_neutral*neutrals_velocity[timestep]*neutrals_density[timestep]*n_to_nektar
    
    #Units fluid_velocity (nektar_units)
    #=[nektar_momentum]*[1/nektar_mass]*[1/SI_number_density]*[1/conversion_factor_SI_to_nektar_for_density]
    if ion_density[timestep] != 0.0:
        fluid_velocity[timestep]=fluid_momentum[timestep]/(mass_ion*ion_density[timestep]*n_to_nektar)   
    random_fluid_velocity=np.random.normal(fluid_velocity[timestep],fluid_velocity_sd,size=None)
    
    #Units dp_tot
    #=[unitless]*[nektar_mass]*[nektar_velocity] 
    dp_tot = dw*(- mass_neutral*neutrals_velocity[timestep] + mass_ion*random_fluid_velocity);
    
    if timestep < max_timesteps:
        if ion_density[timestep] !=0 and w != 0.0:
            neutrals_velocity[timestep+1]=neutrals_velocity[timestep] + dp_tot/(w*mass_neutral)
            fluid_momentum[timestep+1]=fluid_momentum[timestep] - dp_tot*n_to_nektar/dt_nektar
    timestep=timestep+1

plt.clf()
plt.plot(time, fluid_momentum,'b', label="Fluid Momentum")  
plt.plot(time, neutrals_momentum,'g', label="Neutrals Momentum")  

plt.plot(time, neutrals_momentum+fluid_momentum,'ro', label="Total Momentum of Fluid+Particles")  
plt.plot(time, total_momentum_injected,'y', label="Total injected momentum")  
plt.axvline(x = 100, color = 'k', linestyle='--')
plt.xlim([0, max_timesteps*dt])
plt.ylim(bottom=0)
plt.legend(loc="lower right")
plt.xlabel('Time [Seconds]')
plt.ylabel('Momentum [Nektar Units]')
plt.tight_layout()
plt.savefig('Momentum.png')
