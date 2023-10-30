import numpy as np
import matplotlib.pyplot as plt  # Used for plotting

class IonFluid:
    def __init__(self, mass: float, density: float, momentum: float, temperature: float, Volume: float) -> None:
        self.density = density #Density of the fluid
        self.mom = momentum #Bulk pseudo momentum (density * bulk velocity)
        self.mass = mass #Mass of a single ion in the fluid
        #self.kinetic_energy = density*Volume*(0.5*temperature) #0.5*(Number_of_ions)*(average value of v^2 in terms of temperature)
        self.kinetic_energy = ((momentum*Volume)**2.0)/(2.0*(density*Volume)) + mass*density*Volume*((temperature/mass))/2.0
        self.temperature= temperature
        self.Volume=Volume

    def getRandomMom(self, number_of_samples: int) -> float:
        #Sample velocity with standard deviation = particle thermal velocity = sqrt(KT/M) 
        velocity_single_ion=self.mom/(self.mass*self.density)
        b = np.random.normal(velocity_single_ion,(self.temperature/self.mass)**0.5,number_of_samples)
        b = self.mass*(b - np.mean(b))*(((self.temperature/self.mass)**0.5)/(np.std(b))) + self.mass*velocity_single_ion
        #plt.clf()
        #count, bins, ignored = plt.hist(b, 30, density=True)
        #sigma=(self.temperature/self.mass)**0.5
        #mu=velocity_single_ion
        #plt.plot(bins, 1/(sigma * np.sqrt(2.0 * np.pi)) * np.exp( - (bins - mu)**2.0 / (2 * sigma**2.0) ), linewidth=2, color='r')
        #plt.savefig("./Fluid_Images/Distvel"+str(n)+".png")
        return b

    def updateKineticEnergy(self) -> None:
        self.kinetic_energy = ((self.mom*self.Volume)**2.0)/(2.0*(self.mass*self.density*self.Volume)) + self.mass*self.density*self.Volume*((self.temperature/self.mass))/2.0
        
    def updatetemperature(self,total_energy: float,neutrals_kinetic: float) -> None:
        self.temperature = (((total_energy - neutrals_kinetic - ((self.mom*self.Volume)**2.0)/(2.0*(self.density*self.Volume)) )*2.0)/(self.density*self.Volume))
        
class Plotting:
    def __init__(self,time, mass: float, weight: float, density: float, particle_velocity: float,fluid_momentum: float, Volume: float) -> None:
        self.time = time
        self.momentum_particles = np.zeros((len(time), len(weight)))
        for i in range(len(time)):
            self.momentum_particles[i,:] = mass *weight[:]*particle_velocity[i,:]
            
        self.bulk_momentum =  np.sum(self.momentum_particles,axis=1)
        self.mean = mass * np.mean(particle_velocity,axis=1)
        self.std = mass * np.std(particle_velocity,axis=1)
        self.fluid_momentum = fluid_momentum*Volume
        self.fluid_velocity = fluid_momentum/density

    def plot_bulk_properties(self,filename) -> None:    
        plt.clf()
        plt.plot(self.time, self.bulk_momentum, "g", label="Total Momentum Particles")
        plt.plot(self.time, self.fluid_momentum , "k", label="Total Momentum Fluid")
        plt.plot(self.time, self.fluid_momentum + self.bulk_momentum , "r", label="Total Momentum Fluid + Particles")
        plt.xlim([0, self.time.max()])
        plt.legend(loc="upper right")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig(filename)    
    
    def plot_single_properties(self,filename) -> None:    
        plt.clf()
        plt.fill_between(self.time, self.mean - 3.0*self.std, self.mean + 3.0*self.std, color='g',alpha=0.2)
        plt.fill_between(self.time, self.mean - 2.0*self.std, self.mean + 2.0*self.std, color='g',alpha=0.3)
        plt.fill_between(self.time, self.mean - self.std, self.mean + self.std, color='g',alpha=0.4)
        plt.plot(self.time, self.mean, "g", label="Mean Neutral Momentum")
        plt.plot(self.time, self.fluid_velocity , "k", label="Mean Ion Momentum")
        plt.xlim([0, self.time.max()])
        plt.legend(loc="upper right")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig(filename)

    def plot_log_exp(self,rate,ion_density,filename) -> None:
        plt.clf()
        slope_fluid = np.log((self.fluid_velocity - self.fluid_velocity[0]/2.0 - self.mean[0]/2.0 )/(-self.mean[0]/2.0 + self.fluid_velocity[0]/2.0))
        slope_particles = np.log((self.mean - self.fluid_velocity[0]/2.0 - self.mean[0]/2.0 )/(self.mean[0]/2.0 - self.fluid_velocity[0]/2.0))
        plt.plot(self.time, slope_fluid ,"ko",label="y= - lambda t (fluid)")
        plt.plot(self.time, slope_particles ,"g",label="y= - lambda t (particles)")
        plt.plot(self.time, -2.0*rate*ion_density*self.time ,"b",label="y= - 2*R_CX*n_ions t")
        plt.xlabel("Time [Seconds]")
        plt.legend(loc="upper right")
        plt.savefig(filename)

    def plot_total_energy(self, particles_energy: float, fluid_energy: float, filename) -> None:
        plt.clf() 
        plt.plot(self.time,particles_energy,"r",label="Kinetic Energy Particles")
        plt.plot(self.time,fluid_energy,"g",label="Kinetic Energy Fluid")
        plt.plot(self.time,particles_energy + fluid_energy , "k", label="Total Energy")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Energy [Nektar Units]")
        plt.legend(loc="center right")
        plt.savefig(filename)

    def plot_temperature(self,fluid_temperature: float,particles_temperature: float, filename) -> None:
        plt.clf() 
        plt.plot(self.time,fluid_temperature,"g",label="Fluid Temperature")
        plt.plot(self.time,particles_temperature,"r",label="Particles Temperature")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Temperature")
        plt.legend(loc="center right")
        plt.savefig(filename)

class Particles:
    def __init__(self, mass: float, weight: float, vel: float) -> None:
        self.mass = mass
        self.weight = weight
        self.vel = vel
        self.kinetic_energy=0.5*np.sum(np.multiply(weight,np.multiply(vel,vel)))
        self.temperature=2.0*self.mass*(self.kinetic_energy - 0.5*np.sum(weight*np.mean(vel)*np.mean(vel)))/np.sum(weight)

    def applyCX(self, rate: float, dt: float, ionFluid: IonFluid, Volume: float, ionSingleMom: float) -> None:
        ionSingleMom=np.sort(ionSingleMom)
        self.vel=np.sort(self.vel)
        ionRandomMombulk = ionSingleMom*ionFluid.density
        #mean_vel_before=np.mean(self.vel)
        #vel_temp=np.copy(self.vel)
        ionFluid.mom += ionFluid.density*dt*rate*np.sum(np.multiply(self.weight,self.vel)) * self.mass / Volume
        ionFluid.mom -= dt*rate*np.sum(np.multiply(self.weight,ionRandomMombulk))/Volume
        self.vel += dt*rate*ionFluid.density*(ionSingleMom/self.mass - self.vel)
        #ionSingleMom += dt*rate*ionFluid.density*(self.mass*self.vel - ionSingleMom)
        #self.vel = vel_temp + np.full(len(self.vel),np.mean(self.vel - vel_temp))

    def updateKineticEnergy(self) -> None:
        self.kinetic_energy = (1.0/2.0)*self.mass*np.sum(np.multiply(self.weight,np.multiply(self.vel,self.vel)))

    def updatetemperature(self) -> None:
        self.temperature = 2.0*self.mass*(self.kinetic_energy - 0.5*np.sum(self.weight*np.mean(self.vel)*np.mean(self.vel)))/np.sum(self.weight)

class InitialConditions:
    def __init__(
        self,
        fluid_density: float,
        fluid_momentum: float,
        weight: float,
        neutrals_velocity: float,
    ) -> None:
        self.fluid_density = fluid_density
        self.fluid_momentum = fluid_momentum
        self.weight = weight
        self.neutrals_velocity = neutrals_velocity

class Diagnostics:
    def __init__(self, include_diagnostis: bool, include_sampling: bool, produce_plots: bool) -> None:
        self.include_diagnostics = include_diagnostics
        self.include_sampling = include_sampling
        self.produce_plots = produce_plots

class runner:
    def __init__(
        self,
        ions: IonFluid,
        neutrals: Particles,
        number_of_timesteps,
        dt_SI,
        CXrate,
        Volume
    ):
        self.ions = ions
        self.neutrals = neutrals
        self.number_of_timesteps = number_of_timesteps
        self.dt_SI = dt_SI
        self.CXrate = CXrate
        self.Volume=Volume

    def runCX(self):
        neutralVel = np.zeros((self.number_of_timesteps, len(self.neutrals.weight)))
        ionMom = np.zeros(self.number_of_timesteps)
        particles_energy = np.zeros(self.number_of_timesteps)
        fluid_energy = np.zeros(self.number_of_timesteps)
        temperature_fluid = np.zeros(self.number_of_timesteps)
        temperature_particles = np.zeros(self.number_of_timesteps)
        single_momentum_random = np.zeros(len(self.neutrals.weight))
        time = np.zeros(self.number_of_timesteps)
        
        total_energy = self.ions.kinetic_energy + self.neutrals.kinetic_energy
        for i in range(len(ionMom)):
            ionMom[i] = self.ions.mom
            neutralVel[i, :] = self.neutrals.vel[:]
            time[i] = i * self.dt_SI
            particles_energy[i]= self.neutrals.kinetic_energy
            fluid_energy[i] = self.ions.kinetic_energy
            temperature_fluid[i]=self.ions.temperature
            temperature_particles[i]=self.neutrals.temperature
            single_momentum_random=self.ions.getRandomMom(len(self.neutrals.weight))
            self.neutrals.applyCX(self.CXrate,self.dt_SI, self.ions, self.Volume,single_momentum_random)
            plt.clf()
            if i==0:
                #lower_limit=np.min(self.neutrals.vel)
                #upper_limit=np.max(self.neutrals.vel)
                lower_limit=-2.5
                upper_limit=5.5
            plt.xlim([lower_limit, upper_limit])
            count, bins, ignored = plt.hist(self.neutrals.vel, 300, density=True)
            sigma=np.std(self.neutrals.vel)
            mu=np.mean(self.neutrals.vel)
            plt.plot(bins, 1/(sigma * np.sqrt(2.0 * np.pi)) * np.exp( - (bins - mu)**2.0 / (2 * sigma**2.0) ), linewidth=2, color='r')
            if len(str(i))==1:
                plt.savefig("./Particle_Images/Distvel000"+str(i)+".png")
            if len(str(i))==2:
                plt.savefig("./Particle_Images/Distvel00"+str(i)+".png")
            if len(str(i))==3:
                plt.savefig("./Particle_Images/Distvel0"+str(i)+".png")
            if len(str(i))==4:
                plt.savefig("./Particle_Images/Distvel"+str(i)+".png")
            self.neutrals.updateKineticEnergy() #Updates Neutrals Kinetic Energy
            self.neutrals.updatetemperature() #Updates Temperature of neutrals
            self.ions.updatetemperature(total_energy,self.neutrals.kinetic_energy) #Update fluid temperature
            self.ions.updateKineticEnergy() #Updates Ions Kinetic Energy
            
        newPlot_numerical = Plotting(time,self.neutrals.mass,self.neutrals.weight,self.ions.density,neutralVel,ionMom, self.Volume)
        newPlot_numerical.plot_bulk_properties("Momentum_Numerical_Bulk.png")
        newPlot_numerical.plot_single_properties("Momentum_Numerical_Single.png")
        #newPlot_numerical.plot_log_exp(self.CXrate,self.ions.density,"Momentum_Numerical_Slope.png")
        newPlot_numerical.plot_total_energy(particles_energy,fluid_energy,"Energy.png")
        newPlot_numerical.plot_temperature(temperature_fluid,temperature_particles,"Temperature.png")
