import numpy as np
import matplotlib.pyplot as plt  # Used for plotting
from scipy.stats import norm

class IonFluid:
    def __init__(
        self,
        mass: float,
        density: float,
        momentum: float,
        temperature: float,
        Volume: float,
    ) -> None:
        self.density = density  # Density of the fluid
        self.mom = momentum  # Bulk pseudo momentum (density * bulk velocity)
        self.mass = mass  # Mass of a single ion in the fluid
        # self.kinetic_energy = density*Volume*(0.5*temperature) #0.5*(Number_of_ions)*(average value of v^2 in terms of temperature)
        self.kinetic_energy = ((momentum * Volume) ** 2.0) / (
            2.0 * (density * Volume)
        ) + mass * density * Volume * ((temperature / mass)) / 2.0
        self.temperature = temperature
        self.Volume = Volume

    def getRandomMom(self, number_of_samples: int) -> float:
        # Sample velocity with standard deviation = particle thermal velocity = sqrt(KT/M)
        velocity_single_ion = self.mom / (self.mass * self.density)
        b = np.random.normal(
            velocity_single_ion,
            (self.temperature / self.mass) ** 0.5,
            number_of_samples,
        )
        b = (
            self.mass
            * (b - np.mean(b))
            * (((self.temperature / self.mass) ** 0.5) / (np.std(b)))
            + self.mass * velocity_single_ion
        )
        return b

    def updateKineticEnergy(self) -> None:
        self.kinetic_energy = ((self.mom * self.Volume) ** 2.0) / (
            2.0 * (self.mass * self.density * self.Volume)
        ) + self.mass * self.density * self.Volume * (
            (self.temperature / self.mass)
        ) / 2.0

    def updatetemperature(self, total_energy: float, neutrals_kinetic: float) -> None:
        self.temperature = (
            (
                total_energy
                - neutrals_kinetic
                - ((self.mom * self.Volume) ** 2.0)
                / (2.0 * (self.density * self.Volume))
            )
            * 2.0
        ) / (self.density * self.Volume)


class Plotting:
    def __init__(
        self,
        time,
        mass: float,
        weight: float,
        density: float,
        particle_velocity: float,
        fluid_momentum: float,
        Volume: float,
    ) -> None:
        self.time = time
        self.momentum_particles = np.zeros((len(time), len(weight[0,:])))
        self.mean = np.zeros(len(time))
        #for i in range(len(time)):
            #for j in range(len(particle_velocity[0, :])):
                #self.momentum_particles[i, j] = mass * weight[i,j] * particle_velocity[i, j]
        print("Calculating Momentum of Particles")
        for i in range(len(time)):
            self.momentum_particles[i, :] = mass * np.multiply(weight[i,:],particle_velocity[i, :])
        print("Calculating Bulk Momentum for plots")
        self.bulk_momentum = np.sum(self.momentum_particles, axis=1)

        #self.mean = mass * np.average(particle_velocity,weights=weight, axis=1)
        print("Calculating Mean Momentum of Particles for plots")
        total_weight = np.zeros(len(time))
        for i in range(len(time)):
            total_weight[i]=np.sum(weight[i,:])
            
        for i in range(len(time)):
            self.mean[i]=np.sum(self.momentum_particles[i, :])/total_weight[i]
            
        print("Done Calculating")
        #self.std = mass * np.std(particle_velocity, axis=1)
        self.fluid_momentum = fluid_momentum * Volume
        self.fluid_velocity = fluid_momentum / density

    def plot_bulk_properties(self, filename) -> None:
        plt.clf()
        plt.plot(self.time, self.bulk_momentum, "g", label="Total Momentum Particles")
        plt.plot(self.time, self.fluid_momentum, "k", label="Total Momentum Fluid")
        plt.plot(
            self.time,
            self.fluid_momentum + self.bulk_momentum,
            "r",
            label="Total Momentum Fluid + Particles",
        )
        plt.xlim([0, self.time.max()])
        plt.legend(loc="upper right")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig(filename)

    def plot_single_properties(self, filename) -> None:
        plt.clf()
        #plt.fill_between(
        #    self.time,
        #    self.mean - 3.0 * self.std,
        #    self.mean + 3.0 * self.std,
        #    color="g",
        #    alpha=0.2,
        #)
        #plt.fill_between(
        #    self.time,
        #    self.mean - 2.0 * self.std,
        #    self.mean + 2.0 * self.std,
        #    color="g",
        #    alpha=0.3,
        #)
        #plt.fill_between(
        #    self.time, self.mean - self.std, self.mean + self.std, color="g", alpha=0.4
        #)
        plt.plot(self.time, self.mean, "g", label="Mean Neutral Momentum")
        plt.plot(self.time, self.fluid_velocity, "k", label="Mean Ion Momentum")
        plt.xlim([0, self.time.max()])
        plt.legend(loc="upper right")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig(filename)

    def plot_log_exp(self, rate, ion_density, neutral_density, filename) -> None:
        plt.clf()
        
        slope_fluid=(self.fluid_velocity*(neutral_density+ion_density) - ion_density*self.fluid_velocity[0] - neutral_density*self.mean[0])/((self.fluid_velocity[0]-self.mean[0])*neutral_density)
        slope_fluid=np.log(slope_fluid)
        slope_particles=-(self.mean*(neutral_density+ion_density) - ion_density*self.fluid_velocity[0] - neutral_density*self.mean[0])/((self.fluid_velocity[0]-self.mean[0])*ion_density)
        slope_particles=np.log(slope_particles)

        
        plt.plot(self.time, slope_fluid, "ko", label="y= - lambda t (fluid)")
        plt.plot(self.time, slope_particles, "go", label="y= - lambda t (particles)")
        plt.plot(
            self.time,
            -rate * (ion_density+neutral_density) * self.time,
            "b",
            label="y= - R_CX*(n_ions + n_neutrals) t",
        )
        plt.xlabel("Time [Seconds]")
        plt.legend(loc="upper right")
        plt.savefig(filename)

    def plot_total_energy(
        self, particles_energy: float, fluid_energy: float, filename
    ) -> None:
        plt.clf()
        plt.plot(self.time, particles_energy, "r", label="Kinetic Energy Particles")
        plt.plot(self.time, fluid_energy, "g", label="Kinetic Energy Fluid")
        plt.plot(self.time, particles_energy + fluid_energy, "k", label="Total Energy")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Energy [Nektar Units]")
        plt.legend(loc="center right")
        plt.savefig(filename)

    def plot_temperature(
        self, fluid_temperature: float, filename
    ) -> None:
        plt.clf()
        plt.plot(self.time, fluid_temperature, "g", label="Fluid Temperature")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Temperature")
        plt.legend(loc="upper right")
        plt.savefig(filename)


class Particles:
    def __init__(self, mass: float, weight: float, vel: float,density: float) -> None:
        self.mass = mass
        self.weight = weight
        self.vel = vel
        self.kinetic_energy = 0.5 * np.sum(np.multiply(weight, np.multiply(vel, vel)))
        self.temperature = (
            2.0
            * self.mass
            * (self.kinetic_energy - 0.5 * np.sum(weight * np.mean(vel) * np.mean(vel)))
            / np.sum(weight)
        )
        self.density=density

    def applyCX(
        self,
        rate: float,
        dt: float,
        ionFluid: IonFluid,
        Volume: float
    ) -> None:
        n_samples=100
        vel_temp=np.copy(self.vel)
        weight_temp=np.copy(self.weight)
        single_momentum_random = ionFluid.getRandomMom(n_samples)
        random_integers_temp = np.arange(0,n_samples)
        ionSingleMom = ionFluid.getRandomMom(n_samples)
        length_before=len(self.vel)
        N_CX=dt*rate*ionFluid.density*self.weight
        self.vel=np.concatenate((self.vel, np.zeros(n_samples)), axis=None)
        self.weight=np.concatenate((self.weight, np.zeros(n_samples)), axis=None)
        k=np.floor_divide(len(vel_temp),100)
        
        random_integers=np.tile(random_integers_temp, k)
        self.weight[0:length_before]=weight_temp - N_CX
        self.vel[length_before:]=ionSingleMom
        ionFluid.mom += ionFluid.density * ionFluid.mass*dt*rate*np.sum(np.multiply(weight_temp, vel_temp))/Volume
        commulative_change = np.zeros((n_samples,k))
        commulative_change_2 = np.zeros((n_samples,k))
        
        for j in range(k):
            commulative_change[:,j]=N_CX[n_samples*j:n_samples*j +n_samples]
            commulative_change_2[:,j]=ionFluid.density * ionFluid.mass*dt*rate*np.multiply(weight_temp[n_samples*j:n_samples*j+n_samples], ionSingleMom[:]/ionFluid.mass)/Volume
        
        self.weight[length_before:]=np.sum(commulative_change,axis=1)
        ionFluid.mom -= np.sum(np.sum(commulative_change_2,axis=1),axis=0)

        kinetic_energy_before=self.kinetic_energy
        self.updateKineticEnergy()
        
        ionFluid.kinetic_energy -= self.kinetic_energy - kinetic_energy_before
        com_energy=0.5*ionFluid.density*Volume*((np.mean(ionFluid.mom/ionFluid.density))**2.0)
        thermal_energy=ionFluid.kinetic_energy-com_energy
        
        ionFluid.temperature=(2.0*thermal_energy/(ionFluid.density*Volume))
        

    def updateKineticEnergy(self) -> None:
        self.kinetic_energy = (
            (1.0 / 2.0)
            * self.mass
            * np.sum(np.multiply(self.weight, np.multiply(self.vel, self.vel)))
        )

    def updatetemperature(self) -> None:
        self.temperature = (
            2.0
            * self.mass
            * (
                self.kinetic_energy
                - 0.5 * np.sum(self.weight * np.mean(self.vel) * np.mean(self.vel))
            )
            / np.sum(self.weight)
        )


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
    def __init__(
        self, include_diagnostis: bool, include_sampling: bool, produce_plots: bool
    ) -> None:
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
        Volume,
    ):
        self.ions = ions
        self.neutrals = neutrals
        self.number_of_timesteps = number_of_timesteps
        self.dt_SI = dt_SI
        self.CXrate = CXrate
        self.Volume = Volume

    def runCX(self):
        neutralVel = np.zeros((self.number_of_timesteps, len(self.neutrals.weight)   +100*self.number_of_timesteps  ))
        neutralweight = np.zeros((self.number_of_timesteps, len(self.neutrals.weight)   +100*self.number_of_timesteps ))
        ionMom = np.zeros(self.number_of_timesteps)
        particles_energy = np.zeros(self.number_of_timesteps)
        fluid_energy = np.zeros(self.number_of_timesteps)
        temperature_fluid = np.zeros(self.number_of_timesteps)
        temperature_particles = np.zeros(self.number_of_timesteps)
        single_momentum_random = np.zeros(len(self.neutrals.weight)   +100*self.number_of_timesteps )
        time = np.zeros(self.number_of_timesteps)

        total_energy = self.ions.kinetic_energy + self.neutrals.kinetic_energy
        for i in range(len(ionMom)):
            ionMom[i] = self.ions.mom
            for j in range(len(self.neutrals.vel)):
                neutralVel[i, j] = self.neutrals.vel[j]
                neutralweight[i,j]= self.neutrals.weight[j]
            time[i] = i * self.dt_SI
            particles_energy[i] = self.neutrals.kinetic_energy
            fluid_energy[i] = self.ions.kinetic_energy
            temperature_fluid[i] = self.ions.temperature
            temperature_particles[i] = self.neutrals.temperature
            #plt.clf()
            #plt.title("time = "+str(round(time[i],10)))
            #count, bins  = np.histogram(single_momentum_random, 300)
            #sigma=(self.ions.temperature/self.ions.mass)**0.5
            #mu=self.ions.mom/self.ions.density
            #plt.plot(bins, 1/(sigma * np.sqrt(2.0 * np.pi)) * np.exp( - (bins - mu)**2.0 / (2 * sigma**2.0) ), linewidth=2, color='b', label="Fluid Distribution")
            #if i == 0:
                # lower_limit=np.min(self.neutrals.vel)
                # upper_limit=np.max(self.neutrals.vel)
            #    lower_limit = -3.5
            #    upper_limit = 5.5
            #plt.xlim([lower_limit, upper_limit])
            #plt.ylim(bottom=0)
            #count, bins  = np.histogram(self.neutrals.vel, 300)
            #sigma = np.std(self.neutrals.vel)
            #mu = np.mean(self.neutrals.vel)
            #plt.plot(
            #   bins,
            #   1
            #   / (sigma * np.sqrt(2.0 * np.pi))
            #   * np.exp(-((bins - mu) ** 2.0) / (2 * sigma**2.0)),
            #   linewidth=2,
            #   color="r",
            #   label="Particle Distribution"
            #)
            #plt.legend(loc="upper right")
            #if len(str(i)) == 1:
            #   plt.savefig("./Particle_Images/Distvel0000" + str(i) + ".png")
            #if len(str(i)) == 2:
            #   plt.savefig("./Particle_Images/Distvel000" + str(i) + ".png")
            #if len(str(i)) == 3:
            #   plt.savefig("./Particle_Images/Distvel00" + str(i) + ".png")
            #if len(str(i)) == 4:
            #   plt.savefig("./Particle_Images/Distvel0" + str(i) + ".png") 
            #if len(str(i)) == 5:
            #   plt.savefig("./Particle_Images/Distvel" + str(i) + ".png") 
            #plt.clf()
            #plt.title("time = "+str(round(time[i],10)))
            #plt.hist(self.neutrals.vel,bins=np.floor_divide(len(self.neutrals.vel),100),weights=self.neutrals.weight)
            #if len(str(i)) == 1:
            #   plt.savefig("./Particle_Images/Distvel0000" + str(i) + ".png")
            #if len(str(i)) == 2:
            #   plt.savefig("./Particle_Images/Distvel000" + str(i) + ".png")
            #if len(str(i)) == 3:
            #   plt.savefig("./Particle_Images/Distvel00" + str(i) + ".png")
            #if len(str(i)) == 4:
            #   plt.savefig("./Particle_Images/Distvel0" + str(i) + ".png") 
            #if len(str(i)) == 5:
            #   plt.savefig("./Particle_Images/Distvel" + str(i) + ".png") 
            self.neutrals.applyCX(
                self.CXrate, self.dt_SI, self.ions, self.Volume
            )
        
        #newPlot_numerical = Plotting(
        #    time,
        #    self.neutrals.mass,
        #    neutralweight,
        #    self.ions.density,
        #    neutralVel,
        #    ionMom,
        #    self.Volume,
        #)
        #newPlot_numerical.plot_bulk_properties("Momentum_Numerical_Bulk.png")
        #newPlot_numerical.plot_single_properties("Momentum_Numerical_Single.png")
        #newPlot_numerical.plot_log_exp(self.CXrate,self.ions.density,self.neutrals.density,"Momentum_Numerical_Slope.png")
        #newPlot_numerical.plot_total_energy(
        #    particles_energy, fluid_energy, "Energy.png"
        #)
        #newPlot_numerical.plot_temperature(
        #    temperature_fluid, "Temperature.png"
        #)
