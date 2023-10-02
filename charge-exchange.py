import numpy as np
import matplotlib.pyplot as plt  # Used for plotting
from scipy.linalg import expm  # Used for calculating analytical solution


class IonFluid:
    def __init__(self, mass: float, density: float, momentum: float) -> None:
        self.density = density
        self.mom = momentum
        self.mass = mass
        self.energy = (momentum ** 2.0) / (2.0 * mass)

    def getRandomMom(self) -> float:

        return np.random.normal(self.mom, self.energy)

    def getIonEnergy(self) -> float:

        return (self.mom ** 2.0) / (2.0 * self.mass)


class Particle:
    def __init__(self, mass: float, weight: float, vel: float) -> None:
        self.mass = mass
        self.weight = weight
        self.vel = vel

    def applyCX(self, dw, ionFluid: IonFluid) -> None:

        ionMom = ionFluid.getRandomMom()

        ionFluid.mom += dw * self.vel * self.mass / self.weight
        self.vel -= dw * self.vel / self.weight

        self.vel += dw * ionMom / (self.weight * self.mass)
        ionFluid.mom -= dw * ionMom / self.weight
        ionFluid.energy = ionFluid.getIonEnergy()


class runner:
    def __init__(
        self,
        ions: IonFluid,
        neutrals: Particle,
        number_of_timesteps,
        dt_nektar,
        dt_SI,
        CXrate,
    ):
        self.ions = ions
        self.neutrals = neutrals
        self.number_of_timesteps = number_of_timesteps
        self.dt_SI = dt_SI
        self.dt_nektar = dt_nektar
        self.CXrate = CXrate

    def runCX(self):
        neutralVel = np.zeros(self.number_of_timesteps)
        ionMom = np.zeros(self.number_of_timesteps)

        ionMom_analytical = np.zeros(self.number_of_timesteps)
        neutralVel_analytical = np.zeros(self.number_of_timesteps)
        time = np.zeros(self.number_of_timesteps)

        for i in range(len(ionMom)):
            # This section determines the numerical solution
            dw = self.dt_SI * self.CXrate * self.neutrals.weight * self.ions.density
            if dw > self.neutrals.weight:
                dw = self.neutrals.weight
            self.neutrals.applyCX(dw, self.ions)

            neutralVel[i] = self.neutrals.vel
            ionMom[i] = self.ions.mom
            time[i] = i * self.dt_nektar

            # This section determines the analytical solution
            a = np.matrix("-1 1; 1 -1")
            b = expm(
                (dw * i / self.neutrals.weight) * a
            )  # lambda = dw/(w*dt_nektar) time=i*dt_nektar => lambda*time=dt*i/weight
            c = np.array([initial_neutral_velocity, initial_fluid_momentum])
            d = b.dot(c)
            neutralVel_analytical[i] = d[0]
            ionMom_analytical[i] = d[1]

        plt.clf()
        plt.plot(time, ionMom, "b", label="Fluid Momentum")
        plt.plot(time, self.neutrals.mass * neutralVel, "g", label="Neutrals Momentum")
        plt.plot(
            time,
            self.neutrals.mass * neutralVel + ionMom,
            "r",
            label="Total Momentum of Fluid+Particles",
        )

        plt.plot(
            time,
            self.neutrals.mass * neutralVel_analytical,
            "y",
            label="Anlytical solution Particles",
        )
        plt.plot(time, ionMom_analytical, "k", label="Anlytical Solution Fluid")

        plt.xlim([0, len(ionMom) * dt_nektar])
        plt.ylim(bottom=0)
        plt.legend(loc="lower right")
        plt.xlabel("Time [Nektar units]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig("Momentum.png")


# Constants
# mass_ion=1.0
# mass_neutral=1.0

# Conversion of units
# SI_to_Nektar_n=3.33333e-19

# Initial Conditions
# weight=1.65E+16
# initial_neutral_velocity=0.2
# initial_fluid_momentum=0.8
# initial_fluid_density=3.3e18

# Timestep info
# dt_nektar=0.005
# dt_SI=1.66667e-07

# Rate info
# CXrate = 5E-14

# newPart = Particle(mass_neutral,weight,initial_neutral_velocity)
# newIonFluid = IonFluid(mass_ion,initial_fluid_density,initial_fluid_momentum)

# CXrunner = runner(newIonFluid,newPart,250,dt_nektar,dt_SI,CXrate)
# CXrunner.runCX()
