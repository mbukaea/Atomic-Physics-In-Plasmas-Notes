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


class InitialConditions:
    def __init__(
        self, fluid_density, fluid_momentum, weight, neutrals_velocity
    ) -> None:
        self.fluid_density = fluid_density
        self.fluid_momentum = fluid_momentum
        self.weight = weight
        self.neutrals_velocity = neutrals_velocity


class Diagnostics:
    def __init__(self, include_diagnostis: bool) -> None:
        self.include_diagnostics = include_diagnostics


class runner:
    def __init__(
        self,
        initial_conditions: InitialConditions,
        ions: IonFluid,
        neutrals: Particle,
        number_of_timesteps,
        dt_nektar,
        dt_SI,
        CXrate,
    ):
        self.initial_conditions = initial_conditions
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
            a = np.array([[-1, 1], [1, -1]])
            b = expm(
                (dw * i / self.initial_conditions.weight) * a
            )  # lambda = dw/(w*dt_nektar) time=i*dt_nektar => lambda*time=dt*i/weight
            c = np.array(
                [
                    self.initial_conditions.neutrals_velocity,
                    self.initial_conditions.fluid_momentum,
                ]
            )
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

        plt.xlim([0, len(ionMom) * self.dt_nektar])
        plt.ylim(bottom=0)
        plt.legend(loc="lower right")
        plt.xlabel("Time [Nektar units]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig("Momentum.png")
