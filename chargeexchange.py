import numpy as np
import matplotlib.pyplot as plt  # Used for plotting
from scipy.linalg import expm  # Used for calculating analytical solution


class IonFluid:
    def __init__(self, mass: float, density: float, momentum: float) -> None:
        self.density = density
        self.mom = momentum
        self.mass = mass
        self.energy = self.getIonEnergy()

    def getRandomMom(self) -> float:
        return np.random.normal(self.mom, self.energy)

    def getIonEnergy(self) -> float:
        return (self.mom**2.0) / (2.0 * self.mass)


class Particles:
    def __init__(self, mass: float, weight: float, vel: float) -> None:
        self.mass = mass
        self.weight = weight
        self.vel = vel

    def applyCX(self, dw, ionFluid: IonFluid, k: int) -> None:
        ionMom = ionFluid.getRandomMom()

        ionFluid.mom += dw * self.vel[k] * self.mass / self.weight[k]
        self.vel[k] -= dw * self.vel[k] / self.weight[k]

        self.vel[k] += dw * ionMom / (self.weight[k] * self.mass)
        ionFluid.mom -= dw * ionMom / self.weight[k]

        ionFluid.energy = ionFluid.getIonEnergy()


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
        initial_conditions: InitialConditions,
        ions: IonFluid,
        neutrals: Particles,
        number_of_timesteps,
        dt_SI,
        CXrate,
    ):
        self.initial_conditions = initial_conditions
        self.ions = ions
        self.neutrals = neutrals
        self.number_of_timesteps = number_of_timesteps
        self.dt_SI = dt_SI
        self.CXrate = CXrate

    def runCX(self):
        neutralVel = np.zeros((self.number_of_timesteps, len(self.neutrals.weight)))
        ionMom = np.zeros(self.number_of_timesteps)
        ionMom_analytical = np.zeros(self.number_of_timesteps)
        neutralVel_analytical = np.zeros(
            (self.number_of_timesteps, len(self.neutrals.weight))
        )
        time = np.zeros(self.number_of_timesteps)
        dw = np.zeros(len(self.neutrals.weight))
        a = np.zeros((len(self.neutrals.weight) + 1, len(self.neutrals.weight) + 1))
        for j in range(len(self.neutrals.weight) + 1):
            a[j, j] = -1
            a[j, len(self.neutrals.weight)] = 1
            if j == len(self.neutrals.weight):
                a[j, 0:j] = 1
                a[j, j] = -len(self.neutrals.weight)
        v0 = np.concatenate(
            (
                self.initial_conditions.neutrals_velocity,
                np.array([self.initial_conditions.fluid_momentum]),
            ),
            axis=None,
        )
        for i in range(len(ionMom)):
            dmom = 0.0
            ionMom[i] = self.ions.mom
            neutralVel[i, :] = self.neutrals.vel[:]
            time[i] = i * self.dt_SI
            # This section determines the analytical solution
            b = expm(
                (self.CXrate * self.ions.density * time[i]) * a
            )  # lambda = dw/(w*dt_nektar) time=i*dt_nektar => lambda*time=dt*i/weight
            v = b.dot(v0)
            neutralVel_analytical[i, :] = v[: (len(self.neutrals.weight))]
            ionMom_analytical[i] = v[len(self.neutrals.weight)]
            for j in range(len(self.neutrals.weight)):
                # This section determines the numerical solution
                dw[j] = (
                    self.dt_SI
                    * self.CXrate
                    * self.neutrals.weight[j]
                    * self.ions.density
                )
                if dw[j] > self.neutrals.weight[j]:
                    dw[j] = self.neutrals.weight[j]
                self.neutrals.applyCX(dw[j], self.ions, j)
                dmom += self.ions.mom - ionMom[i]
                self.ions.mom = ionMom[i]
                self.ions.getIonEnergy()
                if j == len(self.neutrals.weight) - 1:
                    self.ions.mom = ionMom[i] + dmom
                    self.ions.getIonEnergy()

        plt.clf()
        plt.plot(
            time,
            self.neutrals.mass * neutralVel_analytical[:, 0],
            "g",
            label="Anlytical Solution Macroparticles",
        )
        for j in range(1, len(dw)):
            plt.plot(time, self.neutrals.mass * neutralVel_analytical[:, j], "g")
        plt.plot(time, ionMom_analytical, "k", label="Anlytical Solution Fluid")
        plt.xlim([0, len(ionMom) * self.dt_SI])
        plt.ylim(bottom=0)
        plt.legend(loc="lower right")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig("MomentumAnalytical.png")

        plt.clf()
        plt.plot(
            time,
            self.neutrals.mass * neutralVel[:, 0],
            "g",
            label="Numerical Solution Macroparticles Momentum",
        )
        for j in range(1, len(dw)):
            plt.plot(time, self.neutrals.mass * neutralVel[:, j], "g")
        plt.plot(time, ionMom, "b", label="Numerical Solution Fluid Momentum")
        plt.xlim([0, len(ionMom) * self.dt_SI])
        plt.ylim(bottom=0)
        plt.legend(loc="lower right")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig("MomentumNumerical.png")

        plt.clf()
        plt.plot(
            time,
            self.neutrals.mass * neutralVel_analytical[:, 0]
            - self.neutrals.mass * neutralVel[:, 0],
            "g",
            label="Difference Macroparticle Momentum",
        )
        for j in range(1, len(dw)):
            plt.plot(
                time,
                self.neutrals.mass * neutralVel_analytical[:, j]
                - self.neutrals.mass * neutralVel[:, j],
                "g",
            )
        plt.plot(
            time, ionMom_analytical - ionMom, "b", label="Difference Fluid Momentum"
        )
        plt.xlim([0, len(ionMom) * self.dt_SI])
        plt.ylim(bottom=0)
        plt.legend(loc="lower right")
        plt.xlabel("Time [Seconds]")
        plt.ylabel("Momentum [Nektar Units]")
        plt.tight_layout()
        plt.savefig("MomentumDifference.png")
