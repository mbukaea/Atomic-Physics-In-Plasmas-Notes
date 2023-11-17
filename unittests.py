import numpy as np
import chargeexchange as CX
import unittest

class Testing(unittest.TestCase):
    # This test checks the initialisation of the IonFluid class and functions therein
    def test_IonFluid_class(self):
        print("\n Testing the IonFluid class")
        # Constants
        mass_ion = 1.0

        # Initial Conditions
        initial_fluid_momentum_single_ion = 0.8
        initial_fluid_density = 3.3e18
        initial_fluid_temperature=1.0
        initial_fluid_bulk_momentum=initial_fluid_momentum_single_ion*initial_fluid_density
        volume=1

        newIonFluid = CX.IonFluid(
            mass_ion, initial_fluid_density, initial_fluid_bulk_momentum,initial_fluid_temperature,volume
        )

        self.assertEqual(newIonFluid.mom, initial_fluid_bulk_momentum)
        self.assertEqual(newIonFluid.mass, mass_ion)
        self.assertEqual(newIonFluid.density, initial_fluid_density)
        #self.assertEqual(
        #    newIonFluid.getIonEnergy(),
        #    (initial_fluid_bulk_momentum**2.0) / (2.0 * mass_ion*initial_fluid_density),
        #)
        self.assertEqual(
            newIonFluid.kinetic_energy, ((initial_fluid_bulk_momentum*volume)**2.0)/(2.0*(initial_fluid_density*volume)) + mass_ion*initial_fluid_density*volume*((initial_fluid_temperature/mass_ion))/2.0
        )

    # This test checks the initialisation of the Particle class
    def test_Particle_class(self):
        print("\n Testing the Particle class")
        # Constants
        mass_neutral = 1.0

        # Initial Conditions
        weight = np.array([1.65e16, 2.0e16])
        initial_neutral_velocity = np.array([0.2, 0.3])

        newPart = CX.Particles(mass_neutral, weight, initial_neutral_velocity)

        self.assertEqual(newPart.mass, mass_neutral)
        self.assertEqual(newPart.vel[0], initial_neutral_velocity[0])
        self.assertEqual(newPart.vel[1], initial_neutral_velocity[1])
        self.assertEqual(newPart.weight[0], weight[0])
        self.assertEqual(newPart.weight[1], weight[1])

    # This test checks the initialisation of the InitialConditions class
    def test_InitialConditions_class(self):
        print("\n Testing the InitialConditions class")
        # Initial Conditions
        weight = 1.65e16
        initial_neutral_velocity = 0.2
        initial_fluid_momentum = 0.8
        initial_fluid_density = 3.3e18

        IC = CX.InitialConditions(
            initial_fluid_density,
            initial_fluid_momentum,
            weight,
            initial_neutral_velocity,
        )

        self.assertEqual(IC.fluid_density, initial_fluid_density)
        self.assertEqual(IC.fluid_momentum, initial_fluid_momentum)
        self.assertEqual(IC.weight, weight)
        self.assertEqual(IC.neutrals_velocity, initial_neutral_velocity)

    # This test checks whether the charge exchange function conserves total momentum
    def test_CX_function(self):
        print("\n Testing the charge exchange function contained within particle class")
        # Constants
        mass_ion = 1.0
        mass_neutral = 1.0
        number_of_macroparticles=200000
        volume=2.0
        # Conversion of units (currently not used)
        SI_to_Nektar_n = 3.33333e-19

        #Fluid properties
        initial_fluid_momentum_single_ion = 1.5
        initial_fluid_density = 3.3e18
        initial_fluid_temperature= 1E-12
        initial_fluid_bulk_momentum_density=initial_fluid_momentum_single_ion*initial_fluid_density

        #Particle properties
        neutrals_density=3.3e+18   
        mu = 1.5
        particle_thermal_velocity = 1.0 # mean and standard deviation
        initial_neutral_velocity = np.random.normal(mu, particle_thermal_velocity, number_of_macroparticles)
        initial_neutral_velocity = (initial_neutral_velocity - np.mean(initial_neutral_velocity))*(particle_thermal_velocity/np.std(initial_neutral_velocity)) + mu
        weight = np.full(number_of_macroparticles,neutrals_density*volume/float(number_of_macroparticles))
        #Timestep info
        number_of_timesteps = 32000
        dt_SI = 1e-9

        #Rate info
        CXrate = 5e-14

        #Define Initial conditions class
        IC = CX.InitialConditions(
            initial_fluid_density,
            initial_fluid_bulk_momentum_density,
            weight,
            initial_neutral_velocity,
        )

        #Define particles and fluid
        newPart = CX.Particles(mass_neutral, weight, initial_neutral_velocity)
        newIonFluid = CX.IonFluid(
            mass_ion, initial_fluid_density, initial_fluid_bulk_momentum_density,initial_fluid_temperature,volume
        )
        #Setup runner class
        CXrunner = CX.runner( newIonFluid, newPart, number_of_timesteps, dt_SI, CXrate,volume)

        momentum_particles_before = mass_neutral*np.sum(np.multiply(weight,initial_neutral_velocity))
        momentum_fluid_before = initial_fluid_bulk_momentum_density*volume
        
        #Peform charge exchange between fluid and particles
        CXrunner.runCX()
  
        momentum_particles_after = CXrunner.neutrals.mass*np.sum(np.multiply(CXrunner.neutrals.weight,CXrunner.neutrals.vel))
        
        #Test that total momentum is conserved
        self.assertLess(
        abs(
              CXrunner.ions.mom*volume
            + momentum_particles_after
            - momentum_fluid_before
            - momentum_particles_before
        )/abs(initial_fluid_bulk_momentum_density*volume+momentum_particles_before),
        1e-14, "Total Momentum not conserved for charge exchange"
        )


if __name__ == "__main__":
    unittest.main()
