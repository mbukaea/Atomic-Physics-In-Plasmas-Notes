import chargeexchange as CX
import unittest

class Testing(unittest.TestCase):
    # This test checks the initialisation of the IonFluid class and functions therein
    def test_IonFluid_class(self):
        print("\n Testing the IonFluid class")
        # Constants
        mass_ion = 1.0

        # Initial Conditions
        initial_fluid_momentum = 0.8
        initial_fluid_density = 3.3e18
        
        newIonFluid = CX.IonFluid(
            mass_ion, initial_fluid_density, initial_fluid_momentum
        )

        self.assertEqual(newIonFluid.mom,initial_fluid_momentum)
        self.assertEqual(newIonFluid.mass,mass_ion)
        self.assertEqual(newIonFluid.density,initial_fluid_density)
        self.assertEqual(newIonFluid.getIonEnergy(),(initial_fluid_momentum**2.0)/(2.0*mass_ion))
        self.assertEqual(newIonFluid.energy,(initial_fluid_momentum**2.0)/(2.0*mass_ion))
        
    # This test checks the initialisation of the Particle class  
    def test_Particle_class(self):
        print("\n Testing the Particle class")
        # Constants
        mass_neutral = 1.0

        # Initial Conditions
        weight = 1.65e16
        initial_neutral_velocity = 0.2

        newPart = CX.Particle(mass_neutral, weight, initial_neutral_velocity)

        self.assertEqual(newPart.mass,mass_neutral)
        self.assertEqual(newPart.vel,initial_neutral_velocity)
        self.assertEqual(newPart.weight,weight)

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
        
        self.assertEqual(IC.fluid_density,initial_fluid_density)
        self.assertEqual(IC.fluid_momentum,initial_fluid_momentum)
        self.assertEqual(IC.weight,weight)
        self.assertEqual(IC.neutrals_velocity,initial_neutral_velocity)
    
    # This test checks whether the charge exchange function conserves momentum
    # after running for 250 timesteps
    def test_CX_function(self):
        print("\n Testing the charge exchange function contained within particle class")
        # Constants
        mass_ion = 1.0
        mass_neutral = 1.0

        # Conversion of units
        SI_to_Nektar_n = 3.33333e-19

        # Initial Conditions
        weight = 1.65e16
        initial_neutral_velocity = 0.2
        initial_fluid_momentum = 0.8
        initial_fluid_density = 3.3e18

        # Timestep info
        number_of_timesteps = 250
        dt_nektar = 0.005
        dt_SI = 1.66667e-07

        # Rate info
        CXrate = 5e-14

        IC = CX.InitialConditions(
            initial_fluid_density,
            initial_fluid_momentum,
            weight,
            initial_neutral_velocity,
        )
        newPart = CX.Particle(mass_neutral, weight, initial_neutral_velocity)
        newIonFluid = CX.IonFluid(
            mass_ion, initial_fluid_density, initial_fluid_momentum
        )

        CXrunner = CX.runner(IC, newIonFluid, newPart, number_of_timesteps, dt_nektar, dt_SI, CXrate)
        CXrunner.runCX()
        self.assertLess(
            abs(
                CXrunner.ions.mom
                + CXrunner.neutrals.mass * CXrunner.neutrals.vel
                - initial_fluid_momentum
                - mass_neutral * initial_neutral_velocity
            ),
            1e-12,
        )


if __name__ == "__main__":
    unittest.main()
