import unittest
from TWINGL.TwinglModel import DirectEvolutionModel as Twigl, Kinetic_Map, Nu_Fission_Map
from glob import glob
import os
import subprocess


class TestTwingl(unittest.TestCase):
    def test_init_maps(self):
        Kinetic_Map()
        Nu_Fission_Map()
        return 0

    os.chdir('C:\\TestD\\')
    MainFile = 'twiglS.cii'
    MODEL = Twigl(file=MainFile,
                  KineticMap=Kinetic_Map,
                  NuFissionMap=Nu_Fission_Map,
                  groups=2,
                  UseDerivative=True,
                  WorkingDirectory='C:\\TestD\\',
                  seccion_4=MainFile.replace('S.cii', '.sc4'),
                  seccion_5=MainFile.replace('S.cii', '.sc5'),
                  seccion_26=MainFile.replace('S.cii', '.sc26'),
                  materiales=MainFile.replace('S.cii', '.mat'),
                  geo_type='XY',
                  black_absorber_ID=10)

    def bring_in_files(self):
        if 'C:\\TestD' not in glob('C:\\*'):
            subprocess.run(['mkdir', 'C:\\TestD'], shell=True)

        RootDirectory = 'C:\\TWINGL\\'
        FilesToCopy = ['twiglS.cii', 'twigl.sc4', 'twigl.sc5',
                       'twigl.sc26', 'twigl.mat', 'twigl0.ci@',
                       'twigl_eq.cdb.meshflux', 'twigl.bib']

        DirectoriesToCopy = ['XSU_MOD', 'Source_To_LU17']

        for file in FilesToCopy:
            subprocess.run(['copy', '/Y', RootDirectory + file, '.'], shell=True)

        for Directory in DirectoriesToCopy:
            if Directory not in glob('*'):
                subprocess.run(['mkdir', Directory], shell=True)
            subprocess.run(['copy', '/Y', RootDirectory + Directory, Directory], shell=True)
        return

    def test_eq_precrs(self):
        self.MODEL.equilibrium_precrs()
        return

    def test_evolve(self):
        dt = 0.01
        self.MODEL.evolve(dt)
        print(self.MODEL.t, self.MODEL.Pt)
        return


if __name__ == '__main__':
    unittest.main()

# for N in np.logspace(3, 5, dtype=int, num=10)[2:]:  # N = 4641
#     dt = tfinal / N
#     Times = []
#     Powers = []
#     for n in range(N + 1):
#         TwiglModel.evolve(dt)
#         print(TwiglModel.t, TwiglModel.Pt)
#         Powers.append(TwiglModel.Pt)
#         Times.append(TwiglModel.t)
#
#     np.savetxt('evol_{}_{}.txt'.format(N, UseDerivative),
#                np.r_['0', np.array([Times]), np.array([Powers])].T)

# def test_group_mismatch(self):
#     self.bring_in_files()
#     MainFile = 'twiglS.cii'
#     NumberOfGroups = 1 # debería ser 2
#     with self.assertRaises(AssertionError):
#         Twigl(FluxFile=MainFile,
#               KineticMap=Kinetic_Map,
#               NuFissionMap=Nu_Fission_Map,
#               groups=NumberOfGroups,
#               UseDerivative=True,
#               WorkingDirectory='C:\\TestD',
#               seccion_4=MainFile.replace('S.cii', '.sc4'),
#               seccion_5=MainFile.replace('S.cii', '.sc5'),
#               seccion_26=MainFile.replace('S.cii', '.sc26'),
#               materiales=MainFile.replace('S.cii', '.mat'),
#               geo_type='XY')
#     return 0
