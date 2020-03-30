import unittest
from TWINGL.TwinglModel import DirectEvolutionModel as Twigl
import numpy as np
from PyCAR.PyCIT.lector_mest_posta import LoadCDP, get_reactivity_pro as grp, Geometry
from PyCAR.PyCIT.FT import MeshFlux
from glob import glob
import os
import subprocess

BLACKABSORBER = 10
DEFAULT_STATE = 0
DEFAULT_BURNUP = 0


class TWIGL_MAP(object):
    def MapMaterial(self, meshmaterial):
        attr = '_reg' + str(meshmaterial)
        try:
            return getattr(self, attr)[0][0.0]
        except AttributeError:
            return getattr(self, '_NULL')


class Nu_Fission_Map(TWIGL_MAP):
    wdir = 'C:\\TWINGL\\'

    def __init__(self, **kwargs):
        if 'wdir' in kwargs:
            setattr(self, 'wdir', kwargs['wdir'])

        wdir = getattr(self, 'wdir')
        self._reg1 = LoadCDP(wdir + 'REGION1.cdp')
        self._reg2 = LoadCDP(wdir + 'REGION2.cdp')
        self._reg3 = LoadCDP(wdir + 'REGION3.cdp')

        self._NULL = {'XS': np.zeros(self._reg3[0][0.0]['XS'].shape),
                      'SM': np.zeros(self._reg3[0][0.0]['SM'].shape)}


class Kinetic_Map(TWIGL_MAP):
    def __init__(self):
        v1 = 10 ** 7
        v2 = 2 * (10 ** 5)
        self._reg1 = self._reg2 = self._reg3 = {0: {0.0: {'Beta': np.array([0.0075]),
                                                          'Decays': np.array([0.08]),
                                                          'Neutron Velocity': np.array([1 / v1, 1 / v2])}}}

        self._NULL = {'Beta': [0],
                      'Decays': [0],
                      'Neutron Velocity': [0, 0]}


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

    def test_eq_precrs(self):
        self.MODEL.equilibrium_precrs()

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
#         Twigl(file=MainFile,
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
