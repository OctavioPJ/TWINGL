from TWINGL.TwinglModel import DirectEvolutionModel as Twigl, Kinetic_Map, Nu_Fission_Map
from CAREM.Citvap import CitvapModel
from CAREM.Poinki import PointKineticsModel
from PyCAR.PyCIT.FT import MeshFlux
import os
import numpy as np
import subprocess


class Cross(CitvapModel):
    def __init__(self, ModelFile, CitvapFile, *args, **kwargs):
        super().__init__(file=CitvapFile, *args, **kwargs)
        self.Model = Twigl(file=ModelFile, *args, **kwargs)

    pass  # Cross


if __name__ == '__main__':
    os.chdir('C:\\TestD\\')
    MainFile = 'twiglS.cii'
    MODEL = Cross(ModelFile=MainFile,
                  CitvapFile=MainFile.replace('S.cii', '.cii'),
                  KineticMap=Kinetic_Map,
                  NuFissionMap=Nu_Fission_Map,
                  groups=2,
                  UseDerivative=True,
                  WorkingDirectory='C:\\TestD\\',
                  seccion_4=MainFile.replace('S.cii', '.sc4'),
                  seccion_5=MainFile.replace('S.cii', '.sc5'),
                  seccion_26=MainFile.replace('S.cii', '.sc26'),
                  materiales=MainFile.replace('S.cii', '.mat'),
                  GeometricType='XY',
                  black_absorber_ID=10,
                  FluxConvergence=1.0E-6)

    MODEL.Calculate()

    tfinal = 0.5
    dt = 0.001
    Nsteps = int(tfinal/dt)
    Powers = []
    Times = []
    Reactivity = []
    for n in range(Nsteps+1):
        MODEL.Model.evolve(dt)

        MODEL.run(executable='pre_cit')
        MODEL.Model._xsu_mod__(NoDerivative=True, OutsideFile='twigl.XSU')
        MODEL.run(executable='citvap.exe')
        subprocess.run(['copy', '/Y', 'twigl.cdb', 'DBad\\time{}.cdb'.format(MODEL.Model.t)], shell=True)

        Powers.append(MODEL.Model.Pt)
        Times.append(MODEL.Model.t)
        Reactivity.append(MODEL.React[-1]+750)
        print(MODEL.Model.t, MODEL.Model.Pt, MODEL.React[-1]+750)

        subprocess.run(['copy', '/Y', 'twiglS.cdb', 'DBRes1\\time{}.cdb'.format(MODEL.Model.t)], shell=True)
        subprocess.run(['copy', '/Y', 'twiglS.cio', 'DBRes1\\time{}.cio'.format(MODEL.Model.t)], shell=True)

    np.savetxt('rampa_react.txt', np.array([Times, Powers]))
