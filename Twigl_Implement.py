from TWINGL.TwinglModel import DirectEvolutionModel as Twigl, Kinetic_Map, Nu_Fission_Map
import os
import numpy as np
import subprocess

if __name__ == '__main__':
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
                  black_absorber_ID=10,
                  FluxConvergence=1.0E-6)

    tfinal = 0.5
    dt = 0.001
    Nsteps = int(tfinal/dt)
    Powers = []
    Times = []

    for n in range(Nsteps+1):
        MODEL.evolve(dt)
        Powers.append(MODEL.Pt)
        Times.append(MODEL.t)
        print(MODEL.t, MODEL.Pt)
        subprocess.run(['copy', '/Y', 'twiglS.cdb', 'DBRes\\time{}.cdb'.format(MODEL.t)], shell=True)
        subprocess.run(['copy', '/Y', 'twiglS.cio', 'DBRes\\time{}.cio'.format(MODEL.t)], shell=True)

        MODEL._xsu_mod__()
    np.savetxt('rampa_clase.txt', np.array([Times, Powers]))
