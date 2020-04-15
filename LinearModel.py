import numpy as np
from CAREM.Poinki import PointLinearModel as LinearModel, PointKineticsModel
from PyCAR.PyCIT.lector_mest_posta import get_reactivity_pro as grp
import os
from functools import reduce
from operator import add
import matplotlib.pyplot as plt
os.chdir('C:\\TestD')


if __name__ == '__main__':
    alf = LinearModel(beta_info=0.0075, decay=0.08, Lambda=31.9157, nprc=1)
    realf = np.array(reduce(add, grp('twigl.PLT', last=False)))
    k = 1/(1-realf/100000)/(1-0.0075)
    Ro = 1 - 1/k
    yp = Ro[-1]*1E5/0.2

    notalf = PointKineticsModel(beta_info=0.0075, decay=0.08, Lambda=31.9157, nprc=1)

    dt = 0.001
    tfinal = 0.5
    Powers = []
    Times = []
    Nsteps = int(tfinal / dt)
    t = 0
    EqFLAG = True

    for n in range(Nsteps + 1):
        notalf.evolve(2*dt, dt, Ro[n]*1E5, Equilibrium=EqFLAG, P0=1.0)
        Powers.append(notalf.Ut[-1, 0])
        t += dt
        Times.append(t)
        if EqFLAG:
            EqFLAG = False

    P1 = alf.evolve(0.2+dt, dt, 0, yp, P0=1.0)[:, 0]
    t1 = alf.t
    P2 = alf.evolve(0.3+dt, dt, Ro[-1]*1E5, 0)[:, 0]
    t2 = alf.t + t1[-1]

    DATA_RAMP = np.array([[0, 1, 1],
                          [0.1, 1.311, 1.308],
                          [0.2, 1.952, 1.961],
                          [0.3, 2.045, 2.076],
                          [0.4, 2.069, 2.093],
                          [0.5, 2.092, 2.111]])

    plt.plot(np.r_['0', t1, t2], np.r_['0', P1, P2], '-o',
             DATA_RAMP[:, 0], DATA_RAMP[:, 1], 'o',
             DATA_RAMP[:, 0], DATA_RAMP[:, 2], 'o',
             Times, Powers, '-')
