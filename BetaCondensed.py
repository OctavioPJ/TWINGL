import numpy as np
from PyCAR.PyCIT.FT import MeshFlux
from TWINGL.TwinglModel import Nu_Fission_Map, Kinetic_Map
from PyCAR.PyCIT.lector_mest_posta import Geometry


def BetaCalculator(FluxFile, Nx, Ny, Nz, Ng, GeometryFile):
    Flux = MeshFlux(FluxFile, Nx, Ny, Nz, Ng)
    Geo = Geometry(GeometryFile)
    chi_p = [1.0, 0.0]
    NfMap = Nu_Fission_Map()
    KnMap = Kinetic_Map()
    FissionRates = np.empty((Nx, Ny, Nz, Ng))
    VelocityRates = np.empty((Nx, Ny, Nz, Ng))

    OneGroupFlux = np.empty((Nx, Ny, Nz))
    OneGroupFiss = np.empty((Nx, Ny, Nz))
    DelayedFiss = np.empty((Nx, Ny, Nz))
    OneGroupVel = np.empty((Nx, Ny, Nz))

    for x in range(Nx):
        for y in range(Ny):
            for z in range(Nz):
                for group in range(Ng):
                    FissionRates[x, y, z, group] = \
                        NfMap.MapMaterial(Geo.sc5[x, y, z])['XS'][group, 3] * Flux[0, x, y, z, group] \
                        * Geo.Vmesh()

                    VelocityRates[x, y, z, group] = \
                        KnMap.MapMaterial(Geo.sc5[x, y, z])['Neutron Velocity'][group] * Flux[0, x, y, z, group] \
                        * Geo.Vmesh()

                OneGroupFlux[x, y, z] = Flux[0, x, y, z, :].sum()
                OneGroupFiss[x, y, z] = FissionRates[x, y, z, :].sum()
                DelayedFiss[x, y, z] = 0.0075 * FissionRates[x, y, z, :].sum() * Flux[0, x, y, z, :].sum()
                # segun otto, esta formula estaría bien, dios sabrá
                OneGroupVel[x, y, z] = VelocityRates[x, y, z, :].sum()

    Lambda = ((OneGroupVel * OneGroupFlux).sum() /
              (OneGroupFiss * OneGroupFlux).sum())
    Beta = DelayedFiss.sum() /\
           (OneGroupFiss * OneGroupFlux).sum()
    return Beta, Lambda


if __name__ == '__main__':
    FluxFile = 'twigl_eq.cdb.meshflux'
    GeometryFile = 'twigl0.ci@'
    Beta, Lambda = BetaCalculator(FluxFile, 48, 48, 1, 2, GeometryFile)
