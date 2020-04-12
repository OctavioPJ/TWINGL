from CAREM.Poinki import DirectEvInterface
from CAREM.Citvap import CitvapSourceModel
from PyCAR.PyCIT.FT import MeshFlux
from PyCAR.PyCIT.lector_mest_posta import AddableDict, Geometry, LoadCDP
import re
import numpy as np
import os
import time
import subprocess
from abc import ABC


class TWIGL_MAP(ABC):
    def MapMaterial(self, meshmaterial):
        if meshmaterial in [1, 2, 3]:
            attr = '_reg' + str(meshmaterial)
            return getattr(self, attr)[0][0.0]
        else:
            attr = '_NULL'
        return getattr(self, attr)


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


class DirectEvolutionModel(DirectEvInterface):

    def run_tr(self, Convergence_Criterion):
        try:
            self._CV.run(executable='cittrs.exe', epi_1=Convergence_Criterion)
        except AssertionError:
            time.sleep(2.0)
            self._CV.run(executable='cittrs.exe', epi_1=Convergence_Criterion)
        subprocess.run(['caremdb', '-opt:export', '-val:meshflux', self.__database])
        return

    def evolve(self, dt, *args, **kwargs):
        if self.Equilibrium:
            self.equilibrium_precrs()
            self._Ct = self._C0.copy()
            self._Ct_1 = self._C0.copy()
        else:
            self._Ct_1 = self._Ct.copy()
            self.calculate_precrs(dt)
            self.t += dt

        if self._UseDerivative:
            self.dt = dt

        self.calculate_source(dt)
        self.source_to_file()
        self._send_to_lu_17__()
        self._xsu_mod__()
        self.run_tr(self.FluxConvergence)
        self._get_power__()

        if self.Equilibrium:
            self.Equilibrium = False

        self.Flux_t_1 = self.Flux_t.copy()
        self.Flux_t = MeshFlux(self.__database+'.meshflux', self.Nx, self.Ny, self.Nz, self._NOG)
        return

    def _init_flux(self, groups):
        self._NOG = groups
        self.Flux_t = self.Flux_t_1 = \
            MeshFlux(self.__database.replace('S.cdb', '_eq.cdb.meshflux')
                                               , *self.Geom.Cantidad_de_Nodos(), Ng=self._NOG)
        self.Nx, self.Ny, self.Nz = self.Flux_t.shape[1:4]
        return

    def _init_citvap_model(self, file, seccion4, seccion5, materials, BA, geometry_type, seccion26):
        self._CV = CitvapSourceModel(
            file=file,
            sec4=seccion4,
            sec5=seccion5,
            mat=materials,
            sec26=seccion26,
            GeometricType=geometry_type,
            black_absorber=BA)

        self._CV.PrintMeshMap()
        self._CV.Calculate()
        return

    def _init_mesh_parameters(self):
        self.Geom = Geometry(self._file.replace('S.cii', '0.ci@'))  # HARDCODEADO el *S.cii
        self.Vmesh = self.Geom.Vmesh()
        return

    def _init_parameters(self, _NuFissionMap, _KineticMap):
        self.NuFis = _NuFissionMap()
        self.KM = _KineticMap()
        self.BetaM = np.empty((1, self.Nx, self.Ny, self.Nz, self.nprc))
        self.LambdaM = np.empty((1, self.Nx, self.Ny, self.Nz, self.nprc))
        self.NuFisM = np.empty((1, self.Nx, self.Ny, self.Nz, self._NOG))
        self.VelocityM = np.empty((1, self.Nx, self.Ny, self.Nz, self._NOG))
        self.NuFissionRates = np.empty((1, self.Nx, self.Ny, self.Nz, self._NOG))

        for _x in range(self.Nx):
            for _y in range(self.Ny):
                for _z in range(self.Nz):
                    meshmaterial = self.Geom.sc5[_x][_y][_z]
                    KinParam = self.KM.MapMaterial(meshmaterial)
                    self.BetaM[self.state][_x][_y][_z][:] = KinParam['Beta']
                    self.LambdaM[self.state][_x][_y][_z][:] = KinParam['Decays']
                    self.VelocityM[self.state][_x][_y][_z][:] = KinParam['Neutron Velocity']
                    self.NuFisM[self.state][_x][_y][_z][:] = self.NuFis.MapMaterial(meshmaterial)['XS'][:, 3]

        self.react = self.__R0 = -9251.1E-5
        self.keff = 1 / (1 - self.__R0)
        return 0

    def __init__(self, file, NuFissionMap, KineticMap,
                 WorkingDirectory,
                 seccion_26, materiales, seccion_5, seccion_4,
                 state=0, nprc=1, groups=1,
                 Equilibrium=True,
                 UseDerivative=False,
                 GeometricType='TZ',
                 black_absorber_ID=10,
                 FluxConvergence=None, *args, **kwargs):
        """

        :param file:
        :param NuFissionMap:
        :param KineticMap:
        :param WorkingDirectory:
        :param seccion_26:
        :param materiales:
        :param seccion_5:
        :param seccion_4:
        :param state:
        :param nprc:
        :param groups:
        :param Equilibrium:
        :param UseDerivative:
        :param GeometricType:
        :param black_absorber_ID:
        """
        if FluxConvergence is not None:
            self.FluxConvergence = FluxConvergence
        else:
            self.FluxConvergence = 1.0E-6
        self._geo_type = GeometricType
        self._sec26 = seccion_26
        self._mat = materiales
        self._sec5 = seccion_5
        self._sec4 = seccion_4
        self._BA = black_absorber_ID
        self.t = 0.0
        self.dt = 0.0
        assert isinstance(state, int)
        self.state = state
        self.chi_g = [1.0, 0.0]
        self._file = file
        self.Equilibrium = Equilibrium
        self.nprc = nprc
        self.WorkingDirectory = WorkingDirectory

        self._rootfile = self._file.replace('.cii', '')
        self.__database = self._file.replace('.cii', '.cdb')
        self.__p_parser = re.compile(r'POWER\(WATTS\) +([0-9]\.[0-9]{5}E\+[0-9]{2})')
        self.Pt = 1
        self._Q = {}

        self._C0 = self._Ct = self._Ct_1 = {self.state: {}}

        self._UseDerivative = UseDerivative
        self.__datfile = r'source.dat'

        self._init_mesh_parameters()
        self._init_flux(groups)
        self._init_citvap_model(self._file, self._sec4, self._sec5, self._mat, self._BA, self._geo_type, self._sec26)
        self._init_parameters(NuFissionMap, KineticMap)
        return

    @property
    def Ct_1(self):
        return AddableDict(self._Ct_1).as_array()

    @property
    def C0(self):
        return AddableDict(self._C0).as_array()

    @property
    def Ct(self):
        return AddableDict(self._Ct).as_array()

    @property
    def Q(self):
        return AddableDict(self._Q).as_array()

    def equilibrium_precrs(self, *args, **kwargs):
        if not self._C0[self.state]:
            for nx in range(self.Flux_t.shape[1]):
                self._C0[self.state][nx] = {}
                for ny in range(self.Flux_t.shape[2]):
                    self._C0[self.state][nx][ny] = {}
                    for nz in range(self.Flux_t.shape[3]):
                        FluxL = [self.Flux_t[self.state, nx, ny, nz, group]
                                 for group in range(self.Flux_t.shape[-1])]
                        NuFisL = [self.NuFisM[self.state][nx][ny][nz][group] / self.keff
                                  for group in range(self.Flux_t.shape[-1])]
                        Nu_FluxM = [NuFisL[group] * FluxL[group] * self.Geom.Vmesh()
                                    for group in range(self.Flux_t.shape[-1])]
                        Bet_k = self.BetaM[self.state][nx][ny][nz]
                        Lamb_k = self.LambdaM[self.state][nx][ny][nz]

                        self._C0[self.state][nx][ny][nz] = [Bet_k[prec] / Lamb_k[prec] * sum(Nu_FluxM)
                                                            if Lamb_k[prec] != 0 else 0.0 for prec in range(self.nprc)]
        return

    def calculate_precrs(self, dt):
        for nx in range(self.Nx):
            # self._Ct[self.state][nx] = {}
            for ny in range(self.Ny):
                # self._Ct[self.state][nx][ny] = {}
                for nz in range(self.Nz):
                    FluxL = [self.Flux_t[self.state, nx, ny, nz, group] for group in range(self._NOG)]
                    NuFisL = [self.NuFisM[self.state][nx][ny][nz][group] for group in range(self._NOG)]

                    Bet_k = self.BetaM[self.state][nx][ny][nz]
                    Lamb_k = self.LambdaM[self.state][nx][ny][nz]
                    _c_t_1 = self._Ct_1[self.state][nx][ny][nz]

                    Nu_Flux = [NuFisL[group] * FluxL[group] * self.Geom.Vmesh()
                               for group in range(self._NOG)]

                    self._Ct[self.state][nx][ny][nz] = \
                        [_c_t_1[prec] * np.exp(-Lamb_k[prec] * dt) +
                         Bet_k[prec] / Lamb_k[prec] * sum(Nu_Flux) * (1 - np.exp(-Lamb_k[prec] * dt))
                         if Lamb_k[prec] != 0 else 0.0 for prec in range(self.nprc)]
        return

    def calculate_source(self, dt):
        for group in range(self.Flux_t.shape[-1]):
            self._Q[group] = {self.state: {}}
            for nx in range(self.Flux_t.shape[1]):
                self._Q[group][self.state][nx] = {}
                for ny in range(self.Flux_t.shape[2]):
                    self._Q[group][self.state][nx][ny] = {}
                    for nz in range(self.Flux_t.shape[3]):
                        _Lmk = self.LambdaM[self.state][nx][ny][nz]
                        _C = self._Ct[self.state][nx][ny][nz]

                        inv_V = self.VelocityM[self.state, nx, ny, nz, group]

                        self._Q[group][self.state][nx][ny][nz] = \
                            self.chi_g[group] * sum([_Lmk[prc] * _C[prc] for prc in range(self.nprc)])
                        if self._UseDerivative:
                            self._Q[group][self.state][nx][ny][nz] +=\
                                inv_V / dt * self.Flux_t_1[self.state, nx, ny, nz, group] * self.Vmesh
        return

    def source_to_file(self):
        with open(self.__datfile, 'w') as fod:
            for group in range(self._NOG):
                for nz in range(self.Nz):
                    for ny in range(self.Ny):
                        for nx in range(self.Nx):
                            fod.write('{:15.7E}'.format(self._Q[group][self.state][nx][ny][nz]))
                fod.write('\n')
        return

    def _send_to_lu_17__(self):
        subprocess.run(['copy', '/Y', self.__datfile, 'Source_To_LU17'], shell=True)
        os.chdir(self.WorkingDirectory+'Source_To_LU17')
        msg = subprocess.run(['main.exe', *map(str, [self._NOG, self.Nx, self.Ny, self.Nz])]
                             , shell=True, capture_output=True)
        msg2 = subprocess.run(['copy', '/Y', 'fort.17', '..\\' + self._rootfile + '.src'],
                              shell=True, capture_output=True)
        os.chdir(self.WorkingDirectory)
        return msg

    def _xsu_mod__(self):
        XSU_FILE = '{}'.format(self._CV.RootFile + '.XSU')
        try:
            self._CV.run(executable='pre_cit.exe')
        except AssertionError:
            time.sleep(2.0)
            self._CV.run(executable='pre_cit.exe')

        subprocess.run(['move', XSU_FILE, 'XSU_MOD\\'], shell=True)
        os.chdir(self.WorkingDirectory + 'XSU_MOD')
        if self._UseDerivative:
            time_step = self.dt
        else:
            time_step = 0.0
        msg = subprocess.run(['neutyp', '{}'.format(XSU_FILE), str(self.Equilibrium),
                              *map(lambda a: '{:12.5E}'.format(round(a, ndigits=5)), [time_step, self.t])],
                             shell=True, capture_output=True)
        assert msg.stderr in ['', b''], 'Error en la sobreescritura del archivo ROOT.XSU\n' \
                                        '{}'.format(msg.stderr.decode('utf-8'))
        os.chdir(self.WorkingDirectory)
        msg2 = subprocess.run(['move', '/Y', self.WorkingDirectory+'\\XSU_MOD\\{}'.format(XSU_FILE), XSU_FILE],
                              shell=True, capture_output=True)
        assert msg2.stderr in ['', b''], 'Error en el movimiento del archivo ROOT.XSU\n' \
                                         '{}'.format(msg2.stderr.decode('utf-8'))
        return msg

    def _get_power__(self, file=None):
        if not file:
            fid = open(self._rootfile + '.cio', 'r')
        else:
            fid = open(file, 'r')
        self.Pt = float(self.__p_parser.findall(fid.read())[0])
        fid.close()
        return
