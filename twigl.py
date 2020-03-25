# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:33:21 2019

@author: CNEA
"""
import sys
if 'D:\\' not in sys.path: sys.path.append('D:\\')
from PyCAR.PyCIT.lector_mest_posta import LoadCDP,get_reactivity_pro as grp
from PyCAR.PyCIT.FT import MeshFlux
from functools import reduce
from operator import add
import matplotlib.pyplot as plt

import os

BLACKABSORBER = 10
DEFAULT_STATE = 0
DEFAULT_BURNUP = 0


class TWIGL_MAP(object):
    def MapMaterial(self,meshmaterial,**kwargs):
        attr = '_reg'+str(meshmaterial)
        try:
            return getattr(self,attr)[0][0.0]
        except AttributeError as err:
            return getattr(self,'_NULL')


class Nu_Fission_Map(TWIGL_MAP):
    wdir = 'D:\\TWINGL\\'

    def __init__(self,**kwargs):
        if 'wdir' in kwargs:
            setattr(self,'wdir',kwargs['wdir'])
            
        wdir = getattr(self,'wdir')
        self._reg1 = LoadCDP(wdir+'REGION1.cdp')
        self._reg2 = LoadCDP(wdir+'REGION2.cdp')
        self._reg3 = LoadCDP(wdir+'REGION3.cdp')
        
        self._NULL = {'XS': np.zeros(self._reg3[0][0.0]['XS'].shape),\
                      'SM': np.zeros(self._reg3[0][0.0]['SM'].shape)}
                       

class Kinetic_Map(TWIGL_MAP):
    def __init__(self,**kwargs):
        v1 = 10**7
        v2 = 2*(10**5)
        self._reg1 = self._reg2 = self._reg3 = {0:{0.0:{'Beta':np.array([0.0075]),\
                                                        'Decays':np.array([0.08]),\
                                                        'Neutron Velocity':np.array([1/v1,1/v2])}}}
        
        self._NULL = {'Beta':[0],'Decays':[0],'Neutron Velocity':[0,0]}
        

if __name__ == '__main__':
    
    from citvappru.SourceCAREM2_1_1 import Geometry,PrecCalc
    import numpy as np
    import re
    os.chdir('D:\\TWINGL\\')
    FILENAME = 'twigl.cii'
    ROOTFILE = FILENAME[:-4]
    DATABASE = ROOTFILE+'_eq.cdb'
    
    Geom = Geometry(FILENAME.replace('.cii','0.cii'))
    NuFis = Nu_Fission_Map()   
    KM    = Kinetic_Map()   
    _NOG = 2
    Flux = MeshFlux(DATABASE+'.meshflux',*Geom.Cantidad_de_Nodos(),_NOG)
    Nx,Ny,Nz = Flux.shape[1:4]

    BetaM = np.empty((1,Nx,Ny,Nz,1))
    LambdaM = np.empty((1,Nx,Ny,Nz,1))
    NuFisM = np.empty((1,Nx,Ny,Nz,_NOG))
    VelocityM = np.empty((1,Nx,Ny,Nz,_NOG))
    
    state = 0
       
    Vmesh = Geom.Vmesh()
    
    react = -9251.1#reactividad inicial inicial, no inicial
    keff = 1/(1-react/100000)

    for _x in range(Nx):
        for _y in range(Ny):
            for _z in range(Nz):
                
                meshmaterial = Geom.sc5[_x][_y][_z]
                
                KinParam = KM.MapMaterial(meshmaterial)
                
                BetaM[state][_x][_y][_z][:] = KinParam['Beta']
                
                LambdaM[state][_x][_y][_z][:] = KinParam['Decays']
                                 
                VelocityM[state][_x][_y][_z][:] = KinParam['Neutron Velocity']
                #CASO ESPECIAL DONDE  CONVIENE DIVIDIR AQUI POR EL KEFF
                NuFisM[state][_x][_y][_z][:] =  NuFis.MapMaterial(meshmaterial)['XS'][:,3]/keff

                   
    C0 = {}
    NPRC = BetaM.shape[-1]
    C0[state]={}
    for nx in range(Flux.shape[1]):
        C0[state][nx] = {}
        for ny in range(Flux.shape[2]):
            C0[state][nx][ny] = {}
            for nz in range(Flux.shape[3]):
                
                FluxL = [Flux[state,nx,ny,nz,group] for group in range(Flux.shape[-1])]
                NuFisL= [NuFisM[state][nx][ny][nz][group] for group in range(Flux.shape[-1])]
                Nu_FluxM=[NuFisL[group]*FluxL[group]*Vmesh for group in range(Flux.shape[-1])]
    
                Bet_k = BetaM[state][nx][ny][nz]
                Lamb_k = LambdaM[state][nx][ny][nz]
                Nu_Flux = sum(Nu_FluxM)
    
                C0[state][nx][ny][nz] = [Bet_k[prec]*Nu_Flux/Lamb_k[prec]\
                                             if Lamb_k[prec] != 0 else 0.0   \
                                             for prec in range(NPRC) ]
                
    p=re.compile(r'POWER\(WATTS\)\s+([0-9]\.[0-9]{5}E\+[0-9]{2})')
#    powers=[[1e+06]]    
    powers = []

    Q={}
    EQUILIBRIUM=True
    chi_g = [1.0,0.0]

#    HARDCODEADO
    v1 = 10**7
    v2 = 2*(10**5)
    dt=0.01
    tfinal=0.5
    DATABASE=ROOTFILE+'S.cdb'
    Times = np.arange(0,tfinal+2*dt,dt)
    FAILED_COMPLETE_TEST = False
    DERIVATIVE = True
    for t in Times:
        if EQUILIBRIUM:
            C_t_1 = C0.copy()
            C_t = C0.copy()
            Flux_t_1 = Flux
        else:
            TFlux = MeshFlux(DATABASE+'.meshflux',Nx,Ny,Nz,_NOG)
            C_t = PrecCalc(C_t_1,TFlux,NuFisM,Vmesh,dt\
                           ,LambdaM = LambdaM\
                           ,BetaM = BetaM)
            
            C_t_1 = C_t.copy()
            Flux_t_1 = TFlux
            
        np.savetxt('ramp_fast_{}.dat'.format(t),Flux_t_1[0,:,:,0,0])            
        np.savetxt('ramp_term_{}.dat'.format(t),Flux_t_1[0,:,:,0,1])            
        
        for group in range(Flux.shape[-1]):
            Q[group]={}
            for state in range(Flux.shape[0]):
                Q[group][state]={}
                for nx in range(Flux.shape[1]):
                    Q[group][state][nx] = {}
                    for ny in range(Flux.shape[2]):
                        Q[group][state][nx][ny] = {}
                        for nz in range(Flux.shape[3]):
                            _Lmk= LambdaM[state][nx][ny][nz]
                            _C  = C_t[state][nx][ny][nz]
                            
                            _invV = VelocityM[state,nx,ny,nz,group]
                            T_1Flux = Flux_t_1[state,nx,ny,nz,group]
                            
                            Q[group][state][nx][ny][nz] = \
                            chi_g[group]*sum([ _Lmk[prc]*_C[prc] for prc in range(NPRC)])\
                            +     _invV/dt * T_1Flux *Vmesh * DERIVATIVE
    
        with open('source.dat','w') as fod:
            for group in range(_NOG):
                for state in range(Flux.shape[0]):
                    for nz in range(Flux.shape[3]):
                        for ny in range(Flux.shape[2]):
                            for nx in range(Flux.shape[1]):
                                fod.write('{:15.7E}'.format(Q[group][state][nx][ny][nz])) 
                fod.write('\n')
                    
        OsArgs = (FILENAME.replace('.cii','S.cii'),\
                  _NOG,Nx,Ny,Nz,EQUILIBRIUM,\
                  *['{:12.5E}'.format(_d) for _d in [dt*DERIVATIVE,t]])
                  #PARA HACER UNA PERTURBACION ESCALOR,COLOCAR t=1
                  #PARA CALCULAR SIN DERIVADA TEMPORAL dt=0
        try:
            os.system('ExTwigl.bat '+' '.join(map(str,OsArgs)))
            fid = open(ROOTFILE+'S.cio','r')
        #        powers.append([float(pi.groups()[0]) for pi in p.finditer(fid.read())])
            powers.append(float(next(p.finditer(fid.read())).groups()[0]))
            print(powers[-1])
            fid.close()
            if EQUILIBRIUM: EQUILIBRIUM = False     
        except StopIteration:
            FAILED_COMPLETE_TEST = True
            break
            pass
        
    Normalized_Pow = np.array(powers)/1E+06
    if FAILED_COMPLETE_TEST: Times = [k*dt for k in range(len(powers))]
    
    DATA_STEP = np.array([[0,1,1],
                          [0.1,2.023,2.061],
                          [0.2,2.051,2.080],
                          [0.3,2.074,2.097],
                          [0.4,2.096,2.114],
                          [0.5,2.117,2.132]])
    
    DATA_RAMP = np.array([[0,1,1],
                          [0.1,1.311,1.308],
                          [0.2,1.952,1.961],
                          [0.3,2.045,2.076],
                          [0.4,2.069,2.093],
                          [0.5,2.092,2.111]])                
    
    
    plt.plot(Times,Normalized_Pow,DATA_STEP[:,0],DATA_STEP[:,1],'ro',DATA_STEP[:,0],DATA_STEP[:,2],'go')
        
#    plt.plot(Times,Normalized_Pow,DATA_RAMP[:,0],DATA_RAMP[:,1],'ro',DATA_RAMP[:,0],DATA_RAMP[:,2],'go')