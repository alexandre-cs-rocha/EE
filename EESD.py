import numpy as np
import random
from pathlib import Path
from dss import DSS as dss_engine

def InitializeDSS() -> tuple:
    DSSObj = dss_engine
    flag = DSSObj.Start(0)
    if flag:
        print('OpenDSS COM Interface initialized succesfully.')
        
    else:
        print('OpenDSS COMInterface failed to start.')
        
    #Set up the interface variables - Comunication OpenDSS with Python
#    DSSActiveCktElement = DSSObj.ActiveCircuit.ActiveCktElement
    DSSText = DSSObj.Text
    DSSCircuit = DSSObj.ActiveCircuit

#    DSSBus          =   DSSCircuit.ActiveBus
#    DSSTransformers =   DSSCircuit.Transformers
#    DSSSolution     =   DSSCircuit.Solution
#    DSSMeters       =   DSSCircuit.Meters
#    DSSLines        =   DSSCircuit.Lines
#    DSSLoads        =   DSSCircuit.Loads
#    DSSVsources     =   DSSCircuit.Vsources
#    DSSPVSystems    =   DSSCircuit.PVSystems
    DSSMonitors     =   DSSCircuit.Monitors
            
    return DSSCircuit, DSSText, DSSObj, DSSMonitors

def desvio_padrao(vet_med: np.array) -> np.array:
    pr = 0.002
    dp = (pr * vet_med)/(3*1000)
    return dp

path = Path(__file__)
CurrentFolder = path.parent
MasterFile = CurrentFolder / '13Bus' / 'IEEE13Nodeckt.dss'

num_med = 5

DSSCircuit, DSSText, DSSObj, DSSMonitors = InitializeDSS()
DSSObj = dss_engine
DSSCircuit = DSSObj.ActiveCircuit
DSSMonitors = DSSCircuit.Monitors
DSSText.Command = f'Compile {MasterFile}'


elementos = random.sample(DSSCircuit.Lines.AllNames, num_med)

for i, elem in enumerate(elementos):
    DSSText.Command = f'New Monitor.{i} element=Line.{elem}, terminal=1, mode=32'


DSSText.Command = 'Solve'
#print(DSSCircuit.Generators.AllNames)
DSSMonitors.SampleAll()
DSSMonitors.SaveAll()

DSSMonitors.First
vet_med = np.zeros(num_med)
for i in range(DSSCircuit.Monitors.Count):
    DSSCircuit.SetActiveElement(DSSMonitors.Element)
    DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
    vet_med[i] = DSSMonitors.Channel(1) / (DSSCircuit.ActiveBus.kVBase*1000)

    DSSMonitors.Next

x = np.array([[1, 0] for _ in range(DSSCircuit.NumBuses-1)])

des_pad = desvio_padrao((vet_med * 0.002)/300)
