addpath('../');
addpath('../fluenceSimulation');

filename = 'simulationinput.h5';
outputfilename = 'simulationoutput.h5';

dev = DeviceInfo('sim_handheldAcuity');
dev.angle_sensor = (( -17.5 : -145/255 : -162.5 ))* pi/180;


%%
load('simulated_lightfluence.mat')

asim = AcousticSim(dev);
asim = asim.setSystem_Handheld();
asim = asim.setSensorMask(1);
tempdt = 1.3072e-08;
temptidx = 1800;
asim.kgrid.t_array = (0:temptidx-1)*tempdt;
[asim, res, time_kwave] = asim.runSimulation(p0);

%%
