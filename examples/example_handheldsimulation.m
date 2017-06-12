%% simple Matlab based simulation of acoustic signals 
% example for simulating acoustic signals recorded by a curved array of 256 transducers
% which originate from an intial pressure distribution p0 (here loaded from
% a mat file). With the fluence simualtion package (not provided by the
% github repository here) new p0 can be created.
% 
% 2017, KPY

addpath('../');

dev = DeviceInfo('sim_handheldAcuity');
dev.angle_sensor = (( -17.5 : -145/255 : -162.5 ))* pi/180;

%%
load('simulated_lightfluence.mat')
p0 = imresize(p0,0.5);

% create acoustic simulation object 
asim = AcousticSim(dev);
asim = asim.setSystem_Handheld();
asim = asim.setSensorMask(1);

% specify the time array to limit number of time points for simulation
tempdt = 1.3072e-08;
temptidx = 1800;
asim.kgrid.t_array = (0:temptidx-1)*tempdt;
[asim, res, time_kwave] = asim.runSimulation(p0);

% res is the signal matrix

