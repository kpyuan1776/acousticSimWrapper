rng('shuffle');

filename = 'siminput_3.h5'; %filename shouldn't be too long: 'simulationinput2.h5' is already too long
outputfilename = 'simoutput_3.h5';

vesselimg = 'image003.jpg';

acousticsim_filename = 'oadata_exampl_3.mat';
reconFilename = 'recon_oaexampl_3.mat';
recon_nnFilename = 'recon_nn_oaexampl_3.mat';
specAngleFilename = 'spectral_angle_3.mat';

dev = DeviceInfo('sim_handheldAcuity');


%%

vessel_img = 1-loadImage(vesselimg);
vessel_img(:,86:88) = 0; %make it quadratic
figure;imagesc(vessel_img);colormap gray;

%% the 512x512x32 grid takes 180s to simulate on GPU and 512x512x20 around 100s

asim = AcousticSim(dev);

asim = asim.setSystem_Handheld();
asim.Nz = 20; %the simulation needs a perfectly matched layer also in z-direction to perform a 3D simulation
asim = makeSimGrid(asim,asim.Nx,asim.Ny,asim.Nz);

asim = asim.setSensorMask_forGPU_mod();
asim.path_to_file = '';
asim.pathToGPUsim = '..\';
%assign time array to speed up simulation
tempdt = 1.3072e-08;
temptidx = 1800;
asim.kgrid.t_array = (0:temptidx-1)*tempdt;

asim = asim.prepareSimulationGPU(vessel_img,1,filename);

[asim, sensor_data, time_kwave, info] = asim.runSimulation_GPU(filename,outputfilename);

%% reconstruction
% needs: calculate_matrix_highpass_1.m
%  Calculate_MatrixMB_Luis.m
%  lsqr_b.m
%  nnls_conjgrad_armijo.m

num_sensor_points = size(sensor_data,1);
time_kwave_pad = time_kwave;

% forward model parameters
fs = 40000000; % sampling frequency
len = 2030; % number of samples
ts = 0:1/fs:(len-1)/fs; % sampling instants
time_res = 2;
c0 = 1530; % speed of sound (ORG: 1520)
image_width = 20e-3;
n = 150;
limits(1) = dev.r_sensor-(image_width)*sqrt(2)/2; % limits for the signal
limits(2) = dev.r_sensor+(image_width)*sqrt(2)/2; % limits for the signal
dx = image_width/n; % increment in x
dt = dx/(time_res*c0); % increment in t employed to make the model-based reconstruction
fac = fs/c0;
pos_start = max(1,int32((limits(1))*fac));
pos_end = min(len,int32((limits(2))*fac));
t = ts(pos_start):dt:ts(pos_end);           % downsampled (& cut) time vector  (less than ts)
sizeT = length(t);
n_angles = 2*n;
angle_sensor_amat = dev.angle_sensor(end:-1:1);
n_iter = 50;
lambda_reg = 1e6;

A_mat = Calculate_MatrixMB_Luis(c0,n,image_width,t,dt,dev.r_sensor,angle_sensor_amat,n_angles);

% now reshape and reconstruct the signals:
sigMat_l = sensor_data';
t_sig = time_kwave_pad;
for proj=1:num_sensor_points
    
    x_sig = squeeze(sigMat_l(:,proj)); % get one signal sequence from one projection
    sigMat_l_cut(:,proj) = interp1(t_sig, x_sig, t);
    %                             sigMat_clean_noise(:,proj,l) = awgn( x_sig_full, SNR, 20*log10( norm(x_sig_full(x_sig_full~=0) ) ) );
    
end
sigMat_rec = sigMat_l_cut;

p_vec = reshape(sigMat_rec, sizeT*num_sensor_points,1);           % reshape raw data matrix to column-major format


nn = 150*150;%n*n;
FILT = lambda_reg*calculate_matrix_highpass_1(150);%n);
L = sparse(FILT);
cell_mat{1,1} = A_mat;
cell_mat{2,1} = L;
A_mat_2 = cell2mat(cell_mat);
p_vec = [p_vec; zeros(nn,1)];
h_vec = lsqr_b(A_mat_2,double(p_vec),n_iter);
h_vec_nn = nnls_conjgrad_armijo(A_mat_2,double(p_vec),zeros(nn,1),0.001,5,3);

Recon_tmp = reshape(h_vec(:,end), 150,150);
Recon_tmp_nn = reshape(h_vec_nn(:,end), 150,150);

%% plot results

figure;
subplot(1,3,1);
imagesc(vessel_img);
title('original image');
subplot(1,3,2);
imagesc(Recon_tmp);
title('Model-based Tik. regularized')
subplot(1,3,3);
imagesc(Recon_tmp_nn);
title('Model-based non-neg. constraint');