classdef AcousticSim
    %ACOUSTICSIM contains the basic parameters for doing an acoustic
    %simulation given a pressure distribution
    %   How to run a simulation?:
    %       dev = DeviceInfo('sim_mousecart');or 'sim_handheldAcuity'
    %       asim = AcousticSim(dev);
    %       asim = asim.makeSimGrid(150,256,256,32);
    %       asim = asim.setSensorMask(1);
    %       [asim, res, time_kwave] = asim.runSimulation(p0);
    
    properties
        devInfo; %of class DeviceInfo
        x = 30.4e-3;            % total grid size [m]
        y = 30.4e-3;           % total grid size [m]
        z = 4*5e-3;
        dx = 0.00013333;
        dy = 0.00013333;
        dz = 2*0.00013333;
        Nx = 256;
        Ny = 256;
        Nz = 96;
        sensor_pos = [0, 0];        % [m]
        slice_thickness = 2;
        kgrid;
        medium;
        source;
        sensor;
        smooth_trg = 1;
        path_to_file='';
        pathToGPUsim='';
    end
    
    methods
        function ASobj = AcousticSim(deviceInfo)
            ASobj.devInfo = deviceInfo;
            ASobj.medium.sound_speed = 1530;
            ASobj.medium.density = 1000;
        end
        
        function obj = setSystem_Handheld(obj)
           obj.Nx = 512;
           obj.Ny = 512;
           obj.Nz = 32;
           obj.dx = 0.00013333/2;
           obj.dy = 0.00013333/2;
           obj.dz = 4*0.00013333;
           obj.kgrid = makeGrid(obj.Nx, obj.dx, obj.Ny, obj.dy, obj.Nz, obj.dz);
        end
        
        %Nx=Ny=256, Nz=96
        function obj = makeSimGrid(obj,Nx,Ny,Nz)

            if mod(Nx,2)
                Nx = Nx+1;
            end
            %pad_x = (Nx - n)/2;
    
            if mod(Ny,2)
                Ny = Ny+1;
            end
            %pad_y = (Ny - n)/2;

            if mod(Nz,2)
                Nz = Nz+1;
            end
            %pad_z = (Nz-obj.slice_thickness)/2;
            
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.Nz = Nz;
            obj.kgrid = makeGrid(obj.Nx, obj.dx, obj.Ny, obj.dy, obj.Nz, obj.dz);
        end
        
        %         function obj = makeSimDomain()
        %            obj.medium.sound_speed = 1530;
        %            obj.medium.density = 1000;
        %         end
        
        function obj = setSensorMask(obj, varargin)
            
            if nargin>1
                nzDets = varargin{1};
            else
                nzDets = 24;%Nz/2;
            end
            
            sensor_radius = obj.devInfo.r_sensor;
            num_sensor_points = obj.devInfo.proj;
            sensor_pos = obj.sensor_pos;
            ang_ini =  obj.devInfo.angle_sensor(1);
            ang_end = obj.devInfo.angle_sensor(end);
            sensor_angle = abs(ang_end - ang_ini);
            
            
            cartsensors_all = [];
            cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points, sensor_pos, sensor_angle);
            cartsensors_temp = [cart_sensor_mask; zeros(1, numel(cart_sensor_mask)./2)];
            cartsensors_all = [cartsensors_all, cartsensors_temp];
            if nzDets>1
                for nz_dist=1:(nzDets-1)
                    cart_sensor_mask = makeCartCircle(sqrt(sensor_radius^2-(nz_dist*obj.dz)^2), num_sensor_points, sensor_pos, sensor_angle);
                    cartsensors_temp = [cart_sensor_mask; zeros(1, numel(cart_sensor_mask)./2)+nz_dist*obj.dz];
                    cartsensors_all = [cartsensors_all, cartsensors_temp];
                    cartsensors_temp = [cart_sensor_mask; zeros(1, numel(cart_sensor_mask)./2)-nz_dist*obj.dz];
                    cartsensors_all = [cartsensors_all, cartsensors_temp];
                end
                
                cart_sensor_mask = makeCartCircle(sqrt(sensor_radius^2-(nz_dist*obj.dz)^2), num_sensor_points, sensor_pos, sensor_angle);
                cartsensors_temp = [cart_sensor_mask; zeros(1, numel(cart_sensor_mask)./2)-nz_dist*obj.dz];
                cartsensors_all = [cartsensors_all, cartsensors_temp];
            end
            obj.sensor.mask = cartsensors_all;
        end
        
        function obj = setSensorMask_forGPU(obj, varargin)
            sensor_radius = obj.devInfo.r_sensor;%11e-3;     % [m]
            ang_ini = obj.devInfo.angle_sensor(1);%-0.7853981633;
            ang_end = obj.devInfo.angle_sensor(end);%3.926990817079;
            % ang_step = dev.angle_sensor(2)-dev.angle_sensor(1);%0.0184799567858;
            % angle_sensor = dev.angle_sensor;%[ang_ini:ang_step:ang_end];
            sensor_angle = abs(ang_end - ang_ini);
            
            sensor_pos = [0, 0];        % [m]
            num_sensor_points = 256;
            
            % make central detector slice
            cartsensors_all = [];
            cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points, sensor_pos, sensor_angle);
            cartsensors_temp = [cart_sensor_mask; zeros(1, numel(cart_sensor_mask)./2)];
            cartsensors_all = [cartsensors_all, cartsensors_temp];
            
            obj.sensor.mask = cartsensors_all;
            
            [data_mask, order, reorder] = cart2grid(obj.kgrid, obj.sensor.mask);
            %voxelPlot(data_mask,'Color',[1 0 0],'Transparency', 0.5);
            obj.sensor.mask = data_mask;
        end
        
        function obj = setSensorMask_forGPU_mod(obj, varargin)
            sensor_radius = obj.devInfo.r_sensor;%11e-3;     % [m]
            ang_ini = obj.devInfo.angle_sensor(1);%-0.7853981633;
            ang_end = obj.devInfo.angle_sensor(end);%3.926990817079;
            sensor_angle = abs(ang_end - ang_ini);
            
            sensor_pos = [0, 0];        % [m]
            num_sensor_points = 256;
            
            sensormask = zeros(obj.kgrid.Nx,obj.kgrid.Ny,obj.kgrid.Nz);
            middle_x = round(obj.kgrid.Nx/2);
            middle_y = round(obj.kgrid.Ny/2);
            middle_z = round(obj.kgrid.Nz/2)+1;
            
            % make central detector slice

            for i=1:num_sensor_points
                x_pos = middle_x+round(sensor_radius*sin(obj.devInfo.angle_sensor(i))/obj.kgrid.dx);
                y_pos = middle_y+round(sensor_radius*cos(obj.devInfo.angle_sensor(i))/obj.kgrid.dy);
                z_pos = middle_z;
                sensormask(x_pos,y_pos,z_pos) = 1;
            end

            obj.sensor.mask = sensormask;
        end
        
        %function obj = rotateSensorMask(obj, varargin)
        
        function obj = prepareSimulationGPU(obj,p0,smooth_trg,filename)
            
             nx = size(p0, 1);
             ny = size(p0, 2);
            
            pad_x = floor((obj.Nx - nx)/2);
            pad_y = floor((obj.Ny - ny)/2);
            pad_z = floor((obj.Nz-obj.slice_thickness)/2);
            
            p0_init = zeros(obj.Nx,obj.Ny,obj.Nz);
            p0_init = repmat(p0, [1, 1, obj.slice_thickness]);
            p0_init = padarray(p0_init, [pad_x, pad_y, pad_z], 'both');
            p0_init = p0_init;
            obj.source.p0 = p0_init;
            
            % input arguments
            input_args = {'PlotLayout', false, 'PMLInside', true, 'PlotPML', false, ...
                'DataCast', 'single', 'CartInterp', 'nearest', 'Smooth', logical(smooth_trg)};
            
            % run the simulation
            kspaceFirstOrder3D(obj.kgrid, obj.medium, obj.source, obj.sensor, 'SaveToDisk', filename);
        end
        
        function [obj, res, time_kwave, varargout] = runSimulation_GPU(obj, inputfilename,outputfilename,varargin)                      
            
            if nargin>4
                pad_num = varargin{1};
            else
                pad_num = 0;%round(((r_sensor - sensor_radius)./1530)./dt_kwave);
            end
            
            cmd = [[obj.pathToGPUsim 'kwavegpusims\kspaceFirstOrder3D-CUDA.exe'], ...
                ' -i ', [obj.path_to_file inputfilename], ...
                ' -o ', [obj.path_to_file outputfilename]];
            exec = @(cmd) system(cmd);
            [status] = exec(cmd);
            assert(status == 0);
            
            sensor_data2.p = h5read(outputfilename,'/p');
            
            tempdt = h5read(outputfilename,'/dt');
            temptidx = h5read(outputfilename,'/t_index');
            time_kwave = 0:tempdt:(tempdt*double(temptidx-1));
            
            dt_kwave = time_kwave(2) - time_kwave(1);
            
            
            time_kwave_pad = [time_kwave, (time_kwave(end)+dt_kwave):dt_kwave:(time_kwave(end)+pad_num*dt_kwave)];
            
            if strcmp(obj.devInfo.devname,'sim_mousecart')
                numreorder = reorderSensorNumbering();
            elseif strcmp(obj.devInfo.devname,'sim_handheldAcuity')
                numreorder = reorderSensorNumbering_hh();
            end
%             for i=1:256
%                 res(numreorder(i),:) = sensor_data2.p((i),:);
%             end
            res = sensor_data2.p;
            
            if nargout>3
                info.dt_kwave = dt_kwave;
                info.time_kwave_pad = time_kwave_pad;
               varargout{1} = info;
            end
        end
        
        
        function plotSensorArray(obj)
            [data_mask, order, reorder] = cart2grid(obj.kgrid, obj.sensor.mask);
            voxelPlot(data_mask,'Color',[1 0 0],'Transparency', 0.5);
        end
        
        function [obj, res, time_kwave] = runSimulation(obj, img_init,varargin)
            
            if nargin>2
                flag = varargin{1};
            else
                flag = 'matlabsim'
            end
            
            if strcmp(flag,'matlabsim')
                for wav = 1:size(img_init,3)
                    p0 = img_init(:,:,wav);
                    [obj, sensor_data, time_kwave] = obj.runSimulationSingleSlice(p0);
                    res(:,:,wav) = sensor_data;
                end
            end
            
        end
        
        function [obj, sensor_data, time_kwave] = runSimulationSingleSlice(obj, p0)
            
            n = size(p0, 1);
            
            pad_x = (obj.Nx - n)/2;
            pad_y = (obj.Ny - n)/2;
            pad_z = (obj.Nz-obj.slice_thickness)/2;
            
            p0_init = repmat(p0, [1, 1, obj.slice_thickness]);
            p0_init = padarray(p0_init, [pad_x, pad_y, pad_z], 'both');
           % p0_init = p0_init;
            obj.source.p0 = p0_init;
            
            % input arguments
            input_args = {'PlotLayout', false, 'PMLInside', true, 'PlotPML', false,'PlotSim',false, ...
                'DataCast', 'single', 'CartInterp', 'nearest', 'Smooth', logical(obj.smooth_trg)};
            
            % run the simulation
            sensor_data = kspaceFirstOrder3D(obj.kgrid, obj.medium, obj.source, obj.sensor, input_args{:});
            [time_kwave, ~] = makeTime(obj.kgrid, obj.medium.sound_speed);
        end
    end
    
end

