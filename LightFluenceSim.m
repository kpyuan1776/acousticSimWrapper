classdef LightFluenceSim
    %LIGHTFLUENCESIM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n;
        wavelengths = 700:10:900;
        spectra;
        datapath;
        voxelsize_in_cm = 0.0133; 
    end
    
    methods
        function lobj = LightFluenceSim()
            lobj.n = 150;
            lobj.datapath = 'spectra';
            lobj.wavelengths = 700:10:900;
        end
        
        function obj = loadSpectra(obj, varargin)
            if nargin>1
                fatspecfile = varargin{1};
                waterspecfile = varargin{2};
                obj.spectra = LoadSpectra(obj.datapath,obj.wavelengths,fatspecfile,waterspecfile); 
            else
                obj.spectra = LoadSpectra(obj.datapath,obj.wavelengths);
            end
        end
        
        function [obj, MSPimg] = runSimulationRndsO2(obj)
            
            wavs = obj.wavelengths(16);
            oxy_stds = 0.4;
            
            iter_leaning = 1;
            leaning_mat = [-0.40:0.05:0.40];
            iter_gran = 1;
            granulation = [15];
            
            iter_mua = 1;
            mua_mat = [0.3];
            iter_mus = 1;
            mus_mat = [10];
            
            deviations = [ oxy_stds(1), 3, 0.1];
            [U_d_container U_d_container_norm MSPimg Gold_StandardSO2, m_a, mu_s] =...
                createRandomOPMapN_check_revisions(obj.datapath, leaning_mat(iter_leaning), granulation(iter_gran), 'fvm', ...
                deviations, mua_mat(iter_mua), mus_mat(iter_mus),wavs);
            MSPimg = imresize(MSPimg,obj.n/150);
        end
        
        function [obj, MSPimg] = runSimulationRndsO2_hh(obj)
            
            wavs = obj.wavelengths(16);
            oxy_stds = 0.4;
            
            iter_leaning = 1;
            leaning_mat = [-0.40:0.05:0.40];
            iter_gran = 1;
            granulation = [15];
            
            iter_mua = 1;
            mua_mat = [0.3];
            iter_mus = 1;
            mus_mat = [10];
            
            deviations = [ oxy_stds(1), 3, 0.1];
            [U_d_container U_d_container_norm MSPimg Gold_StandardSO2, m_a, mu_s] =...
                createSmoothRandomOPMapN_HH(obj.spectra,obj.n, leaning_mat(iter_leaning), granulation(iter_gran), 'fvm', ...
                deviations, mua_mat(iter_mua), mus_mat(iter_mus),wavs);
            %MSPimg = imresize(MSPimg,obj.n/150);
        end
        
        function [simdomain] = defineSimulationDomain(obj)
            % Define a disk tissue geometry
            simdomain = zeros(obj.n,obj.n);
            radius_x = round(obj.n/2)*3.0;
            radius_y = round(obj.n/2)*1.0;
            centerpoint = [round(obj.n/2),round(obj.n/2)*2];
            [X, Y] = meshgrid(1:size(simdomain,1),1:size(simdomain,2));
            ind = find((abs(X-centerpoint(1)).^2+abs(Y-centerpoint(2)).^2  < radius_y*radius_x));

            simdomain(ind) = 1;

        end
        
        function [obj, optprops] = nsimNoisyBloodBckgrd(obj, domain, varargin)
            if nargin>2
                oxy_std = varargin{1};
            else
                oxy_std = 0.4;
                scat_std = 3;
                abs_std = 0.1;
                mus_mat = 10;
                mua_mat = 0.3;
                n_points_r = 15;
                leaning = -0.40;
            end
            
            % In order to avoid
            if abs(n_points_r/2 - round(n_points_r/2))>0
                n_points_r = n_points_r+1;
            end
            I = randn(n_points_r,n_points_r);
            
            % The following code makes the sO2 map symmetric (not used here)
            %     I2 = zeros(size(I));
            %     I2(:,1:round(size(I2,1)/2)) = I(:,1:round(size(I2,1)/2));
            %     I2(:,round(size(I2,1)/2)+1:end) = I(:,round(size(I2,1)/2):-1:1);
            
            rand_image = 0.5 + leaning + oxy_std.*I;
            [Xri Yri] = meshgrid(linspace(1,n_points_r,n_points_r),linspace(1,n_points_r,n_points_r));
            [X150 Y150] = meshgrid(linspace(1,n_points_r,obj.n),linspace(1,n_points_r,obj.n));
            
            image_oxy_level_rand = domain.*interp2(Xri,Yri,rand_image,X150,Y150,'cubic');
            image_oxy_level_rand(image_oxy_level_rand>1) = 1;
            image_oxy_level_rand(image_oxy_level_rand<0) = 0;
            image_deoxy_level_rand = (1-image_oxy_level_rand).*domain;
            
            % scattering map
            % scatter coefficient is random: mean 10, (5-15)
            scatter_mean = mus_mat;
            scatter_max_dev = mus_mat;   %maximum allowedpeak to peak deviation from the mean
            scatter_std = scat_std;  % standart deviation
            scatter_rand = scatter_mean + scatter_std.*randn(n_points_r,n_points_r);
            [Xri Yri] = meshgrid(linspace(1,n_points_r,n_points_r),linspace(1,n_points_r,n_points_r));
            [X150 Y150] = meshgrid(linspace(1,n_points_r,obj.n),linspace(1,n_points_r,obj.n));
            scatter_coef = domain.*interp2(Xri,Yri,scatter_rand,X150,Y150,'cubic');
            % make sure the numbers are reasonable (e.g. not negative coefficients)
            scatter_coef(scatter_coef < scatter_mean - scatter_max_dev./2) = scatter_mean - scatter_max_dev./2;
            scatter_coef(scatter_coef > scatter_mean + scatter_max_dev./2) = scatter_mean + scatter_max_dev./2;
            scatter_coef = domain.*scatter_coef;% sets zeroes around the tissue model
     
            
            % absorption map
            % Define the absorption coefficient mu_a (cm-1) at 800 nm that will be multiplied with the image structure
            m_a_mean = mua_mat;
            m_a_max_dev = mua_mat;  %maximum allowed peak to peak deviation on uniform distribution
            m_a_std = abs_std; % standart deviation
            absorption_rand = m_a_mean + m_a_std.*randn(n_points_r,n_points_r);
            [Xri Yri] = meshgrid(linspace(1,n_points_r,n_points_r),linspace(1,n_points_r,n_points_r));
            [X150 Y150] = meshgrid(linspace(1,n_points_r,obj.n),linspace(1,n_points_r,obj.n));
            m_a = domain.*interp2(Xri,Yri,absorption_rand,X150,Y150,'cubic');
            % make sure the numbers are reasonable
            m_a(m_a < m_a_mean - m_a_max_dev./2) = m_a_mean - m_a_max_dev./2;
            m_a(m_a > m_a_mean + m_a_max_dev./2) = m_a_mean + m_a_max_dev./2;
            m_a = domain.*m_a; % sets zeroes around the tissue model
            
            for i = 1:size(m_a,1)
                for j = 1:size(m_a, 2)
                    deoxy_spectra_spat(i,j,:) = m_a(i,j)*(obj.spectra(1,:));
                    oxy_spectra_spat(i,j,:) = m_a(i,j)*(obj.spectra(2,:));
                end
            end
            
            mu_a_spat = zeros(size(deoxy_spectra_spat));
            for w=1:size(mu_a_spat,3)
                mu_a_spat(:,:,w) = image_deoxy_level_rand(:,:).*deoxy_spectra_spat(:,:,w) + image_oxy_level_rand(:,:).*oxy_spectra_spat(:,:,w);
            end
            
            mu_redscatt = zeros(size(deoxy_spectra_spat));
            for w = 1:size(mu_redscatt,3)
                mu_redscatt(:,:,w) = scatter_coef;
            end
            
            %optprops.deoxy_spectra_spat = deoxy_spectra_spat;
            %optprops.oxy_spectra_spat = oxy_spectra_spat;
            optprops.so2 = image_oxy_level_rand;
            optprops.mua = mu_a_spat;
            optprops.mus = mu_redscatt;
        end
        
        function [obj optprops_out] = addSomeVessels(obj, domain, optprops, strength, varargin)
            
            if nargin>4
                filename_vesselimg = varargin{1};
                sO2 = varargin{2}; %0<sO2<1
            else
                filename_vesselimg = 'vessel_2.png';%'vessel_1.png';
                sO2 = 0.9;
            end
            
            vessel_img = 1-loadImage(filename_vesselimg);
            %vessel_img = imgaussfilt(vessel_img,1)
            %figure;imagesc(B);
            
            [obj, mu_s] = obj.returnRedScatteringCoeff(obj.wavelengths,'blood');
            m_a = strength*0.3.*vessel_img;
            for i = 1:size(m_a,1)
                for j = 1:size(m_a, 2)
                    if vessel_img(i,j)>0
                        %fprintf('test %i %i \n',i,j);
                        %optprops.deoxy_spectra_spat(i,j,:) = (1-sO2).*m_a(i,j)*(obj.spectra(1,:));
                        %optprops.oxy_spectra_spat(i,j,:) = sO2.*m_a(i,j)*(obj.spectra(2,:));
                        optprops.mua(i,j,:) = (1-sO2).*m_a(i,j)*(obj.spectra(1,:))+sO2.*m_a(i,j)*(obj.spectra(2,:));
                        optprops.mus(i,j,:) = mu_s(:);
                        optprops.so2(i,j,:) = sO2;
                    end
                end
            end
            
            optprops_out = optprops;
        end
        
        function [obj optprops_out] = addFatnWater(obj, domain, optprops, strength_fat,strength_water, varargin)
            if nargin>5
                filename_fatimg = varargin{1};
                filename_waterimg = varargin{2};
            else
                filename_fatimg = 'layers_fat2.png';
                filename_waterimg = 'layers_water2.png';
            end
            
            fat_img = 1-loadImage(filename_fatimg);
            water_img = 1-loadImage(filename_waterimg);
            
            [obj, mu_s] = obj.returnRedScatteringCoeff(obj.wavelengths,'fat');
            m_a_fat = strength_fat.*fat_img;
            m_a_water = strength_water.*water_img;
            for i = 1:size(m_a_fat,1)
                for j = 1:size(m_a_fat, 2)
                    if fat_img(i,j)>0
                        %fprintf('test %i %i \n',i,j);
                        optprops.mua(i,j,:) = squeeze(optprops.mua(i,j,:))+squeeze(m_a_fat(i,j)*(obj.spectra(3,:)'));
                        optprops.mus(i,j,:) = mu_s(:);
                    elseif water_img(i,j)>0
                        optprops.mua(i,j,:) = squeeze(optprops.mua(i,j,:))+squeeze(m_a_water(i,j)*(obj.spectra(4,:)'));
                        optprops.mus(i,j,:) = mu_s(:);
                    end
                end
            end
            
            optprops_out = optprops;
        end
        
         function [obj optprops_out] = addNoiseToSpecs(obj, domain, optprops, strength_noise)
  
            for i = 1:size(domain,1)
                for j = 1:size(domain, 2)
                    if domain(i,j)>0
                        noise = abs(randn(size(obj.spectra(1,:))));
                        optprops.mua(i,j,:) = squeeze(optprops.mua(i,j,:))+squeeze(strength_noise.*noise');
                        %optprops.mus(i,j,:) = mu_s(:);
                    end
                end
            end
            
            optprops_out = optprops;
        end
        
        function [obj, p0, fluence_d, fluence_d_norm] = doFluenceSimulation_hh(obj,wavelength,optprops,domain)
            
            wav_idx = find(obj.wavelengths==wavelength);
            if isempty(wav_idx)
                fprintf('wavelength does not match. Reinitialize light fluence sim object with required wavelength \n');
                fluence_d = 0;
                fluence_d_norm = 0;
                p0 = 0;
            else
            mu_s_spat = optprops.mus(:,:,wav_idx);
            mu_a_spat = optprops.mua(:,:,wav_idx);
                
                method = 'fvm';
                binning = 1;
            
                switch method
                    case 'fvm'
                        addpath('2DFiniteVolumeMethod\Fluence_Simulations\fluence\fvm\');
                        %Call fvm light propagation model
                        D=1./(3*(mu_a_spat + mu_s_spat)).*double(domain) ; %kappa
                        surface_mask = 1-imtranslate(domain,[0, 1]);
                        lightSources = getLightSourcesHH(domain, obj.voxelsize_in_cm, binning,surface_mask ) ;
                        %lightSources = getLightSources(mesh_shape, voxelsize_in_cm, binning ) ;
                        sol_image_fvm_rand = getFluenceModel( domain, obj.voxelsize_in_cm, binning, mu_a_spat, D , lightSources ) ;
                        image_mixed_fluence_spat(:,:) = sol_image_fvm_rand;
                        % image_mixed_with_fluence is th product of the simulations
                        image_mixed_with_fluence_spat(:,:) = mu_a_spat.*sol_image_fvm_rand;
                        U_d_spat(:,:) = image_mixed_fluence_spat(:,:);
                end
                
                fluence_d = U_d_spat;
                
                for i=1:size(U_d_spat,1)
                    for j=1:size(U_d_spat,2)
                        if norm(squeeze(fluence_d(i,j)))>0
                            fluence_d_norm(i,j) = fluence_d(i,j)./norm(squeeze(fluence_d(i,j)));
                        else
                            fluence_d_norm(i,j) = fluence_d(i,j);
                        end
                    end
                end
                
                p0 = image_mixed_with_fluence_spat;
            end
        end
        
        function [obj, p0_multiwav, fluence_d, fluence_d_norm] = doFluenceSimulationAllWavelength_hh(obj,optprops,domain)
            
            method = 'fvm';
            binning = 1;
            
            for w=1:length(obj.wavelengths)
                mu_s_spat = optprops.mus(:,:,w);
                mu_a_spat = optprops.mua(:,:,w);
                
                switch method
                    case 'fvm'
                        addpath('2DFiniteVolumeMethod\Fluence_Simulations\fluence\fvm\');
                        %Call fvm light propagation model
                        D=1./(3*(mu_a_spat + mu_s_spat)).*double(domain) ; %kappa
                        surface_mask = 1-imtranslate(domain,[0, 1]);
                        lightSources = getLightSourcesHH(domain, obj.voxelsize_in_cm, binning,surface_mask ) ;
                        %lightSources = getLightSources(mesh_shape, voxelsize_in_cm, binning ) ;
                        sol_image_fvm_rand = getFluenceModel( domain, obj.voxelsize_in_cm, binning, mu_a_spat, D , lightSources ) ;
                        image_mixed_fluence_spat(:,:,w) = sol_image_fvm_rand;
                        % image_mixed_with_fluence is th product of the simulations
                        image_mixed_with_fluence_spat(:,:,w) = mu_a_spat.*sol_image_fvm_rand;
                        U_d_spat(:,:,w) = image_mixed_fluence_spat(:,:,w);
                end
            end
            fluence_d = U_d_spat;
            
            for i=1:size(U_d_spat,1)
                for j=1:size(U_d_spat,2)
                    if norm(squeeze(fluence_d(i,j,:)))>0
                        fluence_d_norm(i,j,:) = fluence_d(i,j,:)./norm(squeeze(fluence_d(i,j,:)));
                    else
                        fluence_d_norm(i,j,:) = fluence_d(i,j,:);
                    end
                end
            end
            
            p0_multiwav = image_mixed_with_fluence_spat;
            
        end
        
        function [obj, mu_s] = returnRedScatteringCoeff(obj,wavelength,material)
            %wavelength in nm 
            switch material
                case 'blood'
                    a = 22; %cm^-1
                    b = 0.660; 
                case 'skin'
                    a = 30; %cm^-1
                    b = 1.1;
                case 'muscle'
                    a = 13; %cm^-1
                    b = 0.93; 
                case 'fat'
                    a = 21.6; %cm^-1
                    b = 0.930;
                otherwise
                    a = 10; %cm^-1
                    b = 0.1;
            end      
            mu_s = a.*(wavelength./500).^(-b);
        end
        
        function checkSpectralMap(obj,specmap,slice)
           h = checkSpectrum(specmap,slice);
           %setappdata(h,'mydata',mydatavalue)
        end
        
    end
    
end

