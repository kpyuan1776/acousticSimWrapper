function [spectra] = LoadSpectra(data_path, wavelengths, varargin) 
% Load the spectra of hemoglobin and possible molecular agents from a
% predefined data path. The data path should contain a file called
% hb_spectra (hemoglobin spectra from Scott Prahl measurements) and
% possibly the spectrum of a target agent which is named: "agent_..."

% If the excitation wavelengths are given as a parameter then the spectra
% interpolated for the specific wavelengths (spline interpolation) are
% returned. Otherwise the excitation wavelengths are loaded from the same
% folder.

% Small modification: varargin{1} = fat spectrum and, varargin{2} = water

keep_path = pwd;

cd(data_path);
if (nargin<2)
    load('wavelengths.mat');
end

 nb_wavelengths = length(wavelengths);
 generaldirectory1 = dir ('*agent_*');
 generaldirectory2 = dir ('*hb_spectra*');
 
 if size(generaldirectory2,1)
    Nb_spectra = size(generaldirectory1,1) + 2;
 else
    Nb_spectra = size(generaldirectory1,1);
 end
 
 switch nargin
     case 3
         Nb_spectra = Nb_spectra +1;
     case 4
         Nb_spectra = Nb_spectra +2;
 end
    
 spectra = ones(Nb_spectra,nb_wavelengths);
 k=1;
 %load all the agents in the folder
 for i=1:size(generaldirectory1,1)
     load(generaldirectory1(i).name);
     var = eval(generaldirectory1(i).name(7:end-4));
     spectra(k,:) = spline(var(:,1),var(:,2),wavelengths);
     k=k+1;
 end
  
 %Load oxy - deoxy
 if size(generaldirectory2,1)
        load(generaldirectory2(1).name);
        
        %Callibrate, at 800 nm Hb should be equal to 1.
      
        hg_spectra(:,2) = hg_spectra(:,2)./761.72;
        hg_spectra(:,3) = hg_spectra(:,3)./761.72;
        
        %Resample
        spectra(k,:) = spline(hg_spectra(:,1),hg_spectra(:,3),wavelengths);
        spectra(k+1,:) = spline(hg_spectra(:,1),hg_spectra(:,2),wavelengths);
        clear hg_spectra;
 end

 switch nargin
     case 3
         if (strcmp(varargin{1},'fat_spectra.mat'))
            load([data_path '/' varargin{1}]);
            idxNormwav = find(fat_spectra(:,1) == 800);
            fat_spectra(:,2) = fat_spectra(:,2)./fat_spectra(idxNormwav,2);
            spectra(k+2,:) = spline(fat_spectra(:,1),fat_spectra(:,2),wavelengths);
         elseif (strcmp(varargin{1},'water_spectra.mat'))
            load([data_path '/' varargin{1}]);
            idxNormwav = find(water_spectra(:,1) == 800);
            water_spectra(:,2) = water_spectra(:,2)./water_spectra(idxNormwav,2);
            spectra(k+2,:) = spline(water_spectra(:,1),water_spectra(:,2),wavelengths);
         end
     case 4
         load([data_path '/' varargin{1}]);
         load([data_path '/' varargin{2}]);
         idxNormwav = find(fat_spectra(:,1) == 800);
         fat_spectra(:,2) = fat_spectra(:,2)./fat_spectra(idxNormwav,2);
         spectra(k+2,:) = spline(fat_spectra(:,1),fat_spectra(:,2),wavelengths);
         idxNormwav = find(water_spectra(:,1) == 800);
         water_spectra(:,2) = water_spectra(:,2)./water_spectra(idxNormwav,2);
         spectra(k+3,:) = spline(water_spectra(:,1),water_spectra(:,2),wavelengths);
 end
 
 cd(keep_path)
 
end