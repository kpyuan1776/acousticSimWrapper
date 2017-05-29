classdef DeviceInfo
    %DEVICEINFO class contains information about the device geometry for
    %acoustic reconstruction
    %   create object: dev = deviceInfo('sim_mousecart')
    
    
    
    properties
        n;
        samples;
        proj;
        time_res;
        r_sensor;
        angle_sensor;
        c0;
        limits;
        imp_resp;
        fsample;
        delayNumber;
        image_width;
        devname;
    end
    
    methods
        function Dobj = DeviceInfo(device)
            if strcmp(device,'handheld_acuity')
                Dobj.n = 100;
                Dobj.samples = 2030;
                Dobj.proj = 256;
                Dobj.time_res = 2;
                Dobj.r_sensor = .061;
                Dobj.angle_sensor = ( -4 : -172/255 : -176 ) * pi/180;
                Dobj.c0 = 1500;
                Dobj.limits = Dobj.r_sensor + [ -.03 .03 ] ;
                Dobj.imp_resp = zeros( Dobj.samples, 1 ) ;
                Dobj.imp_resp( 1015 ) = 1 ;
                Dobj.delayNumber = 460;
                Dobj.image_width = 60e-3;
                Dobj.fsample = 4e7;
                Dobj.devname = 'handheld_acuity';
            elseif strcmp(device,'handheld_eip')
                Dobj.n = 100;
                Dobj.samples = 2030;
                Dobj.proj = 256;
                Dobj.time_res = 2;
                Dobj.r_sensor = .061;
                Dobj.angle_sensor = ( -4 : -172/255 : -176 ) * pi/180;
                Dobj.c0 = 1500;
                Dobj.limits = Dobj.r_sensor + [ -.03 .03 ] ;
                Dobj.imp_resp = zeros( Dobj.samples, 1 ) ;
                Dobj.imp_resp( 1015 ) = 1 ;
                Dobj.delayNumber = 560;
                Dobj.image_width = 60e-3;
                Dobj.fsample = 4e7;
                Dobj.devname = 'handheld_eip';
            elseif strcmp(device,'handheld_acuityNew')
                Dobj.n = 100;
                Dobj.samples = 2030;
                Dobj.proj = 256;
                Dobj.time_res = 2;
                Dobj.r_sensor = .060;
                Dobj.angle_sensor = ( -17.5 : -145/255 : -162.5 ) * pi/180;
                Dobj.c0 = 1500;
                Dobj.limits = Dobj.r_sensor + [ -.03 .03 ] ;
                Dobj.imp_resp = zeros( Dobj.samples, 1 ) ;
                Dobj.imp_resp( 1015 ) = 1 ;
                Dobj.delayNumber = 460;
                Dobj.image_width = 60e-3;
                Dobj.fsample = 4e7;
                Dobj.devname = 'handheld_acuityNew';
            elseif strcmp(device,'handheld_muenster')
                Dobj.n = 100;
                Dobj.samples = 2030;
                Dobj.proj = 256;
                Dobj.time_res = 2;
                Dobj.r_sensor = .0405;
                %Dobj.angle_sensor = ( -17.5 : -145/255 : -162.5 ) * pi/180;
                Dobj.angle_sensor = -2.657366:0.00852211548825356:-0.48;%( -27.74 : -124.02/255 : -151.76 ) * pi/180;
                Dobj.c0 = 1500;
                Dobj.limits = Dobj.r_sensor + [ -.02 .02 ] ;
                Dobj.imp_resp = zeros( Dobj.samples, 1 ) ;
                Dobj.imp_resp( 1015 ) = 1 ;
                Dobj.delayNumber = 0.0;
                Dobj.image_width = 60e-3;
                Dobj.fsample = 4e7;
                Dobj.devname = 'handheld_muenster';
            elseif strcmp(device,'sim_mousecart')
                Dobj.n = 100;
                Dobj.samples = 2030;
                Dobj.proj = 256;
                Dobj.time_res = 2;
                Dobj.r_sensor = 15e-3;
                ang_ini = -0.7853981633;
                ang_end = 3.926990817079;
                ang_step = 0.0184799567858;
                Dobj.angle_sensor = [ang_ini:ang_step:ang_end];
                Dobj.c0 = 1530;
                Dobj.limits = Dobj.r_sensor + [ -.02 .02 ] ;
                Dobj.imp_resp = zeros( Dobj.samples, 1 ) ;
                Dobj.imp_resp( 1015 ) = 1 ;
                Dobj.delayNumber = 0.0;
                Dobj.image_width = 60e-3;
                Dobj.fsample = 4e7;
                Dobj.devname = 'sim_mousecart';
            elseif strcmp(device,'sim_handheldAcuity')
                Dobj.n = 150;
                Dobj.samples = 2030;
                Dobj.proj = 256;
                Dobj.time_res = 2;
                Dobj.r_sensor = .015;
                Dobj.angle_sensor = ( -17.5 : -145/255 : -162.5 ) * pi/180;
                Dobj.c0 = 1500;
                Dobj.limits = Dobj.r_sensor + [ -.0075 .0075 ] ;
                Dobj.imp_resp = zeros( Dobj.samples, 1 ) ;
                Dobj.imp_resp( 1015 ) = 1 ;
                Dobj.delayNumber = 0;
                Dobj.image_width = 10e-3;
                Dobj.fsample = 4e7;
                Dobj.devname = 'sim_handheldAcuity';
            else
                fprintf('device unknown');
            end
        end
        
        function rparams = getParameters(obj)
            rparams = struct( 'n', obj.n, ...
                'proj', obj.proj, ...
                'time_res', obj.time_res, ...
                'r_sensor', obj.r_sensor , ...
                'angle_sensor', obj.angle_sensor, ...
                'c', obj.c0, ...
                'limits', obj.limits, ...
                'imp_resp', obj.imp_resp, ...
                'fs', obj.fsample) ;
        end
        
    end
    
end
