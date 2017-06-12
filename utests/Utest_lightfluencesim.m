classdef Utest_lightfluencesim < matlab.unittest.TestCase
    % UTEST_LIGHTFLUENCESIM checks light fluence simulations
    %   Detailed explanation goes here
    
    properties (TestParameter)
        sizeP_0_Methods = struct('small',150,'medium',300);
        
    end
    
    methods (Test)
        %test if preprocessing methods work
        function test_simulateSingleWavelengthRnd(testCase,sizeP_0_Methods)
            
            addpath('../');
            addpath('../fluenceSimulation');
            
            dev = DeviceInfo('sim_handheldAcuity');
            dev.angle_sensor = (( -17.5 : -145/255 : -162.5 ))* pi/180;
            lightsim = LightFluenceSim();
            lightsim.n = sizeP_0_Methods;
            lightsim.datapath = '../spectra';
            lightsim = lightsim.loadSpectra();
            
            [lightsim, p0] = lightsim.runSimulationRndsO2_hh();
            
            fprintf('simulate p0 size: %i \n',sizeP_0_Methods);
            verifyEqual(testCase,size(p0,1),lightsim.n);
        end
        
        
    end
    
end
