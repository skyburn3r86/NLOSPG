classdef SimResults
    %SIMRESULTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wl = linspace(1100, 2300, 25); 
        neff = []
        alpha = []
        hWG = 100e-9; 
        hOrganic = 100e-9; 
    end
    
    methods
        function obj = SimResults(wl, neff, alpha, hWG, hOrganic)
            obj.wl = wl; 
            obj.neff = neff; 
            obj.alpha = alpha; 
            obj.hWG = hWG; 
            obj.hOrganic = hOrganic;
            
        end
        
    end
end

