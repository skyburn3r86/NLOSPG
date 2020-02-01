% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

classdef QuantumField
    %QuantumField Summary of this class goes here
    %   Detailed explanation goes here
    properties
        hbar = 1.0545718e-34;
        eps0 = 8.854e-12;
        mu0 = 4*pi*1e-7;
        c_const = 0;
        x = [];
        y = [];
        ux = [];
        uy = [];
        uz = [];
        vx = [];
        vy = [];
        vz = [];
        neff = 0;
        chi1i = 0;
        chi2r = [];
        wl = 0;
        omega = 0
        psi = [];
        ad = [];
        a = [];
        NPhotons = 10;
    end

    methods
        function obj = QuantumField(x, y, ux, uy, uz, vx, vy, vz, neff, chi2r, wl)
            %FIELD Construct an instance of this class
            %   Detailed explanation goes here
            obj.c_const = 1/sqrt(obj.eps0*obj.mu0); 
            obj.x = x; 
            obj.y = y; 
            obj.ux = ux;
            obj.uy = uy;
            obj.uz = uz;
            obj.vx = vx;
            obj.vy = vy;
            obj.vz = vz;
            obj.neff = real(neff);
            obj.chi1i = 0.5*imag(sqrt(neff))/sqrt(obj.neff);
            obj.chi2r = chi2r;
            obj.wl = wl;
            obj.omega = 2*pi*obj.c_const/obj.wl;
            
        end

        function [obj] = initializeState(obj, psiI)
            n = linspace(1, obj.NPhotons, obj.NPhotons); 
            obj.a = diag(sqrt(n), 1);
            obj.ad = obj.a';
            
            obj.psi = zeros(obj.NPhotons + 1, 1);
            obj.psi(psiI) = 1;
            
            % Make matrices sparse
            obj.a = sparse(obj.a); 
            obj.ad = sparse(obj.ad); 
            obj.psi = sparse(obj.psi); 
        end

        function obj = normalizeField(obj)
            dx = obj.x(2)-obj.x(1);
            dy = obj.y(2)-obj.y(1);
            A = sqrt(sum(sum(obj.ux*conj(obj.ux) + obj.uy*conj(obj.uy) + obj.uz*conj(obj.uz)))*dx*dy);
            obj.ux = 1/A*obj.ux;
            obj.uy = 1/A*obj.uy;
            obj.uz = 1/A*obj.uz;
            obj.vx = 1/A*obj.vx;
            obj.vy = 1/A*obj.vy;
            obj.vz = 1/A*obj.vz;
        end

    end
end

