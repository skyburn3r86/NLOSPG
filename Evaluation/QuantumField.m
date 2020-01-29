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
        c_const = 1/(eps0*mu0)
        ux = [];
        uy = [];
        uz = [];
        vx = [];
        vy = [];
        vz = [];
        neff = 0;
        chi1i = 0;
        chi2r = [];
        wl = [];
        psi = [];
        ad = [];
        a = [];
        NPhotons = 10;
    end

    methods
        function obj = Field(ux, uy, uz, vx, vy, vz, neff, chi2r, wl)
            %FIELD Construct an instance of this class
            %   Detailed explanation goes here
            obj.ux = ux;
            obj.uy = uy;
            obj.uz = uz;
            obj.vx = vx;
            obj.vy = vy;
            obj.vz = vz;
            obj.neff = real(neff);
            obj.chi1i = 0.5*imag(sqrt(neff))/sqrt(obj.neff)
            obj.chi2r = chi2r;
            obj.wl = wl;
        end

        function [obj] = initializeState(obj)
            a = diag(sqrt(n), 1);
            ad = a';

            psi = zeros(obj.NPhotons + 1, 1);
            psi(1) = 1;
        end

    end
end

