% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

classdef ClassicalField
    %ClassicalField Summary of this class goes here
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
        omega = 0;
        Pin = 0.0;
        flux = 0.0;
        tau = 0.0;
        N = 0.0;
    end

    methods
        function obj = ClassicalField(x, y, ux, uy, uz, vx, vy, vz, neff, chi2r, wl)
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
            obj.chi1i = 2*imag(neff)*(obj.neff);
            obj.chi2r = chi2r;
            obj.wl = wl;
            obj.omega = 2*pi*obj.c_const/obj.wl;
        end
        function obj = setPower(obj, Pin)
            obj.Pin = Pin;
            obj.flux = obj.Pin/(obj.hbar*obj.omega);
            if obj.tau ~= 0
                obj.N = obj.flux*obj.tau;
            end
        end
        function obj = setTime(obj, tau)
            obj.tau = tau;
            obj.N = obj.flux*tau;
        end
        function obj = normalizeField(obj, L)
            dx = obj.x(2)-obj.x(1);
            dy = obj.y(2)-obj.y(1);

            A = sqrt(sum(sum(obj.ux*conj(obj.ux) + obj.uy*conj(obj.uy) + obj.uz*conj(obj.uz)))*dx*dy);
            obj.ux = 1/sqrt(L)*1/A*obj.ux;
            obj.uy = 1/sqrt(L)*1/A*obj.uy;
            obj.uz = 1/sqrt(L)*1/A*obj.uz;
            obj.vx = 1/sqrt(L)*1/A*obj.vx;
            obj.vy = 1/sqrt(L)*1/A*obj.vy;
            obj.vz = 1/sqrt(L)*1/A*obj.vz;
        end
        function obj = normalizeClassic(obj)
            dx = obj.x(2)-obj.x(1);
            dy = obj.y(2)-obj.y(1);

            A = sqrt(sum(sum(obj.ux*conj(obj.ux) + obj.uy*conj(obj.uy) + obj.uz*conj(obj.uz)))*dx*dy);
            B = 1/sqrt(obj.c_const*obj.neff*obj.eps0);
            obj.ux = B/A*obj.ux;
            obj.uy = B/A*obj.uy;
            obj.uz = B/A*obj.uz;
            obj.vx = B/A*obj.vx;
            obj.vy = B/A*obj.vy;
            obj.vz = B/A*obj.vz;
        end
    end
end

