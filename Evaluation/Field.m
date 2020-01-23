% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
classdef Field
    %FIELD Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Ex = []
        Ey = []
        Ez = []
        omega = 0;
        neff = 0;
    end

    methods
        function obj = Field(Ex, Ey, Ez, omega, neff)
            %FIELD Construct an instance of this class
            %   Detailed explanation goes here
            obj.Ex = Ex;
            obj.Ey = Ey;
            obj.Ez = Ez;
            obj.omega = omega;
            obj.neff = neff;
        end

        function [obj] = normalizeFields(obj, x, y)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            m = length(x);
            n = length(y);
            [nx, mx] = size(obj.Ex);
            assert(nx == n && mx == m, 'Fields Dimensions do not equal the grid dimensions')
            dx = abs(x(2) - x(1));
            dy = abs(y(2) - y(1));
            I = sum(sum(dx*dy*conj(obj.Ex).*obj.Ex)) + sum(sum(dx*dy*conj(obj.Ey).*obj.Ey)) + sum(sum(dx*dy*conj(obj.Ez).*obj.Ez));
            A = 1/sqrt(abs(I));                                         % Taking Absolute, because imaginary parts will be residual numerical errors
            obj.Ex = A*obj.Ex;
            obj.Ey = A*obj.Ey;
            obj.Ez = A*obj.Ez;
        end
    end
end

