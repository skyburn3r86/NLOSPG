% Generated through IntelliJ IDEA
% Author:           Michael Doderer (Killian Keller)
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
% Last Edit:        26.03.2020 (Killian)
% Code commented and adapted by Killian, Original Version by El Dodo
%
% Source: The amazing O&P Script from the IEF group, where the original source is not specified. Probably Saleh.
% Input: model and oldField.
% Returns the Classical Field Struct of the evaluated Field, and further the overlap.

% TODO: Implement the stretching Factor for the adaptation of the field.

function [overlap, field] = overlapFields(model, oldField, varargin)

    warning('off', 'all')

    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('Arguments needs propertyName/propertyValue pairs')
    end

    dset = 'dset1';
    solnum = 1;
    outersolnum = 1;
    N = 400;
    % transfroming variable inputs into variables
    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'dset'
                dset = varargin{ii + 1};
            case 'solnum'
                solnum = varargin{ii + 1};
            case 'outersolnum'
                outersolnum = varargin{ii + 1};
            case 'N'
                N = 400;
        end
    end

    % creating the meshgrid to interpolate the mph model onto. This is similar to using mpheval, and then mapping onto
    % a regular meshgrid. This seems to be less time consuming, More detailed analysis should be made.

    % loads information about coordinate system from the simulation
    tempx = mpheval(model, 'x', 'dataset', dset, 'outersolnum', outersolnum, 'solnum', solnum);
    x = tempx.p(1, :);
    y = tempx.p(2, :);

    % Creates regular Coordinate System
    x_edge = linspace(0, max(x), N);  % x Coordinates
    y_edge = linspace(0, max(y), N);  % y coordinates
    [X,Y] = meshgrid(x_edge,y_edge);
    coord = [X(:),Y(:)]';
    
    % Create origin and target field structs
    t = struct();
    o = struct(); 
    
    % Loads Fields from Simulation. Unfortunately, returns one single row vector. . Remapping using matlab functions
    [t.Ex, t.Ey, t.Hx, t.Hy, t.Ez, t.Hz] = mphinterp(model ,{'ewfd.Ex','ewfd.Ey','ewfd.Hx','ewfd.Hy','ewfd.Ez','ewfd.Hz'}, ...
    'coord',coord,'solnum',solnum,'outersolnum',outersolnum);
    o.Ex = reshape(oldField.ux, [1, N^2]); 
    o.Ey = reshape(oldField.uy, [1, N^2]);
    o.Ez = reshape(oldField.uz, [1, N^2]);
    o.Hx = reshape(oldField.vx, [1, N^2]);
    o.Hy = reshape(oldField.vy, [1, N^2]);
    o.Hz = reshape(oldField.vz, [1, N^2]);
    
    % Calculate normalization, i.e. the total input power from both the input mode and the tb Coupled mode.
    Ptot = 1/2 * sum(real(o.Ex.*conj(o.Hy) - o.Ey.*conj(o.Hx)));
    PTarget = 1/2 * sum(real(t.Ex.*conj(t.Hy) - t.Ey.*conj(t.Hx)));

    % Calculate the overlap
    overlap = abs(1/4*sum(o.Ex.*conj(t.Hy) - o.Ey.*conj(t.Hx) + o.Hy.*conj(t.Ex) - o.Hx.*conj(t.Ey))) .^2;
    overlap = overlap/Ptot/PTarget;

    % Return the classicalfield
    field = ClassicalField(x_edge, y_edge, t.Ex, t.Ey, t.Ez, t.Hx, t.Hy, t.Hz, 0, 0, 0, 1310e-9);
    end