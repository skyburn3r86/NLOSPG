% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
%
% NELDERMEAD: Nelder-Mead Optimization. NelderMead(@function, varargin)
%   Takes as argument the function to optimize. The function must return the value to minimize.
%
%   Variable Input Arguments: Names, Units, mins and maxs. Order must be exact. Names and Units are Cell arrays, mins and maxs must be arrays
%
%   Example to Optimize a function:
%   x = NelderMead(@FunctionToOptimize, 'Names', {'height', 'width', 'ratio'}, 'Units', {'[nm]', '[nm]', ''} 'mins', [0, 0, 0], 'maxs', [1e6, 1e6, 1e6])
%
%   Additional optional arguments:
%   - alpha, gamma, rho, sigma parameters for the optimization.
%   - TerminationSTD: Termination of Algorithm is determined by the Standard Deviation of the arguments, if converging,
%   the function arguments will have a small std. Argument consists of an array, first entry is the maximal STD for the parameters,
%   The second entry is the mean STD for all the parameters. If the std of the arguments is smaller than the given values, algorithm terminates.
%
%   See https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method for more details about these parameters and the Optimization

function [x] = NelderMead(func, varargin)
    % Initial Values for the Optimization
    options = struct(...
            'alpha', 1, ...
            'gamma', 2, ...
            'rho', 0.5, ...
            'sigma', 0.5, ...
            'Names', {'r', 'hWG', 'hOrganic', 'wWG'}, ...
            'Units', {'', '[nm]', '[nm]', '[nm]'}, ...
            'mins', [0.1, 300, 50, 500], ...
            'maxs', [0.9, 500, 200, 1000], ...
            'TerminationSTD', [2, 1]);
    % ** Parse the Input varargins.

    optionNames=fieldnames(options);

    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
       if any(strcmp(inpName,optionNames))
          options(1).(inpName) = pair{2}{1};
          options(2).(inpName) = pair{2}{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    % Done Parsing

    % define parameters for Nelder Mead
    alpha = options.alpha;
    gamma = options.gamma;
    rho = options.rho;
    sigma = options.sigma;

    % Define Bounds for Nelder Mead
    mins = options.mins;
    maxs = options.maxs;

    % Define startpoints, random Points within acceptance Boundaries
    Names = {options(:).Names};
    Units = {options(:).Units};
    N = length(Names);                  % Number of Variables. For nelder mead, N + 1 startpoints are required.
    x = zeros(N+1, N);                  % Start Vector
    for ii = 1:N+1
        for kk = 1:N
            x(ii, kk) = rand()*(maxs(kk)-mins(kk)) + mins(kk);
        end
    end

    % Calculate the first N + 1 iterations
    y = zeros(N+1, 1);
    for ii = 1:N+1
        args = cell(2*N, 1);        % Creating varargin for the target functions
        for kk = 1:N
            args{2*kk-1} = Names{kk};
            args{2*kk} = {x(ii, kk), Units{kk}};
        end

        dk = func(args{:});                      % Makes the Simulation.
        p = 0;
        for kk = 1:N
            p = p + punish(x(ii, kk), mins(kk), maxs(kk));
        end
        y(ii) = dk + p;
    end
    yBucket = y;
    term = options.TerminationSTD
    while true
        % 1. Order accoring to the values at the vertices
        [~, I] = sort(y);
        x = x(I, :);
        y = y(I);

        % Termination Condition, Check spread of Starting variables, i.e.
        % Convergence to a point
        stdV = zeros(N, 1);
        for ii = 1:N
           stdV(ii) = std(x(:, ii));
        end
        if max(stdV) < term(1) && mean(stdV) < term(2)
            x = x(1, :);
            disp(num2str(length(yBucket)));
            disp(num2str(y(1)));
            break;
        end

        % 2. Calculate Centroid
        xo = zeros(1, N);
        for ii = 1:N
            for kk = 1:N
                xo(ii) = xo(ii) + x(kk, ii);
            end
            xo(ii) = xo(ii)/N;
        end

        % 3. Calculate the reflection
        xr = zeros(1, N);
        for ii = 1:N
            xrTemp = xo(ii) + alpha*(xo(ii) - x(N+1, ii));
            alphaTemp = alpha;
            while xrTemp <= 0
                alphaTemp = alphaTemp/2;
                xrTemp = xo(ii) + alphaTemp*(xo(ii) - x(N+1, ii));
            end
            xr(ii) = xrTemp;
        end

        args = cell(2*N, 1);            % Create New args
        for kk = 1:N
            args{2*kk-1} = Names{kk};
            args{2*kk} = {xr(kk), Units{kk}};
        end

        dk = func(args{:});
        p = 0;
        for kk = 1:N
            p = p + punish(xr(kk), mins(kk), maxs(kk));
        end
        yr = dk + p;

        if yr < y(N) && yr > y(1)              % Check for Performance of new point
            y(N+1) = yr;
            yBucket(end + 1) = yr;
            x(N+1, :) = xr;
            continue;                           % If yr lies in between the second worst and best, go to 1
        elseif yr < y(1)                        % Check if performed best
            % 4. Expansion
            xe = zeros(1, N);
            for ii = 1:N
                xeTemp = xo(ii) + gamma*(xr(ii) - xo(ii));
                gammaTemp = gamma;
                while xeTemp <= 0
                    gammaTemp = gammaTemp/2;
                    xeTemp = xo(ii) + gammaTemp*(xr(ii) - xo(ii));
                end
                xe(ii) = xeTemp;
            end

            args = cell(2*N, 1);            % Create New args
            for kk = 1:N
                args{2*kk-1} = Names{kk};
                args{2*kk} = {xe(kk), Units{kk}};
            end

            dk = func(args{:});
            p = 0;
            for kk = 1:N
                p = p + punish(xe(kk), mins(kk), maxs(kk));
            end
            ye = dk + p;
            if ye < yr
                y(N + 1) = ye;
                x(N + 1, :) = xe;
                yBucket(end + 1) = ye;
                continue;                       % if expanded point better than reflected point, replcae and go to 1
            else
                y(N + 1) = yr;
                yBucket(end + 1) = yr;
                x(N + 1, :) = xr;
                continue;                        % if reflected  point better than reflected point, replcae and go to 1
            end
        elseif yr >= y(N)
            % 5. Contraction
            xc = zeros(1, N);
            for ii = 1:N
                xcTemp = xo(ii) + rho*(x(N+1, ii) - xo(ii));
                rhoTemp = rho;
                while xcTemp <= 0
                    rhoTemp = rhoTemp/2;
                    xcTemp = xo(ii) + rho*(x(N+1, ii) - xo(ii));
                end
                xc(ii) = xcTemp;
            end

            args = cell(2*N, 1);            % Create New args
            for kk = 1:N
                args{2*kk-1} = Names{kk};
                args{2*kk} = {xc(kk), Units{kk}};
            end

            dk = func(args{:});
            p = 0;
            for kk = 1:N
                p = p + punish(xc(kk), mins(kk), maxs(kk));
            end
            yc = dk + p;
            if yc < y(N+1)
                y(N+1) = yc;
                yBucket(end + 1) = yc;
                x(N+1, :) = xc;
                continue;
            end
        end

        for kk = 2:N+1
            for ll = 1:N
                xTemp = x(1, ll) + sigma*(x(1, ll) - x(kk, ll));
                sigmaTemp = sigma;
                if xTemp <= 0
                    sigmaTemp = sigmaTemp/2;
                    xTemp = x(1, ll) + sigmaTemp*(x(1, ll) - x(kk, ll));
                end
                x(kk, ll) = xTemp;
            end
        end
        for ii = 2:N+1
            args = cell(2*N, 1);
            for kk = 1:N
                args{2*kk-1} = Names{kk};
                args{2*kk} = {x(ii, kk), Units{kk}};
            end
            dk = func(args{:});                      % Makes the Simulation, k1: Pump wavevector, k2: Signal wavevector
            p = 0;
            for kk = 1:N
                p = p + punish(x(ii, kk), mins(kk), maxs(kk));
            end
            y(ii) = dk + p;
            yBucket(end + 1) = y(ii);
        end
    end
end


function y = punish(x, xmin, xmax)
    k = 1e4;
    L = 1e10;
    y = (1/(1+exp(-k*(x-xmax))) + 1/(1+exp(k*(x-xmin))))*L;
end
