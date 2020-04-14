% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Compute(model,varargin)
%COMPUTE Starts the Computation of the Model in this function. Takes the
%model as argument and can be given the 'mat' as parameter. Si and SiNx
%implemented
    global old_neff
    options = struct(...
            'mat', 'Si', ...
            'wl', {[1310,2620], '[nm]'});

    optionNames=fieldnames(options);

    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
       if any(strcmp(inpName,optionNames))
           if isa(pair{2}, 'cell')
                options(1).(inpName) = pair{2}{1};
                options(2).(inpName) = pair{2}{2};
           else
               options(1).(inpName) = pair{2};
               options(2).(inpName) = pair{2}; 
           end

       else
          error('%s is not a recognized parameter name',inpName)
       end
    end

    % Create Study
    model.study.create('std1'); 
    % Mode Analysis
    model.study('std1').create('mode', 'ModeAnalysis');
    model.study('std1').feature('mode').set('ngen', '5');
    model.study('std1').feature('mode').activate('ewfd', true);
    model.study('std1').feature('mode').set('modeFreq', 'c_const/wl');
    model.study('std1').feature('mode').set('shiftactive', true);
    if strcmp(options(1).('mat'), 'Si')
        model.study('std1').feature('mode').set('shift', num2str(old_neff+0.3));
    elseif strcmp(options(1).('mat'), 'SiNx')
        model.study('std1').feature('mode').set('shift', num2str(old_neff+0.2));
    end
    % Generate the string of relevant wavelengths
    lambdas = ''; 
    for ii = 1:length(options(1).('wl'))
        lambdas = [lambdas num2str(options(1).('wl')(ii)) ', '];
    end
    lambdas = lambdas(1:end-2);

    model.study('std1').create('param', 'Parametric');
    model.study('std1').feature('param').setIndex('pname', 'wl', 0);
    model.study('std1').feature('param').setIndex('plistarr', lambdas, 0);
    model.study('std1').feature('param').setIndex('punit', options(2).('wl')(2:end-1), 0);

    model.sol.create('sol1');
    model.sol('sol1').study('std1');

    model.study('std1').feature('mode').set('notlistsolnum', 1);
    model.study('std1').feature('mode').set('notsolnum', '1');
    model.study('std1').feature('mode').set('listsolnum', 1);
    model.study('std1').feature('mode').set('solnum', '1');

    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').feature('st1').set('study', 'std1');
    model.sol('sol1').feature('st1').set('studystep', 'mode');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').feature('v1').set('control', 'mode');
    model.sol('sol1').create('e1', 'Eigenvalue');
    model.sol('sol1').feature('e1').set('neigs', 10);
    model.sol('sol1').feature('e1').set('shift', '1');
    model.sol('sol1').feature('e1').set('control', 'mode');
    model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
    model.sol('sol1').feature('e1').create('d1', 'Direct');
    model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'mumps');
    model.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (ewfd)');
    model.sol('sol1').attach('std1');

    model.batch.create('p1', 'Parametric');
    model.batch('p1').study('std1');
    model.batch('p1').create('so1', 'Solutionseq');
    model.batch('p1').feature('so1').set('seq', 'sol1');
    model.batch('p1').feature('so1').set('store', 'on');
    model.batch('p1').feature('so1').set('clear', 'on');
    model.batch('p1').feature('so1').set('psol', 'none');
    model.batch('p1').set('pname', {'wl'});
    model.batch('p1').set('plistarr', {lambdas});
    model.batch('p1').set('sweeptype', 'sparse');
    model.batch('p1').set('probesel', 'all');
    model.batch('p1').set('probes', {});
    model.batch('p1').set('plot', 'off');
    model.batch('p1').set('err', 'on');
    model.batch('p1').attach('std1');
    model.batch('p1').set('control', 'param');

    model.sol.create('sol2');
    model.sol('sol2').study('std1');
    model.sol('sol2').label('Parametric Solutions 1');

    model.batch('p1').feature('so1').set('psol', 'sol2');
    model.batch('p1').run;

end

