% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Compute(model,varargin)
%COMPUTE Starts the Computation of the Model in this function. Takes the
%model as argument and can be given the 'mat' as parameter. Si and SiNx
%implemented

    options = struct(...
            'mat', 'Si');

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
           end

       else
          error('%s is not a recognized parameter name',inpName)
       end
    end

    % Create Study
    % Electro-Statics
    model.study.create('std1');
    model.study('std1').setGenConv(true);
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').activate('ewfd', false);
    model.study('std1').feature('stat').activate('es', true);
    model.study('std1').feature('stat').setIndex('activate', false, 1);

    % Mode Analysis
    model.study('std1').create('mode', 'ModeAnalysis');
    model.study('std1').feature('mode').set('ngen', '5');
    model.study('std1').feature('mode').activate('ewfd', true);
    model.study('std1').feature('mode').set('modeFreq', 'c_const/wl');
    model.study('std1').feature('mode').set('shiftactive', true);
    if strcmp(options.('mat'), 'Si')
        model.study('std1').feature('mode').set('shift', '3.2');
    elseif strcmp(options.('mat'), 'SiNx')
        model.study('std1').feature('mode').set('shift', '2.0');
    end

    % Create Solution
    model.sol.create('sol1');
    model.sol('sol1').study('std1');

    model.study('std1').feature('stat').set('notlistsolnum', 1);
    model.study('std1').feature('stat').set('notsolnum', '1');
    model.study('std1').feature('stat').set('listsolnum', 1);
    model.study('std1').feature('stat').set('solnum', '1');
    model.study('std1').feature('mode').set('notlistsolnum', 1);
    model.study('std1').feature('mode').set('notsolnum', 'auto');
    model.study('std1').feature('mode').set('listsolnum', 1);
    model.study('std1').feature('mode').set('solnum', 'auto');

    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').feature('st1').set('study', 'std1');
    model.sol('sol1').feature('st1').set('studystep', 'stat');
    % Create Variables for Static Solution
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').feature('v1').set('control', 'stat');
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
    model.sol('sol1').feature('s1').feature.remove('fcDef');
    model.sol('sol1').create('su1', 'StoreSolution');
    model.sol('sol1').create('st2', 'StudyStep');
    model.sol('sol1').feature('st2').set('study', 'std1');
    model.sol('sol1').feature('st2').set('studystep', 'mode');

    % Create Variables for Mode Solution
    model.sol('sol1').create('v2', 'Variables');
    model.sol('sol1').feature('v2').set('initmethod', 'sol');
    model.sol('sol1').feature('v2').set('initsol', 'sol1');
    model.sol('sol1').feature('v2').set('initsoluse', 'su1');
    model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
    model.sol('sol1').feature('v2').set('notsol', 'sol1');
    model.sol('sol1').feature('v2').set('control', 'mode');

    % Create Eigenvalue thingy
    model.sol('sol1').create('e1', 'Eigenvalue');
    model.sol('sol1').feature('e1').set('neigs', 6);
    model.sol('sol1').feature('e1').set('shift', '1');
    model.sol('sol1').feature('e1').set('control', 'mode');
    model.sol('sol1').feature('e1').set('linpmethod', 'sol');
    model.sol('sol1').feature('e1').set('linpsol', 'sol1');
    model.sol('sol1').feature('e1').set('linpsoluse', 'su1');
    model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
    model.sol('sol1').feature('e1').create('d1', 'Direct');
    model.sol('sol1').feature('v2').set('notsolnum', 'auto');
    model.sol('sol1').feature('v2').set('notsolvertype', 'solnum');
    model.sol('sol1').feature('v2').set('notlistsolnum', {'1'});
    model.sol('sol1').feature('v2').set('notsolnum', 'auto');
    model.sol('sol1').feature('v2').set('notlistsolnum', {'1'});
    model.sol('sol1').feature('v2').set('notsolnum', 'auto');
    model.sol('sol1').feature('v2').set('control', 'mode');
    model.sol('sol1').attach('std1');

    model.study('std1').create('param', 'Parametric');
    model.study('std1').feature('param').setIndex('plistarr', '', 0);
    model.study('std1').feature('param').setIndex('punit', 'm', 0);
    model.study('std1').feature('param').setIndex('pname', 'wl', 0);
    model.study('std1').feature('param').setIndex('plistarr', 'range(lmin,(lmax-lmin)/(Nl-1),lmax)', 0);

    model.batch.create('p1', 'Parametric');
    model.batch('p1').study('std1');
    model.batch('p1').create('so1', 'Solutionseq');
    model.batch('p1').feature('so1').set('seq', 'sol1');
    model.batch('p1').feature('so1').set('store', 'on');
    model.batch('p1').feature('so1').set('clear', 'on');
    model.batch('p1').feature('so1').set('psol', 'none');
    model.batch('p1').set('pname', {'wl'});
    model.batch('p1').set('plistarr', {'range(lmin,(lmax-lmin)/(Nl-1),lmax)'});
    model.batch('p1').set('sweeptype', 'sparse');
    model.batch('p1').set('probesel', 'all');
    model.batch('p1').set('probes', {});
    model.batch('p1').set('plot', 'off');
    model.batch('p1').set('err', 'on');
    model.batch('p1').attach('std1');
    model.batch('p1').set('control', 'param');

    model.sol.create('sol3');
    model.sol('sol3').study('std1');
    model.sol('sol3').label('Parametric Solutions 1');
    model.batch('p1').feature('so1').set('psol', 'sol3');
    model.batch('p1').run;
end

