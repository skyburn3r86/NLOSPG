% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Compute(model)
%COMPUTE Starts the Computation of the Model in this function. Takes the
%model as argument.

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
    
    model.study('std1').feature('mode').setIndex('activate', false, 3);
    % defines the value around which solutions are searched
    model.study('std1').feature('mode').set('shift', 'n_start');
    % defines the number of solutions searched for    
    model.study('std1').feature('mode').set('neigsactive', true);
    str_read = char(model.param.get('nr_modes'));
    try
        nr_modes = str2num(str_read(1:strfind(str_read, '[')-2));
        model.study('std1').feature('mode').set('neigs', nr_modes);
    catch
        nr_modes = str2num(str_read);
        model.study('std1').feature('mode').set('neigs', nr_modes);
    end
    %
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
    model.sol('sol1').feature('e1').set('neigs', nr_modes);
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
    
    model.sol('sol1').runAll;
end

