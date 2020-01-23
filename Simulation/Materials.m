% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Materials(model,Selections)
%MATERIALS Creates the Materials for the Model.

    % Define Material Parameters
    options = struct(...
        'Au', [],...
        'SiO2', [],...
        'Si', [],...
        'SiNx', [],...
        'Organics', [],...
        'Al2O3', []);

    % Define Models for ElectroDynamics
    models = struct(...
        'n_Au', 'real(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))',...
        'k_Au', 'imag(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))',...
        'n_SiO2', 'sqrt(1+ASiO2*wl^2/(wl^2-l0SiO2^2) + BSiO2*wl^2/(wl^2-l1SiO2^2) + CSiO2*wl^2/(wl^2-l2SiO2^2))',...
        'k_SiO2', '0',...
        'n_Si', 'ASi + BSi*1/(wl^2-l0Si^2) + CSi*(1/(wl^2-l0Si^2))^2 + DSi*wl^2 + ESi*wl^4',...
        'k_Si', '0',...
        'n_SiNx', 'sqrt(1+ASiNx*wl^2/(wl^2-l0SiNx^2) + BSiNx*wl^2/(wl^2-l1SiNx^2))', ...
        'k_SiNx', '0',...
        'n_Al2O3', 'sqrt(1+AAl2O3*wl^2/(wl^2-l0Al2O3^2) + BAl2O3*wl^2/(wl^2-l1Al2O3^2) + CAl2O3*wl^2/(wl^2-l2Al2O3^2))',...
        'k_Al2O3', '0', ...
        'n_Organics', 'real(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))-0.5*real(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))^3*r33*es.Ey', ...
        'k_Organics', 'imag(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))');
    % Define Parameters for ElectroStatics
    ES = struct(...
        'eps_Al2O3', '9.8',...
        'eps_Organics', '3');

    % Read selections
    optionNames=fieldnames(options);
    selectionNames = fieldnames(Selections);
    for ii = 1:length(selectionNames)
        inpName = selectionNames{ii};
       if any(strcmp(inpName,optionNames))
          options.(inpName) = Selections.(inpName);
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    counter = 1;
    % Writes Selections to the model
    for ii = 1:length(optionNames)
        name = optionNames{ii};
        selection = options.(name);
        n = models.(['n_' name]);
        k = models.(['k_' name]);
        model.component('comp1').material.create(['mat', num2str(counter)], 'Common');
        model.component('comp1').material(['mat', num2str(counter)]).selection.set(selection);
        model.component('comp1').material(['mat', num2str(counter)]).propertyGroup.create('RefractiveIndex', 'Refractive index');
        model.component('comp1').material(['mat', num2str(counter)]).label(name);
        model.component('comp1').material(['mat', num2str(counter)]).propertyGroup('RefractiveIndex').set('n', '');
        model.component('comp1').material(['mat', num2str(counter)]).propertyGroup('RefractiveIndex').set('ki', '');
        model.component('comp1').material(['mat', num2str(counter)]).propertyGroup('RefractiveIndex').set('n', {n '0' '0' '0' n' '0' '0' '0' n});
        model.component('comp1').material(['mat', num2str(counter)]).propertyGroup('RefractiveIndex').set('ki', {k '0' '0' '0' k '0' '0' '0' k});
        if strcmp(name, 'Al2O3') || strcmp(name, 'Organics')
            eps = ES.(['eps_' name]);
            model.component('comp1').material(['mat', num2str(counter)]).propertyGroup('def').set('relpermittivity', {eps});
        end
        counter = counter + 1;
    end

end

