clear all;

% Select the function to be optimized. Needs to output a single FoM
func = @Master;
% Define Parameters to be optimized
parameters= struct(...
            'Names', {'r', 'hOrganic', 'wWG'}, ...
            'Units', {'',  '[nm]', '[nm]'}, ...
            'mins', [0.1, 50, 1200], ...
            'maxs', [0.95, 500, 2000], ...
            'TerminationSTD', [2, 1]);

optionNames = fieldnames(parameters);
newArgs = cell(1, length(optionNames)*2);
for ii = 1:length(optionNames)
    newArgs{2*ii -1} = optionNames{ii};
    x = {parameters(:).(optionNames{ii})};
    if isa(x{1}, 'double')
        x = x{1}; 
    end
    newArgs{2*ii} = x;
end

x = NelderMead(func, newArgs{:})