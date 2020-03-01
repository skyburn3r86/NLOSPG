% Generated by/Last edit by
% Author:           Christian Haffner/Christian Haffner
% E-Mail:           chrisitan.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [param_list] = combParameterSweep(para_sweep)   
% Funciton: generting the array containing all parameter combiniation onto linear
% list (defined by the rows of parameters_list_values    
% input parameters:
% para_sweep = structures that contains array .values and paramter .str

param_list_values = [];
    for dimenison = 1:length(para_sweep)
        param_list.str{dimenison} = (para_sweep{dimenison}.str);  
        param_list.unit{dimenison} = (para_sweep{dimenison}.unit);         
    end
    
    % parameters_list_values(idx_row,idx_colm) row are first index, column second index
    switch length(para_sweep)
        case 1
            param_lis.values = combvec(para_sweep{1}.values)';
        case 2
            param_list.values = combvec(para_sweep{1}.values, para_sweep{2}.values)';
        case 3
            param_list.values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values)';
        case 4
            param_list_values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values, para_sweep{4}.values)';
        case 5
            param_list_values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values, para_sweep{4}.values,...
                para_sweep{5}.values)';
        case 6
            param_list_values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values, para_sweep{4}.values,...
                para_sweep{5}.values, para_sweep{6}.values)';
        case 7
            param_list_values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values, para_sweep{4}.values,...
                para_sweep{5}.values, para_sweep{6}.values,...
                para_sweep{7}.values)';
        case 8
            param_list_values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values, para_sweep{4}.values,...
                para_sweep{5}.values, para_sweep{6}.values,...
                para_sweep{7}.values, para_sweep{8}.values)';
        case 9
            param_list_values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values, para_sweep{4}.values,...
                para_sweep{5}.values, para_sweep{6}.values,...
                para_sweep{7}.values, para_sweep{8}.values,...
                para_sweep{9}.values)';
        case 10
            param_list_values = combvec(para_sweep{1}.values, para_sweep{2}.values,...
                para_sweep{3}.values, para_sweep{4}.values,...
                para_sweep{5}.values, para_sweep{6}.values,...
                para_sweep{7}.values, para_sweep{8}.values,...
                para_sweep{9}.values, para_sweep{10}.values)';
        otherwise
            error('Addapt the function combParameterSweep to handle the desired number of parameters');
    end
            
    for idx_row = 1:size(param_list.values,1)
        param_list.print{idx_row} = [];
        for idx_colm = 1:length(para_sweep)
             param_list.print{idx_row}  = [param_list.print{idx_row}...
               ' >_< ' strrep(param_list.str{idx_colm},'_array', '') ...
               ' ' strrep(num2str(param_list.values(idx_row,idx_colm)),'.' , 'pt') ...
               '' param_list.unit{idx_colm}];
        end
%         param_list_print{idx_row}; for testing
    end
        
end