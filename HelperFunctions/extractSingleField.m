function [field] = extractSingleField(model, varargin)
    warning('off','all');

    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'OuterSolNums'
                OuterSolNums = varargin{ii+1};
            case 'SolNums'
                SolNums = varargin{ii+1};
        end
    end
    

    iter = OuterSolNums;

    temp = mpheval(model, 'ewfd.Ex', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums);
    Ex = temp.('d1');
    cordsx = temp.('p');

    temp = mpheval(model, 'ewfd.Ey', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums);
    Ey = temp.('d1');
    cordsy = temp.('p');

    temp = mpheval(model, 'ewfd.Ez', 'dataset', 'dset2', 'outersolnum', iter, 'solnum', SolNums);
    Ez = temp.('d1');
    cordsz = temp.('p');

    x = cordsx(1, :);
    y = cordsx(2, :);

    n = 400;

    x_edge=linspace(min(x),max(x),n);
    y_edge=linspace(min(y),max(y),n);
    [X,Y]=meshgrid(x_edge,y_edge);

    % Map Fields onto symmetric coordinate system.

    x = cordsx(1, :);
    y = cordsx(2, :);
    ux=griddata(x,y,Ex,X,Y);
    II = isnan(ux);
    ux(II) = 0;

    x = cordsy(1, :);
    y = cordsy(2, :);
    uy=griddata(x,y,Ey,X,Y);
    II = isnan(uy);
    uy(II) = 0;

    x = cordsz(1, :);
    y = cordsz(2, :);
    uz=griddata(x,y,Ez,X,Y);
    II = isnan(uz);
    uz(II) = 0;

    field = ClassicalField(x_edge(1,:), y_edge(1,:), ux, uy, uz, 0, 0, 0, 0, 0, 0, 1310e-9);

    warning('on','all');
end

