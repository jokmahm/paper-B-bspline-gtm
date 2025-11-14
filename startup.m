function startup
    format compact

    % Root folder automatically detected
    PLATEROOT = fileparts(mfilename('fullpath'));

    % Add subfolders dynamically
    folders = { ...
        'Bspline', 'Comparision', 'Data', 'DataFun', 'Fecundity', ...
        'generalFun', 'GTM', 'Maturity', 'Mortality', 'Plots', ...
        'Recruit', 'Valnerability', 'simulation', 'Kmeans', ...
        'MarkovProperties', 'Write', 'GTM_Normal', 'GTM_Gamma', ...
        'GTM_LogNormal' ...
    };

    for i = 1:length(folders)
        addpath(fullfile(PLATEROOT, 'src', folders{i}));
    end

end
