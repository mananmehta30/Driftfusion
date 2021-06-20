function EA_script_plot(type, dir_file_name, varargin)
%EA_SCRIPT_PLOT - Plots a list of Stark spectroscopy a.k.a. electroabsorbance (EA) results with explicit line makeup
% The plots can be saved to files if the dir_file_name option is provided
%
% Syntax:  IS_script_plot(type, dir_file_name, varargin)
%
% Inputs:
%   TYPE - char array, which quantity to plot, can be either
%     '1h', '2h' or 'phase'
%   DIR_FILE_NAME - char array, images with this prefix will be created in
%     a directory with this name. To avoid the saving, an empty char array
%     can be specified: ''.
%   VARARGIN - many arguments, the script will group the arguments in
%     groups of three:
%       the first of each group as the EA struct with the IS simulation
%         data;
%       the second of each group as the legend entry;
%       the third, which has to be a {cell}, as the options to be passed to
%         the plot command (e.g. specifying the color of the line).
%
% Example:
%   EA_script_plot('2h', '20210409_spiro',...
%       EA_spiro_10000x_sc_800mV_1sun, '10000x 1 sun',{':r','LineWidth',3},...
%       EA_spiro_1000x_sc_800mV_1sun, '1000x 1 sun',{'-r'},...
%       EA_spiro_sc_800mV_1sun, ' 1 sun',{'--r'})
%     plot 3 different EA simulations, all in red but with different line
%     style and different legends. And save the graphics as files.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also EA_script.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

fig = figure('Name', 'Amplitude of EA second harmonic E_{AC}^2', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
hold off
%ymin = Inf;
%ymax = -Inf;
for i = 1:(length(varargin)/3)
    isol = i*3-2;
    legendarr{i} = varargin{isol+1};
    options = varargin{isol+2};
    switch type
        case '1h'
            plot(varargin{isol}.Freq, varargin{isol}.AC_ExDC_E_amp, options{:})
        case '2h'
            plot(varargin{isol}.Freq, varargin{isol}.AC_Efield2_amp, options{:})
        case 'phase'
            phase_n_deg = rad2deg(wrapTo2Pi(varargin{isol}.AC_ExDC_E_phase));
            plot(varargin{isol}.Freq, phase_n_deg, options{:})
        otherwise
            error([mfilename(1) ' - *type* option needed'])
    end
    %ymin = min(ymin, min());
    %ymax = max(ymax, max());
    %range = ymax-ymin;
    %ylim([ymin-0.03*range, ymax+0.03*range])
    %xlim([min(min(EA_results.Freq)), max(max(EA_results.Freq))])

    hold on
end

ax = gca;
ax.XScale = 'log'; % for putting the scale in log
ax.YScale = 'log'; % for putting the scale in log
xlabel('Frequency [Hz]');
ylabel('Abs(E_{AC}^2) [V^2/cm^2]');

legend(legendarr, 'Location', 'southeast')
legend boxoff

valid_name = true;
try
    java.io.File(dir_file_name).toPath;
catch
    valid_name = false;
end
if ~isempty(dir_file_name) && valid_name
    if 7~=exist(dir_file_name,'dir')
        mkdir(dir_file_name)
    end
    saveas(fig, [char(dir_file_name) filesep char(dir_file_name) char(['-EA_' type '.fig'])])
    saveas(fig, [char(dir_file_name) filesep char(dir_file_name) char(['-EA_' type '.png'])])
end

end

%------------- END OF CODE --------------