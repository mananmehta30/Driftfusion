function structCell = genIntStructs(struct_eq, struct_light, startInt, endInt, points)
%GENINTSTRUCTS - Generates a cell containing structures of solutions at various light intensities
%
% Syntax:  structCell = genIntStructs(struct_eq, struct_light, startInt, endInt, points)
%
% Inputs:
%   STRUCT_EQ - a solution struct as created by PINDRIFT in dark conditions, or a logic false for not having the dark.
%   STRUCT_LIGHT - a solution struct as created by PINDRIFT with illumination.
%   STARTINT - higher requested illumination.
%   ENDINT - lower requested illumination.
%   POINTS - number of illumination requested between STARTINT and ENDINT, including extrema, except dark.
%
% Outputs:
%   STRUCTCELL - a cell containing structs of solutions at various light
%     intensities
%
% Example:
%   genIntStructs(ssol_i_eq, ssol_i_light, 100, 0.1, 4)
%     prepare solutions at 100, 10, 1, 0.1 and 0 illumination intensities
%   genIntStructs(false, ssol_i_light, 100, 0.1, 4)
%     as above but without the solution at 0 illumination (dark)
%
% Other m-files required: changeLight
% Subfunctions: none
% MAT-files required: none
%
% See also changeLight, pindrift.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com  
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

if endInt > startInt
    warning('endInt should be smaller than startInt')
end

struct_Int = struct_light;

% define light intensity values, always decreasing
Int_array = logspace(log10(startInt), log10(endInt), points);
Int_array = unique(Int_array, 'stable'); % remove duplicates (just in case of startInt equal to endInt) and do not change ordering

% if stabilized dark solution is provided in input, calculate the point at dark
% conditions using that solution
if isa(struct_eq, 'struct') % dark calculation can be disabled using "false" as first command argument
    Int_array = [Int_array, 0];
    % decrease annoiance by figures popping up
    struct_eq.params.figson = 0;
end

% pre allocate
structCell = cell(2, length(Int_array));

existingVars = evalin('base', 'who');
changeLight_tmax = false; % automatic tmax for the first run

for i = 1:length(Int_array)
    disp([mfilename ' - illumination intensity ' num2str(Int_array(i))])
    name = matlab.lang.makeValidName([inputname(2) '_Int_' num2str(Int_array(i))]);
    if any(strcmp(existingVars, name)) % check if a structure with the same name already exists
        struct_Int = evalin('base', name);
    elseif Int_array(i) == 1 % in dark use the symstruct_eq solution, needed for no ions case where changeLight cannot reach dark
        struct_Int = struct_light;
    elseif ~Int_array(i)
        struct_Int = struct_eq;
    end
    % decrease annoiance by figures popping up
    struct_Int.params.figson = 0;
    % in any case (even the one not considered by the if), stabilize at the new intensity
    struct_Int = changeLight(struct_Int, Int_array(i), changeLight_tmax); % change light intensity
    changeLight_tmax = struct_Int.params.tmax / 2; % time to use for next iteration
    structCell{1, i} = struct_Int;
    structCell{2, i} = name;
    % restore figson before saving
    struct_Int.params.figson = 1;
    assignin('base', name, struct_Int);
end

%------------- END OF CODE --------------
