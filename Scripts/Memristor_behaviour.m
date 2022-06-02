%% Simulate Ag and Iodine reaction by varying vacancy concentration at the interface and introduce
% a surface recombination velocity


%% Define memristor
par_memristor = pc('Input_files/memristor.csv');
par_memristor_interface_reactions = par_memristor;
%% Get Equilbrium solutions
soleq_memristor = equilibrate(par_memristor);
soleq_par_memristor_interface_reactions = equilibrate(par_memristor_interface_reactions);


