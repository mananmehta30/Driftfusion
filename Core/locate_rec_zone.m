function rec_zone_loc = locate_rec_zone(par)
% A function to locate the recombination zone within the interfacial
% regions based on the minority carrier densities at equilibrium
% intelligent location
int_logical = strcmp(par.layer_type, 'interface');
loc = find(int_logical); % interface layer locations
rec_zone_loc = strings(1, length(par.stack));

for j = 1:length(loc)

i = loc(j);
% Gradient coefficients for surface recombination equivalence
alpha0 = ((par.EA(i-1) - par.EA(i+1))/(par.kB*par.T) + (log(par.Nc(i+1))-log(par.Nc(i-1))))/par.d(i);
beta0 = ((par.IP(i+1) - par.IP(i-1))/(par.kB*par.T) + (log(par.Nv(i+1))-log(par.Nv(i-1))))/par.d(i);
    
if (alpha0 <= 0 && beta0 > 0) || (alpha0 < 0 && beta0 >= 0)
    if par.n0(i+1) > par.p0(i-1)
        rec_zone_loc(i) = 'R';
    elseif par.n0(i+1) < par.p0(i-1)
        rec_zone_loc(i) = 'L';
    else
        rec_zone_loc(i) = 'C';
    end
elseif (alpha0 >= 0 && beta0 < 0) || (alpha0 > 0 && beta0 <= 0)
    if par.p0(i+1) > par.n0(i-1)
        rec_zone_loc(i) = 'R';
    elseif par.p0(i+1) < par.n0(i-1)
        rec_zone_loc(i) = 'L';
    else
        rec_zone_loc(i) = 'C';
    end
elseif (alpha0 <= 0 && beta0 < 0) || (alpha0 < 0 && beta0 <= 0)
    rec_zone_loc(i) = 'L';
elseif (alpha0 >= 0  && beta0 > 0) || (alpha0 > 0  && beta0 >= 0)
    rec_zone_loc(i) = 'R';
else
    rec_zone_loc(i) = 'C';
end

end