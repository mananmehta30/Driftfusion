
par_memristor = pc('Input_files/memristor');

soleq_memristor = equilibrate(par_memristor);


sol_CV = doCV(soleq_memristor.ion, 0, 0, 0, 0, 0, 0, 0);


