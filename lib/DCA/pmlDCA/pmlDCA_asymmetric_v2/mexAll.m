fprintf('Compiling g_rC.c\n');
mex -outdir functions functions/g_rC.c

fprintf('Compiling calc_inverse_weights.c\n');
mex -outdir functions functions/calc_inverse_weights.c
