clear; clc;

% --- Datos del problema ---
V = [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4];
P = [300.2, 242.1, 201.5, 169.8, 151, 134.3, 120.8, 107.9, 99.1, 93.4, 86.2, 79.5, 75.1];
h = 0.25;

% --- VALOR REAL DE REFERENCIA ---
valor_real = 417.8510; 
% 1. CÁLCULO DE MÉTODOS
% --- Trapecio Compuesto ---
W_trap = (h / 2) * (P(1) + 2 * sum(P(2:end-1)) + P(end));

% --- Simpson 1/3 Compuesto ---
idx_cuatro = 2:2:12;
idx_dos_s13 = 3:2:11;
W_simps13 = (h / 3) * (P(1) + 4 * sum(P(idx_cuatro)) + 2 * sum(P(idx_dos_s13)) + P(end));

% --- Simpson 3/8 Compuesto ---
idx_dos_s38 = 4:3:10;
idx_todos_internos = 2:12;
idx_tres = setdiff(idx_todos_internos, idx_dos_s38);
W_simps38 = (3 * h / 8) * (P(1) + 3 * sum(P(idx_tres)) + 2 * sum(P(idx_dos_s38)) + P(end));

% --- Gauss-Legendre (m=3) con Spline ---
pp = spline(V, P); 
t = [-sqrt(3/5), 0, sqrt(3/5)];
w = [5/9, 8/9, 5/9];
a = 1; b = 4;
V_gauss = ((b - a) / 2) * t + ((b + a) / 2);
P_gauss = ppval(pp, V_gauss);
W_gauss = ((b - a) / 2) * sum(w .* P_gauss);
% 2. CÁLCULO DE ERRORES RELATIVOS PORCENTUALES
err_trap   = abs((valor_real - W_trap) / valor_real) * 100;
err_s13    = abs((valor_real - W_simps13) / valor_real) * 100;
err_s38    = abs((valor_real - W_simps38) / valor_real) * 100;
err_gauss  = abs((valor_real - W_gauss) / valor_real) * 100;

fprintf('\n');
fprintf('MÉTODO          | VALOR APROXIMADO | ERROR RELATIVO %%\n');
fprintf('-------------------------------------------------------\n');
fprintf('Trapecio        |   %14.4f |   %14.4f\n', W_trap, err_trap);
fprintf('Simpson 1/3     |   %14.4f |   %14.4f\n', W_simps13, err_s13);
fprintf('Simpson 3/8     |   %14.4f |   %14.4f\n', W_simps38, err_s38);
fprintf('Gauss-Legendre  |   %14.4f |   %14.4f\n', W_gauss, err_gauss);
fprintf('-------------------------------------------------------\n');