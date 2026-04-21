% =========================================================================
% INTERPOLACIÓN PARAMÉTRICA (Mano - 49 Puntos)
% Métodos: Lagrange, Newton y Matricial (Sin gráficas)
% =========================================================================
clear; clc;

% -------------------------------------------------------------------------
% 1. CARGA DE DATOS (Selección de Excel)
% -------------------------------------------------------------------------
fprintf('Seleccione el archivo de Excel con los datos...\n');
[archivo, ruta] = uigetfile({'*.xlsx;*.xls', 'Archivos de Excel'}, 'Seleccione el Excel');

if isequal(archivo, 0)
    disp('Operación cancelada por el usuario.');
    return;
end

datos = readmatrix(fullfile(ruta, archivo));
x = datos(:, 1)'; % Coordenadas X
y = datos(:, 2)'; % Coordenadas Y
n = length(x);

% Parámetros de interpolación
T = 1:n; % Vector de tiempo paramétrico (1 a 49)
tt = linspace(1, n, 500); % Puntos a evaluar

fprintf('\nArchivo cargado exitosamente: %s\n', archivo);
fprintf('Procesando %d puntos...\n', n);

% -------------------------------------------------------------------------
% 2. MÉTODO 1: POLINOMIO DE LAGRANGE
% -------------------------------------------------------------------------
disp('Calculando Lagrange...');
x_lagrange = zeros(1, length(tt));
y_lagrange = zeros(1, length(tt));

for k = 1:length(tt)
    sum_x = 0; sum_y = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if i ~= j
                L = L * (tt(k) - T(j)) / (T(i) - T(j));
            end
        end
        sum_x = sum_x + L * x(i);
        sum_y = sum_y + L * y(i);
    end
    x_lagrange(k) = sum_x;
    y_lagrange(k) = sum_y;
end

% -------------------------------------------------------------------------
% 3. MÉTODO 2: POLINOMIO DE NEWTON (Diferencias Divididas)
% -------------------------------------------------------------------------
disp('Calculando Newton...');
Dx = zeros(n, n); Dx(:,1) = x'; 
Dy = zeros(n, n); Dy(:,1) = y'; 

% Llenado de matrices de diferencias divididas
for j = 2:n
    for i = 1:n-j+1
        Dx(i,j) = (Dx(i+1,j-1) - Dx(i,j-1)) / (T(i+j-1) - T(i));
        Dy(i,j) = (Dy(i+1,j-1) - Dy(i,j-1)) / (T(i+j-1) - T(i));
    end
end

coef_newton_x = Dx(1,:);
coef_newton_y = Dy(1,:);

x_newton = zeros(1, length(tt));
y_newton = zeros(1, length(tt));

% Evaluación del polinomio
for k = 1:length(tt)
    x_newton(k) = coef_newton_x(1);
    y_newton(k) = coef_newton_y(1);
    prod_term = 1;
    for i = 2:n
        prod_term = prod_term * (tt(k) - T(i-1));
        x_newton(k) = x_newton(k) + coef_newton_x(i) * prod_term;
        y_newton(k) = y_newton(k) + coef_newton_y(i) * prod_term;
    end
end

% -------------------------------------------------------------------------
% 4. MÉTODO 3: POLINOMIO MATRICIAL (Matriz de Vandermonde)
% -------------------------------------------------------------------------
disp('Calculando Matricial (Vandermonde)...');

% Construcción de la Matriz de Vandermonde manual
V = zeros(n, n);
for i = 1:n
    for j = 1:n
        V(i, j) = T(i)^(j-1);
    end
end

% Desactivar advertencia de matriz mal condicionada temporalmente (explicación abajo)
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:illConditionedMatrix');

% Resolución de los sistemas de ecuaciones lineales (V * c = x)
coef_mat_x = V \ x';
coef_mat_y = V \ y';

warning('on', 'all'); % Reactivar advertencias

x_matricial = zeros(1, length(tt));
y_matricial = zeros(1, length(tt));

% Evaluación del polinomio matricial
for k = 1:length(tt)
    for j = 1:n
        x_matricial(k) = x_matricial(k) + coef_mat_x(j) * (tt(k)^(j-1));
        y_matricial(k) = y_matricial(k) + coef_mat_y(j) * (tt(k)^(j-1));
    end
end

% -------------------------------------------------------------------------
% 5. REPORTE FINAL
% -------------------------------------------------------------------------
fprintf('\n=========================================');
fprintf('\n       ALGORITMOS COMPLETADOS            ');
fprintf('\n=========================================');
fprintf('\n Matrices resultantes guardadas en memoria:');
fprintf('\n - x_lagrange,  y_lagrange');
fprintf('\n - x_newton,    y_newton');
fprintf('\n - x_matricial, y_matricial');
fprintf('\n=========================================\n');