clear; clc; close all;

% 1. SELECCIÓN DEL ARCHIVO
fprintf('Seleccione el archivo de Excel con los 49 puntos de la mano...\n');
[archivo, ruta] = uigetfile({'*.xlsx;*.xls', 'Archivos de Excel (*.xlsx, *.xls)'}, 'Seleccione el archivo de datos');

if isequal(archivo, 0)
    disp('Usuario canceló la selección.');
    return;
else
    fullPath = fullfile(ruta, archivo);
    datos = readmatrix(fullPath);
    x = datos(:, 1)'; 
    y = datos(:, 2)'; 
end

% 2. CONFIGURACIÓN
n = length(x);
T = 1:n;               % Parámetro de tiempo (puntos de control)
tt = linspace(1, n, 1000); % 1000 puntos para una suavidad total

% 3. MÉTODO DE SPLINES CÚBICOS (La solución definitiva)
% Calculamos el spline para X y para Y por separado respecto al tiempo T
x_spline = spline(T, x, tt);
y_spline = spline(T, y, tt);

% --- REPORTE EN CONSOLA ---
fprintf('\n=========================================');
fprintf('\n       DIBUJO COMPLETO DE LA MANO        ');
fprintf('\n=========================================');
fprintf('\n Método: Splines Cúbicos Paramétricos');
fprintf('\n Puntos procesados: %d', n);
fprintf('\n Estado: Curva suavizada y completa');
fprintf('\n=========================================\n');

% 4. GRÁFICA
figure('Color', 'w', 'Name', 'Mano Completa - Spline Cúbico');
hold on; 

% Dibujar la curva del Spline (Línea continua y suave)
plot(x_spline, y_spline, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Trazo de la Mano');

% Dibujar los puntos originales del Excel
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Puntos Originales');

% Configuración de vista
grid on;
axis equal; % CRUCIAL: Esto evita que la mano se vea gorda o flaca
title(['Dibujo de Mano Completa - Splines desde: ', archivo]);
xlabel('Coordenada X');
ylabel('Coordenada Y');
legend('Location', 'best');

% Añadir etiquetas numeradas a los puntos (opcional, para verificar orden)
% for i = 1:n
%     text(x(i), y(i), num2str(i), 'FontSize', 8, 'VerticalAlignment', 'bottom');
% end

hold off;
fprintf('\n¡Proceso terminado! Ahora deberías ver la mano completa y fluida.\n');