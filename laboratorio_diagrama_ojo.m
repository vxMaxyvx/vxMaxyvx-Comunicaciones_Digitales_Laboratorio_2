%% LABORATORIO 2: Interferencia Intersimbolo (ISI)
%  Simulacion de transmision NRZ-L con filtro coseno alzado y diagrama de ojo

clear; clc; close all;

%% ----------------------------------------------------------
%  Parametros
%  ----------------------------------------------------------
f0              = 10000;   % Frecuencia de corte del filtro (Hz)
factorTiempoMax = 8;       % Extension temporal del filtro (multiplos de Ts)
muestrasTiempo  = 4001;    % Numero de muestras del filtro
alpha           = 0.75;    % Factor de roll-off  (probar tambien 0.25)
SNR_dB          = 25;      % Relacion senal-ruido en dB
N_simbolos      = 10000;   % Numero de simbolos aleatorios
Ts              = 1/(2*f0);% Periodo de simbolo

%% ----------------------------------------------------------
%  Generar y normalizar el filtro coseno alzado
%  ----------------------------------------------------------
he = generar_he(alpha, f0, factorTiempoMax, muestrasTiempo);
he = he / max(abs(he));

%% ----------------------------------------------------------
%  Secuencia NRZ-L aleatoria  (-1 o +1)
%  ----------------------------------------------------------
simbolos = 2*randi([0 1], 1, N_simbolos) - 1;

% Sobremuestreo: insertar ceros entre simbolos
samples_per_symbol            = round(muestrasTiempo / (2*factorTiempoMax));
simbolos_sobremuestreados     = zeros(1, N_simbolos * samples_per_symbol);
simbolos_sobremuestreados(1:samples_per_symbol:end) = simbolos;

%% ----------------------------------------------------------
%  Filtrado con coseno alzado + ruido gaussiano blanco
%  ----------------------------------------------------------
senal_filtrada = conv(simbolos_sobremuestreados, he, 'same');
senal_filtrada = awgn(senal_filtrada, SNR_dB, 'measured');

%% ----------------------------------------------------------
%  Construccion del diagrama de ojo
%  ----------------------------------------------------------
muestras_por_ventana = 2 * samples_per_symbol;

% Descartar transitorios del inicio y del final
offset     = round(factorTiempoMax * samples_per_symbol);
senal_util = senal_filtrada(offset:end-offset);

num_trazos = floor((length(senal_util) - muestras_por_ventana) ...
                   / samples_per_symbol);

eye_data = zeros(num_trazos, muestras_por_ventana);
for i = 1:num_trazos
    inicio = (i-1)*samples_per_symbol + 1;
    fin    = inicio + muestras_por_ventana - 1;
    if fin <= length(senal_util)
        eye_data(i, :) = senal_util(inicio:fin);
    end
end

%% ----------------------------------------------------------
%  Grafico: Diagrama de ojo
%  ----------------------------------------------------------
figure('Position', [100 100 800 600], 'Color', 'w');
t_ojo = linspace(0, 2, muestras_por_ventana);   % eje en simbolos

hold on;
for i = 1:num_trazos
    plot(t_ojo, eye_data(i, :), 'b', 'LineWidth', 0.5);
end
hold off;

ax = gca;
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', 'k', 'GridAlpha', 0.2, 'FontSize', 12);

xlabel('Tiempo [simbolos]', 'Color', 'k');
ylabel('Amplitud',          'Color', 'k');
title(sprintf('Diagrama de Ojo Manual  (\\alpha = %.2f,  f_0 = %d Hz,  SNR = %.1f dB)', ...
              alpha, f0, SNR_dB), 'Color', 'k');
grid on;
xlim([0 2]);

% Marcar punto de muestreo optimo
line([1 1], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
text(1.05, max(ylim)*0.9, 'Punto de muestreo', 'Color', 'r', 'FontSize', 10);

%% ----------------------------------------------------------
%  Grafico: Respuesta al impulso del filtro
%  ----------------------------------------------------------
figure('Position', [100 100 800 400]);
tmax_plot = factorTiempoMax * Ts;
t_filtro  = linspace(-tmax_plot, tmax_plot, muestrasTiempo);

plot(t_filtro/Ts, he, 'LineWidth', 2);
xlabel('Tiempo [simbolos]'); ylabel('Amplitud');
title(['Respuesta al Impulso del Filtro Raised Cosine  (\alpha = ' num2str(alpha) ')']);
grid on;


%% ==========================================================
%%  FUNCION LOCAL: generar filtro coseno alzado he(t)
%% ==========================================================
function salida = generar_he(alpha, f0, factorTiempoMax, muestrasTiempo)
%GENERAR_HE  Calcula la respuesta al impulso del filtro coseno alzado.
%
%   Entradas:
%     alpha           - Factor de roll-off [0, 1]
%     f0              - Frecuencia de corte (Hz)
%     factorTiempoMax - Extension temporal en multiplos de Ts
%     muestrasTiempo  - Numero de muestras
%
%   Salida:
%     salida          - Vector he(t)

    Ts   = 1/(2*f0);
    tmax = factorTiempoMax * Ts;
    t    = linspace(-tmax, tmax, muestrasTiempo);
    fD   = alpha * f0;
    he   = zeros(1, length(t));

    for i = 1:length(t)
        x = t(i);

        if abs(x) < 1e-10
            sinc_val = 1;
        else
            sinc_val = sin(2*pi*f0*x) / (2*pi*f0*x);
        end

        den = 1 - (4*fD*x)^2;

        if abs(den) < 1e-10
            % Casos donde el denominador es cero
            if abs(x) < 1e-10
                he(i) = 2*f0;
            else
                he(i) = (pi*2*f0/4) * sin(pi/(2*alpha));
            end
        else
            frac  = cos(2*pi*fD*x) / den;
            he(i) = 2*f0 * sinc_val * frac;
        end
    end

    salida = he;
end
