%% ACTIVIDAD PREVIA - Laboratorio 2: Interferencia Intersimbolo (ISI)
%  Grafica He(f) y he(t) del filtro coseno alzado
%  para alpha = 0, 0.25, 0.75, 1

clear; close all; clc;

%% Parametros
f0     = 1000;                  % Frecuencia de corte (Hz)
alphas = [0, 0.25, 0.75, 1];   % Factores de roll-off
Ts     = 1/(2*f0);             % Periodo de simbolo sin ISI

% Ejes de frecuencia y tiempo
Bmax = f0*(1+1);                        % B cuando alpha = 1
f    = linspace(-2*Bmax, 2*Bmax, 4001); % eje de frecuencia
tmax = 8*Ts;
t    = linspace(0, tmax, 4001);         % eje de tiempo (t >= 0)

%% ----------------------------------------------------------
%  RESPUESTA EN FRECUENCIA  He(f)
%  ----------------------------------------------------------
figure; set(gcf,'Color','w');
sgtitle('Respuesta en frecuencia del coseno alzado', 'FontSize', 13);

for k = 1:length(alphas)
    a  = alphas(k);
    fD = a*f0;       % f_delta = a*f0
    f1 = f0 - fD;   % borde interior de la transicion
    B  = f0 + fD;   % ancho de banda total

    He = zeros(1, length(f));

    for i = 1:length(f)
        af = abs(f(i));
        if af < f1
            He(i) = 1;                                       % zona plana
        elseif (a > 0) && (af >= f1) && (af < B)
            He(i) = 0.5*(1 + cos(pi*(af - f1)/(2*fD)));    % transicion
        else
            He(i) = 0;                                       % fuera de banda
        end
    end

    subplot(2,2,k);
    plot(f, He, 'LineWidth', 1.3); grid on;
    xlabel('f (Hz)'); ylabel('|He(f)|');
    title(['\alpha = ' num2str(a)]);
    xlim([-2*B 2*B]); ylim([0 1.1]);
end

%% ----------------------------------------------------------
%  RESPUESTA AL IMPULSO  he(t)
%  ----------------------------------------------------------
figure; set(gcf,'Color','w');
sgtitle('Respuesta al impulso del coseno alzado', 'FontSize', 13);

for k = 1:length(alphas)
    a  = alphas(k);
    fD = a*f0;
    he = zeros(1, length(t));

    for i = 1:length(t)
        x = t(i);

        % sinc normalizado: sin(2*pi*f0*x) / (2*pi*f0*x)
        if abs(x) < 1e-12
            sinc_val = 1;          % limite en x = 0
        else
            sinc_val = sin(2*pi*f0*x) / (2*pi*f0*x);
        end

        if a == 0
            he(i) = 2*f0 * sinc_val;   % caso alpha = 0
        else
            % factor cos/den con manejo de singularidades
            den = 1 - (4*fD*x)^2;
            if abs(abs(4*fD*x) - 1) < 1e-6
                frac = pi/4;           % limite cuando 4*fD*t = +-1
            else
                frac = cos(2*pi*fD*x) / den;
            end
            he(i) = 2*f0 * sinc_val * frac;
        end
    end

    subplot(2,2,k);
    plot(t, he, 'LineWidth', 1.3); grid on; hold on;

    % Marcas en k*Ts para verificar ceros de Nyquist
    ks = 0:floor(tmax/Ts);
    stem(ks*Ts, zeros(size(ks)), 'filled', ...
         'MarkerSize', 3, 'LineStyle', 'none', 'Color', 'r');
    hold off;

    xlabel('t (s)'); ylabel('h_e(t)');
    title(['\alpha = ' num2str(a)]);
    xlim([0 tmax]);
end

%% Info en consola
fprintf('\nTs = 1/(2*f0) = %.6f s\n', Ts);
for m = 1:4
    fprintf('Cero en t = %d*Ts = %.6f s\n', m, m*Ts);
end
