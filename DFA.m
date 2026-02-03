% Cargar datos
data = load('k_channel_completa.txt');  % Asegúrate de que el archivo esté en el mismo directorio
t = data(:,1);                 % Tiempo
I = data(:,2);                 % Corriente unitaria

% Paso 1: Integrar la serie (restando la media)
I_mean = mean(I);
Y = cumsum(I - I_mean);

% Paso 2: Definir tamaños de ventana
N = length(Y);
n_vals = round(logspace(log10(10), log10(N/4), 20));  % 20 escalas logarítmicas

F_n = zeros(size(n_vals));  % Fluctuaciones

for i = 1:length(n_vals)
    n = n_vals(i);
    num_segments = floor(N / n);
    fluct = zeros(num_segments, 1);
    
    for j = 1:num_segments
        idx_start = (j-1)*n + 1;
        idx_end = j*n;
        segment = Y(idx_start:idx_end);
        
        % Ajuste lineal (polinomio de grado 1)
        x = (1:n)';
        p = polyfit(x, segment, 1);
        trend = polyval(p, x);
        
        % Fluctuación cuadrática
        fluct(j) = mean((segment - trend).^2);
    end
    
    % Fluctuación promedio para tamaño n
    F_n(i) = sqrt(mean(fluct));
end

% Paso 3: Ajuste log-log para obtener alfa
log_n = log10(n_vals);
log_F = log10(F_n);
coeffs = polyfit(log_n, log_F, 1);
alpha = coeffs(1);

% Mostrar resultado
fprintf('Estimación del exponente de escalamiento α: %.4f\n', alpha);

% Graficar
figure;
loglog(n_vals, F_n, 'o');
hold on;
loglog(n_vals, 10.^(polyval(coeffs, log_n)), '--r');
xlabel('Tamaño de ventana n');
ylabel('Fluctuación F(n)');
title(['Estimación de α = ', num2str(alpha)]);
legend('F(n)', 'Ajuste lineal');
grid on;
