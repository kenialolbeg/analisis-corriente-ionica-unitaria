function MSE_analysis()
    % Cargar los datos del archivo
    data = load('k_channel.txt');
    time = data(:,1);
    current = data(:,2);
    
    % Parámetros para el cálculo de MSE
    m = 5;      % Dimensión de embedding (típicamente 2)
    r = 0.15;   % Radio como fracción de la desviación estándar (típicamente 0.15-0.25)
    tau = 54;    % Factor de escala inicial
    max_scale = 20; % Máximo factor de escala a calcular
    
    % Calcular MSE
    [mse_values, scales] = multiscale_sample_entropy(current, m, r, tau, max_scale);
    
    % Graficar resultados
    figure;
    plot(scales, mse_values, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    xlabel('Factor de Escala (τ)');
    ylabel('Entropía Muestral (SE)');
    title('Entropía Multiescala (MSE) de la Corriente');
    grid on;
    
    % Mostrar valores en consola
    disp('Resultados de Entropía Multiescala:');
    disp(table(scales', mse_values', 'VariableNames', {'Escala', 'Entropia'}));
end

function [mse, scales] = multiscale_sample_entropy(signal, m, r, tau, max_scale)
    % Calcular la desviación estándar de la señal original
    sd = std(signal);
    
    % Inicializar vectores de resultados
    scales = 1:max_scale;
    mse = zeros(size(scales));
    
    % Calcular Sample Entropy para cada escala
    for i = 1:length(scales)
        scale = scales(i);
        
        % Coarse-graining de la señal
        if scale == 1
            y = signal;
        else
            y = coarse_grain(signal, scale);
        end
        
        % Calcular Sample Entropy para la señal coarse-grained
        mse(i) = sample_entropy(y, m, r*sd);
    end
end

function y = coarse_grain(signal, scale)
    % Aplicar coarse-graining a la señal
    n = length(signal);
    y_length = floor(n/scale);
    y = zeros(1, y_length);
    
    for i = 1:y_length
        y(i) = mean(signal((i-1)*scale + 1 : i*scale));
    end
end

function se = sample_entropy(signal, m, r)
    % Implementación de Sample Entropy
    n = length(signal);
    
    % Prevenir división por cero
    if n <= m
        se = 0;
        return;
    end
    
    % Calcular patrones similares
    phi = zeros(1, 2);
    
    for k = m:m+1
        count = 0;
        patterns = zeros(n - k + 1, k);
        
        % Crear patrones de dimensión k
        for i = 1:n - k + 1
            patterns(i,:) = signal(i:i + k - 1);
        end
        
        % Calcular distancia entre patrones
        for i = 1:n - k
            for j = i+1:n - k + 1
                if max(abs(patterns(i,:) - patterns(j,:))) <= r
                    count = count + 1;
                end
            end
        end
        
        % Calcular phi
        phi(k - m + 1) = count / ((n - k + 1)*(n - k)/2);
    end
    
    % Calcular Sample Entropy
    if phi(2) == 0 || phi(1) == 0
        se = -log(1/((n - m + 1)*(n - m)/2));
    else
        se = -log(phi(2)/phi(1));
    end
end