function optimal_tau = calculate_optimal_tau(filename)
    % Cargar los datos del archivo
    data = load("k_channel.txt");
    time = data(:,1);
    current = data(:,2);
    
    % Normalizar la corriente
    current = (current - mean(current)) / std(current);
    
    % Parámetros para el cálculo de la información mutua
    max_tau = 100; % Máximo retardo a considerar (en muestras)
    num_bins = 20; % Número de bins para el histograma
    
    % Calcular información mutua para diferentes retardos
    mutual_info = zeros(max_tau, 1);
    
    for tau = 1:max_tau
        % Crear vectores desplazados
        x = current(1:end-tau);
        y = current(1+tau:end);
        
        % Calcular histograma 2D
        [counts, ~] = hist3([x y], [num_bins num_bins]);
        joint_prob = counts / sum(counts(:));
        
        % Calcular probabilidades marginales
        prob_x = sum(joint_prob, 2);
        prob_y = sum(joint_prob, 1);
        
        % Calcular información mutua
        mi = 0;
        for i = 1:num_bins
            for j = 1:num_bins
                if joint_prob(i,j) > 0 && prob_x(i) > 0 && prob_y(j) > 0
                    mi = mi + joint_prob(i,j) * log2(joint_prob(i,j) / (prob_x(i) * prob_y(j)));
                end
            end
        end
        
        mutual_info(tau) = mi;
    end
    
    % Encontrar el primer mínimo local de la información mutua
    [~, optimal_tau] = findpeaks(-mutual_info, 'NPeaks', 1);
    
    % Si no se encuentra un mínimo claro, usar un valor por defecto
    if isempty(optimal_tau)
        optimal_tau = 1;
    end
    
    % Graficar la información mutua
    figure;
    plot(1:max_tau, mutual_info, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(optimal_tau, mutual_info(optimal_tau), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Retardo \tau (muestras)');
    ylabel('Información Mutua');
    title('Información Mutua vs Retardo Temporal');
    legend('Información Mutua', 'Retardo Óptimo');
    grid on;
    
    fprintf('El retardo temporal óptimo (tau) es: %d muestras\n', optimal_tau);
end