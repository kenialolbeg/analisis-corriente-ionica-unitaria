function [FNN] = knn_deneme(x,tao,mmax,rtol,atol)
%x : time series
%tao : retraso temporal
%mmax : dimensión de incrustación máxima
%rtol: La tolerancia relativa para determinar si un vecino cercano es falso
%atol: La tolerancia absoluta para determinar si un vecino cercano es falso
%reference:M. B. Kennel, R. Brown, and H. D. I. Abarbanel, Determining
%embedding dimension for phase-space reconstruction using a geometrical 
%construction, Phys. Rev. A 45, 3403 (1992). 
%author:"Merve Kizilkaya"
%rtol=15
%atol=2;
N=length(x); % Calculo tamaño serie de tiempo
Ra=std(x,1);  %Caclculo de la desviacion estandar de la serie de tiempo

for m=1:mmax  %Bucle para cada dimension de incrustación
    M=N-m*tao;
    Y=psr_deneme(x,m,tao,M); %Cálculo de la matriz de puntos en el espacio fase. Esto crea una matriz Y ue representa los puntos en el espacio fase para la dimension de incrustacion actual
    FNN(m,1)=0;
    for n=1:M %Bucle para cada punto en el espacio de fases
        y0=ones(M,1)*Y(n,:);
        distance=sqrt(sum((Y-y0).^2,2)); %Calculo de la distancia euclidiana a los vecinos mas cercanos
        [neardis nearpos]=sort(distance);
        
        D=abs(x(n+m*tao)-x(nearpos(2)+m*tao));
        R=sqrt(D.^2+neardis(2).^2);
        if D/neardis(2) > rtol || R/Ra > atol  %Comprueba si el punto actual tiene vecinos que son consideradso "falsos" segun los criterios de tolerancia establecidos
             FNN(m,1)=FNN(m,1)+1;
        end
    end
end
FNN=(FNN./FNN(1,1))*100; %Calcula el porcentaje de vecinos falsos mas cercanos para cada dimension de incrustacion
figure %visualizacion del resultado
figure
plot(1:length(FNN),FNN)
grid on;
title('Mínima dimensión embebida con falsos vecinos cercanos')
xlabel('Dimensión Embebida')
ylabel('Porcentaje de falsos vecinos cercanos')

function Y=psr_deneme(x,m,tao,npoint)
%Phase space reconstruction
%x : time series 
%m : embedding dimension
%tao : time delay
%npoint : total number of reconstructed vectors
%Y : M x m matrix
% author:"Merve Kizilkaya"
N=length(x);
if nargin == 4
    M=npoint;
else
    M=N-(m-1)*tao;
end

Y=zeros(M,m); 

for i=1:m
    Y(:,i)=x((1:M)+(i-1)*tao)';
end
