%-----------------------------------------------------------------
%%%% Detector de Epilepsia en Dominio Frecuencial
%
%%%% Roberto Rueda Galindo
%
%%%% 28 de Noviembre 2016
%
%-----------------------------------------------------------------

clear all
close all
clc

%%Señal
[hdr,eeg]=edfread('file_03.edf'); %Carga la señal 
eeg=eeg.';
Fs=hdr.samples(1); %Frecuencia de muestreo

%Visualización
figure;
eeg2=eeg+repmat([0:hdr.ns-1]*600,hdr.samples(1)*hdr.records,1); % Separación para la visualización
plot([0:size(eeg,1)-1]/Fs,eeg2)
set(gca,'YTick',[0:hdr.ns-1]*600,'YTickLabel',hdr.label,'FontSize',10)
grid on %Activar Cuadricula
grid minor
hold on %Mantener Figura activa
%onset=1467; %Inicio epilepsia
%endtime=1494; %Fin epilepsia
%h1=line([onset;onset],[min(eeg(:)) max(eeg(:))+600*hdr.ns]); %Marca de Inicio
%set(h1,'linewidth',1,'color','k','linestyle','-.'); %Caracteristicas Marca Inicio
%h2=line([endtime;endtime],[min(eeg(:)) max(eeg(:))+600*hdr.ns]);%Marca de Fin
%set(h2,'linewidth',1,'color','k','linestyle','--'); %Caracteristicas Marca Fin
ylabel('Canales');
title('Señal EEG')
axis([0 3600 -1000 17000]) %Ajustar Gráfica


%%Filtrado
%Filtro eliminación de ruido Paso Banda (0.5 a 60Hz)
[b1,a1] = butter(2,[0.5,60]/(Fs/2)); %Filtro Paso Banda
eegpb = filtfilt(b1,a1,eeg); %Señal después del filtrado Paso Banda

%Filtro eliminación de ruido Notch (50Hz)
Wo = 50/(Fs/2); Bw = Wo/35;
[b2,a2] = iirnotch(Wo,Bw); %Filtro Notch
eegnotch = filtfilt(b2,a2,eegpb); %Señal después del filtrado Notch

%%Filtrado de ondas cerebrales
%Ondas Delta (1-4Hz)
[b_delta,a_delta] = butter(2,[1,4]/(Fs/2)); 
Onda_Delta = filtfilt(b_delta,a_delta,eeg); 
%Ondas Theta (4-8Hz)
[b_theta,a_theta] = butter(2,[4,8]/(Fs/2)); 
Onda_Theta = filtfilt(b_theta,a_theta,eeg); 
%Ondas Alpha (8-13Hz)
[b_alpha,a_alpha] = butter(2,[8,13]/(Fs/2)); 
Onda_Alpha = filtfilt(b_alpha,a_alpha,eeg); 
%Ondas Beta (13-30Hz)
[b_beta,a_beta] = butter(2,[13,30]/(Fs/2)); 
Onda_Beta = filtfilt(b_beta,a_beta,eeg); 
%Ondas Gamma (30-100Hz)
[b_gamma,a_gamma] = butter(2,[30,100]/(Fs/2)); 
Onda_Gamma = filtfilt(b_gamma,a_gamma,eeg); 


%% Estudio de la Energía de cada EEG
%Energía
%Ondas Delta (1-4Hz)
epoch_dur=512; % Número de muestras para cada EEG
nepoch=size(Onda_Delta,1)/epoch_dur;
ene_delta=reshape(Onda_Delta,[epoch_dur,nepoch,hdr.ns]);
ene_delta=sum(ene_delta.^2)/nepoch;
ene_delta=squeeze(ene_delta); %Energía
%Energía
%Ondas Theta (4-8Hz)
epoch_dur=512; % Número de muestras para cada EEG
nepoch=size(Onda_Theta,1)/epoch_dur;
ene_theta=reshape(Onda_Theta,[epoch_dur,nepoch,hdr.ns]);
ene_theta=sum(ene_theta.^2)/nepoch;
ene_theta=squeeze(ene_theta); %Energía
%Energía
%Ondas Alpha (8-13Hz)
epoch_dur=512; % Número de muestras para cada EEG
nepoch=size(Onda_Alpha,1)/epoch_dur;
ene_alpha=reshape(Onda_Alpha,[epoch_dur,nepoch,hdr.ns]);
ene_alpha=sum(ene_alpha.^2)/nepoch;
ene_alpha=squeeze(ene_alpha); %Energía
%Energía
%Ondas Beta (13-30Hz)
epoch_dur=512; % Número de muestras para cada EEG
nepoch=size(Onda_Beta,1)/epoch_dur;
ene_beta=reshape(Onda_Beta,[epoch_dur,nepoch,hdr.ns]);
ene_beta=sum(ene_beta.^2)/nepoch;
ene_beta=squeeze(ene_beta); %Energía
%Energía
%Ondas Gamma (30-100Hz)
epoch_dur=512; % Número de muestras para cada EEG
nepoch=size(Onda_Gamma,1)/epoch_dur;
ene_gamma=reshape(Onda_Gamma,[epoch_dur,nepoch,hdr.ns]);
ene_gamma=sum(ene_gamma.^2)/nepoch;
ene_gamma=squeeze(ene_gamma); %Energía


%%Filtro para eliminar falsas alarmas
%Filtro de Mediana Ondas Delta
ums=7;%Ventana Umbral
enem_Delta=medfilt1 (ene_delta,ums);
%Filtro de Mediana Ondas Theta
enem_Theta=medfilt1 (ene_theta,ums);
%Filtro de Mediana Ondas Alpha
enem_Alpha=medfilt1 (ene_alpha,ums);
%Filtro de Mediana Ondas Beta
enem_Beta=medfilt1 (ene_beta,ums);
%Filtro de Mediana Ondas Theta
enem_Gamma=medfilt1 (ene_gamma,ums);

%%Detección de epilepsia por umbral
%Umbral de detección para Ondas Delta
umbral = ceil(prctile(enem_Delta(:), 95));
Det_Delta = zeros(size(enem_Delta));
Det_Delta(enem_Delta > umbral) = 1;
%Visualización
figure;
subplot(2,3,1)
Det_Delta2 = Det_Delta  + repmat([0:hdr.ns-1]*2,nepoch,1); %Separación para la visualización
plot([0:nepoch-1]*epoch_dur/Fs,Det_Delta2);
grid on %Activar Cuadricula
grid minor
hold on %Mantener Figura activa
%ylabel('Canales');
title('Detección Epilepsia en Ondas Delta')
axis([0 3600 -5 55])

%Umbral de detección para Ondas Theta
umbral = ceil(prctile(enem_Theta(:), 95));
Det_Theta = zeros(size(enem_Theta));
Det_Theta(enem_Theta > umbral) = 1;
%Visualización
subplot(2,3,2)
Det_Theta2 = Det_Theta  + repmat([0:hdr.ns-1]*2,nepoch,1); %Separación para la visualización
plot([0:nepoch-1]*epoch_dur/Fs,Det_Theta2);
grid on %Activar Cuadricula
grid minor
hold on %Mantener Figura activa
%ylabel('Canales');
title('Detección Epilepsia en Ondas Theta')
axis([0 3600 -5 55])


%Umbral de detección para Ondas Alpha
umbral = ceil(prctile(enem_Alpha(:), 95));
Det_Alpha = zeros(size(enem_Alpha));
Det_Alpha(enem_Alpha > umbral) = 1;
%Visualización
%figure
subplot(2,3,3)
Det_Alpha2 = Det_Alpha  + repmat([0:hdr.ns-1]*2,nepoch,1); %Separación para la visualización
plot([0:nepoch-1]*epoch_dur/Fs,Det_Alpha2);
grid on %Activar Cuadricula
grid minor
hold on %Mantener Figura activa
%ylabel('Canales');
title('Detección Epilepsia en Ondas Alpha')
axis([0 3600 -5 55])

%Umbral de detección para Ondas Beta
umbral = ceil(prctile(enem_Beta(:), 95));
Det_Beta = zeros(size(enem_Beta));
Det_Beta(enem_Beta > umbral) = 1;
%Visualización
%figure
subplot(2,3,4)
Det_Beta2 = Det_Beta  + repmat([0:hdr.ns-1]*2,nepoch,1); %Separación para la visualización
plot([0:nepoch-1]*epoch_dur/Fs,Det_Beta2);
grid on %Activar Cuadricula
grid minor
hold on %Mantener Figura activa
%ylabel('Canales');
title('Detección Epilepsia en Ondas Beta')
axis([0 3600 -5 55])

%Umbral de detección para Ondas Gamma
umbral = ceil(prctile(enem_Gamma(:), 95));
Det_Gamma = zeros(size(enem_Gamma));
Det_Gamma(enem_Gamma > umbral) = 1;
%Visualización
%figure
subplot(2,3,5)
Det_Gamma2 = Det_Gamma  + repmat([0:hdr.ns-1]*2,nepoch,1); %Separación para la visualización
plot([0:nepoch-1]*epoch_dur/Fs,Det_Gamma2);
grid on %Activar Cuadricula
grid minor
hold on %Mantener Figura activa
%ylabel('Canales');
title('Detección Epilepsia en Ondas Gamma')
axis([0 3600 -5 55])



epi_r=Det_Beta+Det_Theta;
epi = zeros(size(Det_Beta));
epi(epi_r>=1)=1;

%Buscar Canales Activos
epi_sum = sum(epi,2);%Sumo todos los canales para ver los que estan activos
epi_d = zeros(size(epi(:,1)));%Hago una matriz de ceros

%canales=15; %Canales abiertos para considerar epilepsia
%epi_d (epi_sum > canales)=1;%Comparo para ver si el numero de canales abiertos supera el umbral si es asi 1 si no 0
%epi_t=epi_d'; %transpongo la matriz

%Ajuste de canales activos para encontrar la epilepsia
for i=0:1:12
    canales=15-i 
    epi_d (epi_sum > canales)=1
    epi_t=epi_d
    if sum(epi_t)>0
      break
    end
end
epi_t=epi_d'; %transpongo la matriz


%Detección de Inicio y Fin de la epilepsia
start=zeros(size(epi_t));
finish=zeros(size(epi_t));
epi_x=epi_t(1:end-1)-epi_t(2:end);%Busco saltos en mi Matriz
start(epi_x==-1)=1;
[i start]=find(start==1); %Inicio ataque epileptico
finish(epi_x==1)=1;
[i finish]=find(finish==1);%Fin ataque epileptico

%Visualización
figure
[col,fil] = size (epi);
epi2 = epi + repmat([0:hdr.ns-1]*2,col,1); 
plot([0:col-1]*epoch_dur/Fs,epi2);
grid on 
grid minor
xlabel('Tiempo(seg)');
ylabel('Canales');
title('Epilepsia Detectada')
axis([0 3600 -5 55]) 
h1=line([(2*start)-10;(2*start)-10],[min(eeg(:)) max(eeg(:))+600*hdr.ns]); %Marca de Inicio Epilepsia
set(h1,'linewidth',0.5,'color','r','linestyle','-.'); 
h2=line([(2*finish)+10;(2*finish)+10],[min(eeg(:)) max(eeg(:))+600*hdr.ns]);%Marca de Fin Epilepsia
set(h2,'linewidth',0.5,'color','r','linestyle','--'); 

