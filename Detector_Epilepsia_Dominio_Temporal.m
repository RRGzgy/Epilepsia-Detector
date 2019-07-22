%-----------------------------------------------------------------
%%%% Detector de Epilepsia en Dominio Temporal
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
hold on 
%onset=1467; %Inicio epilepsia
%endtime=1494; %Fin epilepsia
%h1=line([onset;onset],[min(eeg(:)) max(eeg(:))+600*hdr.ns]); %Marca de Inicio
%set(h1,'linewidth',1,'color','k','linestyle','-.'); %Caracteristicas Marca Inicio
%h2=line([endtime;endtime],[min(eeg(:)) max(eeg(:))+600*hdr.ns]);%Marca de Fin
%set(h2,'linewidth',1,'color','k','linestyle','--'); %Caracteristicas Marca Fin
xlabel('Tiempo(seg)');
ylabel('Canales');
title('Señal EEG')
axis([0 3600 -1000 17000]) %Ajustar Gráfica


%%Filtrado
%Filtro eliminación de ruido Paso Banda (0.5 a 60Hz)
Wp=0.5/(Fs/2); %Frecuencia de Paso
Ws=60/(Fs/2); %Frecuencia de Rechazo
Rp=3; %Ondulación de paso en decibelios
Rs=60; %Atenuación de rechazo en decibelios
[n,Wn] = buttord(Wp,Ws,Rp,Rs); %Calculo del Orden del Filtro
[b,a] = butter(n,Wn); %Filtro Paso Banda
eegpb = filtfilt(b,a,eeg); %Señal después del filtrado Paso Banda

%Filtro eliminación de ruido Notch (50Hz)
Wo = 50/(Fs/2); 
Bw = Wo/35;
[b2,a2] = iirnotch(Wo,Bw); %Filtro Notch
eegnotch = filtfilt(b2,a2,eegpb); %Señal después del filtrado Notch

%Visualización EEG despúes de pasar los filtros
figure;
eegnotch2 = eegnotch + repmat([0:hdr.ns-1]*600, hdr.samples(1)*hdr.records,1);% Separación para la visualización
plot([0:size(eeg,1)-1]/Fs,eegnotch2)
set(gca,'YTick',[0:hdr.ns-1]*600,'YTickLabel',hdr.label,'FontSize',10)
grid on %Activar Cuadricula
grid minor
hold on
xlabel('Tiempo(seg)');
ylabel('Canales');
title('Filtrado Señal EEG ')
axis([0 3600 -1000 17000]) %Ajustar Grafica


%% Estudio de la Energía 
%Energía
epoch_dur=512; % Número de muestras 
nepoch=size(eegnotch,1)/epoch_dur;
ene=reshape(eegnotch,[epoch_dur,nepoch,hdr.ns]);
ene=sum(ene.^2)/nepoch;%Energía de los Epoch
ene=squeeze(ene); 


%%Filtro para eliminar falsas alarmas
%Filtro de Mediana
ums=7; %Ventana mediana
enemed=medfilt1 (ene,ums);%Filtro de Mediana


%%Detección de epilepsia por umbral
%Umbral de detección
%% I Forma
umbral = floor(prctile(enemed(:), 95));%Percentil 95%
sumb = zeros(size(enemed));
sumb(enemed > umbral) = 1;

%Visualización
%figure
%sumb2 = sumb + repmat([0:hdr.ns-1]*2,nepoch,1); %Separación para la visualización
%plot([0:nepoch-1]*epoch_dur/Fs,sumb2);
%grid on 
%grid minor
%hold on 
%xlabel('Tiempo(seg)');
%ylabel('Canales');
%title('Umbral de Energía para Señal EEG 1')
%axis([0 3600 -5 55])


%% II Forma
%I = sort(enemed);
%umbral2 = (I(floor(size(enemed,1)*0.95),:));
%sumb=bsxfun(@gt,enemed,umbral2);

%figure
%sumb2 = sumb + repmat([0:hdr.ns-1]*2,nepoch,1); %Separación para la visualización
%plot([0:nepoch-1]*epoch_dur/Fs,sumb2);
%grid on %Activar Cuadricula
%grid minor
%hold on %Mantener Figura activa
%xlabel('Tiempo(seg)');
%ylabel('Canales');
%title('Umbral de Energía para Señal EEG 2')
%axis([0 3600 -5 55])

%% III Forma
%[p,n]=size(enemed);
%I = sort(enemed);
%umbral2 = (I(floor(size(enemed,1)*0.95),:))
%sumb3 = zeros(size(enemed));
%for j=1:n
%sumb3(enemed(:,j) > umbral2(1,j)) = 1;
%end


%Filtro Temporal
umbral_t = 15; %Umbral temporal
sconv = conv2(sumb,ones(umbral_t,1)); %Convolución de la señal en los dos umbrales
epi = zeros(size(sconv));
epi(sconv >= umbral_t) = 1;


%Buscar Canales Activos
epi_sum = sum(epi,2);%Sumo todos los canales para ver los que estan activos
epi_d = zeros(size(enemed(:,1)));%Hago una matriz de ceros

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

