clear all; close all; clc

% Cargamos las datos de aceleraci�n. Hay 4 filas distribuidas de la
% siguiente manera:
%       1. Vector de tiempo en s.
%       2. Aceleraci�n longitudinal (eje x) en g's.
%       3. Aceleraci�n vertical (eje y) en g's.
%       4. Aceleraci�n horizantal (eje z) en g's.
load ('ABHigh_Seatpost_Accel.mat')
t  = ADC(1,:)';
ax = ADC(2,:)';
ay = ADC(3,:)';
az = ADC(4,:)';

% Realizamos el procesamiento de los datos, centramos en 0 las
% aceleraciones y convertimos a unidades del sistema internacional.
ax = (ax-mean(ax))*9.81;
ay = (ay-mean(ay))*9.81;
az = (az-mean(az))*9.81;

% Obtenemos la frecuencia de muestreo del sensor. Miramos cu�nto tiempo se
% demora en tomar dos datos contiguos.
fs=1/mean(diff(t));

%% Procesamiento en el dominio del tiempo

% Frecuencias de la funci�n de transferencia Wk de las ponderaciones de las
% frecuencias principales - Tabla A.1 NTC 5436-1
freqs = [0.4 100 12.5 12.5 2.37 3.35];
ws = 2*pi*freqs;
w1 = ws(1); w2 = ws(2); w3 = ws(3); w4 = ws(4); w5 = ws(5); w6 = ws(6);

% Factores de calidad resonante - Tabla A.1 NTC 5436-1
Qs = [0.63 0.91 0.91];
Q4 = Qs(1); Q5 = Qs(2); Q6 = Qs(3);

% Funciones de transferencia: transici�n a-v y etapa ascedente
Ht = tf([0,1/w3,1],[(1/w4)^2,1/(Q4*w4),1]);
Hs = tf([1,w5/Q5,w5^2],[1,w6/Q6,w6^2]);
Hts=Ht*Hs;                              % Teorema de convoluci�n
[bts,ats] = tfdata(Hts,'v');
% M�todo de transformaci�n bilineal para la conversi�n de filtros 
% an�logos a digitales
[Bts,Ats] = bilinear(bts,ats,fs);
ats = filter(Bts,Ats,ay);               % Filtramos ay con la funci�n de transferencia
                                        % cuyo denominador es Ats y
                                        % numerador Bts

% Funciones de transferencia: paso alto y paso bajo
Hl = tf([0,0,1],[(1/w2)^2,sqrt(2)/w2,1]);
Hh = tf([1,0,0],[1,sqrt(2)*w1,w1^2]);
Hlh=Hl*Hh;                              % Limitaci�n en banda
[blh,alh] = tfdata(Hlh,'v');
[Blh,Alh] = bilinear(blh,alh,fs);
aw = filter(Blh,Alh,ats);               % Filtramos ats con la funci�n de transferencia
                                        % cuyo denominador es Alh y
                                        % numerador Blh

% C�lculo de la aceleraci�m rms filtrada (en el dominio del tiempo)
a_rms_time = sqrt(trapz(t,aw.^2)/t(end));
fprintf('La aceleraci�n rms es = %1.2f m/s^2\n', aw)

% Factor de cresta: M�dulo de la relaci�n entre el valor pico instant�neo
% m�ximo de la se�al de aceleraci�n de frecuencia ponderada y su valor rms.
% El factor de cresta NO indics la severidad de la vibraci�n.
CF = max(abs(aw))/a_rms_time;

% Dosis de vibraci�n a la cuarta potencia
VDV = (trapz(t,aw.^4))^(1/4);
fprintf('La dosis de vibraci�n a la cuarta potencia es = %1.2f m/s^1.75\n', VDV)

%% Procesamiento en el dominio de la frecuencia

% Estimaci�n de la densidad espectral de potencia de Welch. En el caso 
% de la densidad espectral, el resultado debe ser escalado por la frecuencia 
% de muestreo, que no se realiza por psd.
[Y,f] = pwelch(ay,window(@rectwin,round(fs)),[],[],fs);

% Ponderaci�n de la densidad espectral por medio de las magnitudes de las
% funciones de transferencia en dB. Es decir, densidad espectral de
% potencia de la se�al ponderada.
weighPSD = ISOseatVER(Y,f);

% C�lculo de la aceleraci�n rms en el dominio de la frecuencia.
a_rms_freq = sqrt(trapz(f,weighPSD));