function [weighPSD] = ISOseatVER(my_psd, my_freq)

    % ISOseatVER Ponderación (en frecuencia) de la señal mediante la NTC
    % 5436-1. Para obtener la densidad espectral de la señal ponderada

    %   weighPSD = ISOseatVER(Y,f); donde   f = vector de frecuencia, 
    %                                       Y = Estimación de la densidad espectral de potencia de Welch.

    % La función retorna el PSD ponderado mediante las magnitudes en dB de las
    % funciones de transferencia que definen el filtro Wk.

    % Parámetros de la función de transferencia. Tabla A.1 NTC 5436-1.
    f1 = 0.4;
    f2 = 100;
    f3 = 12.5;
    f4 = 12.5;
    f5 = 2.37;
    f6 = 3.35;
    q4 = 0.63;
    q5 = 0.91;
    q6 = 0.91;

    % Magnitudes de las funciones de transferencia. Eqs A.1, A.2, A.3 y A.4 de
    % la NTC 5436-1.
    Hh = ((my_freq.^4)./(my_freq.^4+f1^4)).^0.5;
    Hl = ((f2^4)./(my_freq.^4+f2^4)).^0.5;
    Ht = ((my_freq.^2+f3^2)/(f3^2)).^0.5.*((f4^4*q4^2)./(my_freq.^4*q4^2+my_freq.^2*f4^2*(1-2*q4^2)+f4^4*q4^2)).^0.5;
    Hs = (q6/q5)*((my_freq.^4*q5^2+my_freq.^2*f5^2*(1-2*q5^2)+f5^4*q5^2)./(my_freq.^4*q6^2+my_freq.^2*f6^2*(1-2*q6^2)+f6^4*q6^2)).^0.5;
    W  = (Hh.*Hl.*Ht.*Hs);

    % Ponderación de la señal. Densidad espectral.
    weighPSD = abs(my_psd).*(W.^2);

end