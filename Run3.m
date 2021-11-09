%----------------------------------------------------------------------------------------------------
% MODELO CINÉTICO BIODIGESTOR TIPO BOLSA
% Corrida de simunlación de un biodigestor tipo bolsa
% La simuación es basada en la división del biodigestor en varios reactores agitados 
% Universidad de Ibaguè - Colombia
% Autores:
% Liliana Delgadillo Mirquez; liliana.delgadillo@unibague.edu.co
% Julian Alberto Avila Navarro; juian.avila@unibague.edu.co
% Laura Daniela Cardona Acuña: 2320162088@estudiantesunibague.edu.co
% Ibagué, Mayo 2021
% %----------------------------------------------------------------------------------------------------

% Run: Inicia la simulación. Carga el programa con las variables y
% parametros cinéticos y estequiometricos necesarios. Al mismo tiempo
% incluye los datos eperimentales. Run3 abre el simulink Biodigetor_3 

clear all;
close all;
clc;

global k1 k2 k3 k4 umax1_1 ks1_1 umax2_1 ks2_1 ki2_1 umax1_2 ks1_2 umax2_2 ks2_2 ki2_2 umax1_3 ks1_3 umax2_3 ks2_3 ki2_3 V Pt RT F
    
%#############################
% Coeficientes de rendimiento 
%#############################

k1 = 1-0.3;%0.5;		%Rendimiento S1->X1  (gCOD/g X1)
k2 = 2-0.1;%0.6;		%Rendimiento S1->S2  (mmol S2/g X1)	
k3 = 10-7;%0.8; 		%Rendimiento S2->X2  (mmol S2/ g X2)
k4 = 0.17+0.1;   		%Rendimiento S2->CH4 (mmol CH4/g X2)

%############################
% Parámetros estequiométricos
%############################
% Reactor 1
umax1_1  = 1.2-0.0;%1.2;		%Máxima tasa especifica de crecimiento para X1 (d-1)
ks1_1    = 55+0;%7.1;       %Constante de saturación media para S1 (g/L)
umax2_1  = 0.75;%0.74;      %Máxima tasa especifica de crecimiento para X2 	(d-1)
ks2_1    = 28;%9.28;      %Constante de saturación media para S2 (mmol/L)
ki2_1    = 50+0;%256;		%Constante de inhibición para S2 (mmol/L)
% Reactor 2
umax1_2  = 1.9+0;%1.2;		%Máxima tasa especifica de crecimiento para X1 (d-1)
ks1_2    = 59+0;%7.1;       %Constante de saturación media para S1 (g/L)
umax2_2  = 0.75;%0.74;      %Máxima tasa especifica de crecimiento para X2 	(d-1)
ks2_2    = 28;%9.28;      %Constante de saturación media para S2 (mmol/L)
ki2_2    = 50+100;%256;		%Constante de inhibición para S2 (mmol/L)
% Reactor 3
umax1_3  = 1.7;%1.2;		%Máxima tasa especifica de crecimiento para X1 (d-1)
ks1_3    = 38;%7.1;       %Constante de saturación media para S1 (g/L)
umax2_3  = 0.75;%0.74;      %Máxima tasa especifica de crecimiento para X2 	(d-1)
ks2_3    = 28;%9.28;      %Constante de saturación media para S2 (mmol/L)
ki2_3    = 50;%256;		%Constante de inhibición para S2 (mmol/L)

% #######################
% Datos experimentales
% #######################

data = xlsread('DatosExp1.xlsx'); % Carga los datos experimentales desde excel
DATAin = [data(:,1),data(:,2)*1000,data(:,3)/1000,data(:,4)]; % Tiempo, Qin, S1in, S2in


%#########################################
% Características del biodigestor de bolsa
%#########################################

V = 85*1000;        % Volumen del bodigestor de bolsa (L)
n = 3;              % Numero de reactores en que se divide el biodigestor
F = 1/n;            % Fraccion de Volumen de cada reactor (adimensional)

%############################################
% Constantes realcionadas con la fase gaseosa
%############################################

T = 273 + 28;      % Temperatura en (K). T = 273 + °C
RT = 8.32*T/100;   % Producto Temp x Const. (Pa.m3)/(K.mol)*K = (Pa.m3)/mol
Vm = 8.32*(273+35);% Volumen molar del gas (T=38C;P=1Patm) 	(L/mol)
Pt = 1.013; %1.013 * (1 + P_dome); % Presión absoluta en la cabeza delreactor (Pa)


open('Biodigestor_3');
sim(gcs);

% ###############
% ### FIGURAS ###
% ###############

figure(1)
%Plot carga orgánica de entrada
subplot(1,3,1)
plot(data(:,1),DATAin(:,2).*DATAin(:,3)); hold on;    %in blue: measurements from the flowmeter
grid on;
hold off;
xlabel('Tiempo (d)')
ylabel('Carga orgánica (g/d)')

%Plot de sustrato s1 
subplot(1,3,2)
plot(tout(:,:),out1(:,3),'b'); hold on;
grid on;
plot(tout(:,:),out2(:,3),'r');
plot(tout(:,:),out3(:,3),'k');
plot(data(:,1),data(:,5)/1000,'-or'); 
hold off;
legend('R1','R2','R3','S1exp')
xlabel('Tempo (d)')
ylabel('Materia orgánica (gCOD/L)')

%Plot de biomasa X1 
subplot(1,3,3)
plot(tout(:,:),out1(:,1),'b'); hold on;
grid on;
plot(tout(:,:),out2(:,1),'r');
plot(tout(:,:),out3(:,1),'k');
hold off;
legend('R1','R2','R3')
xlabel('Tiempo (d)')
ylabel('X1 (gCOD/L)')

figure(2)
%Plot carga orgánica de entrada
subplot(1,3,1)
plot(data(:,1),DATAin(:,2).*DATAin(:,3)); hold on;    %in blue: measurements from the flowmeter
grid on;
hold off;
xlabel('Tiempo (d)')
ylabel('Carga orgánica (g/d)')

%Plot de sustrato s2 
subplot(1,3,2)
plot(tout(:,:),out1(:,4),'b'); hold on;
grid on;
plot(tout(:,:),out2(:,4),'r');
plot(tout(:,:),out3(:,4),'k');
hold off;
legend('R1','R2','R3')
xlabel('Tiempo (d)')
ylabel('AGV (mmol/L)')

%Plot de biomasa X2 
subplot(1,3,3)
plot(tout(:,:),out1(:,2),'b'); hold on;
grid on;
plot(tout(:,:),out2(:,2),'r');
plot(tout(:,:),out3(:,2),'k');
hold off;
legend('R1','R2','R3')
xlabel('Tiempo (d)')
ylabel('X2 (gCOD/L)')


figure(3)
%Plot carga organica de entrada
subplot(1,2,1)
plot(data(:,1),DATAin(:,2).*DATAin(:,3)); hold on;    %in blue: measurements from the flowmeter
grid on;
hold off;
xlabel('Tiempo (d)')
ylabel('Carga organica (g/d)')

%Plot de flujo de metano 
subplot(1,2,2)
plot(tout(:,:),out1(:,5),'b'); hold on;
grid on;
plot(tout(:,:),out2(:,5),'r'); 
plot(tout(:,:),out3(:,5),'k');
hold off;
legend('R1','R2','R3')
xlabel('Tiempo (d)')
ylabel('CH4 (L/d)')

figure(4)
%Plot flujo de entrada
subplot(1,2,1)
plot(data(:,1),data(:,3)/1000,'b'); hold on;    %in blue: measurements from the flowmeter
plot(data(:,1),data(:,5)/1000,'-or'); 
grid on;
hold off;
legend('S1 in','S1 out')
xlabel('Tiempo (d)')
ylabel('Carga ogànica (g/d)')

%Plot de flujo de metano 
subplot(1,2,2)
plot(tout(:,:),out3(:,3),'k'); hold on;
plot(data(:,1),data(:,5)/1000,'-*r'); 
hold off;
legend('S1sim','S1exp')
xlabel('Tiempo (d)')
ylabel('DQO (g/L)')
grid


outR3L=out3;
outR23L=out2;
outR13L=out1;
TR3L=tout;
DATAR3L=data