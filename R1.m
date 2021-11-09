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

% R1: Reactor 1. Alimentación del biodigestor. En la alimentación Biomasas
% (X1 y X2) son cero. X1in = X2in = 0

function [sys, X0]=R1(t,X,U,flag);

global k1 k2 k3 k4 umax1_1 ks1_1 umax2_1 ks2_1 ki2_1 V Pt RT F

% Return state derivatives (Calculation of derivatives)
if flag==1
    
% Variables de entrada
Qin = U(1);   % Tasa de flujo a la entrada (L/d)
S1in = U(2);  % Materia orgánica a a entrada (gCOD S1/L)
S2in = U(3);  % Acidos grasos volátiles (mmol/L)


% Tasas de crecimiento específico 
mu1 = umax1_1*X(3)/(ks1_1+X(3));                % Materia orgánica DQO - Monod (1/d)
mu2 = umax2_1*X(4)/(ks2_1+X(4)+(X(4)^2)/ki2_1);   % Acidos grasos volátiles - Haldane (1/d)

r(1) = mu1*X(1);   % (gCOD_X1/L*d)
r(2) = mu2*X(2);   % (gCOD_X2/L*d)

% Tasa de Dilución 
D = Qin/V;      % (1/d)

%#################################################
% Ecuaciones diferenciales Reactor de alimentación
%#################################################

dX(1) = (mu1*F - D)*X(1);                               % Biomasa 1 (gCOD_X1/L)
dX(2) = (mu2*F - D)*X(2);                               % Biomasa 2 (gCOD_X2/L)
dX(3) = D*(S1in - X(3)) - k1*mu1*X(1)*F;                % Sustrato 1 - DQO (gCOD/L)
dX(4) = D*(S2in - X(4)) + k2*mu1*X(1)*F - k3*mu2*X(2)*F;% Sustrato 2 - VFA (mmol/L)

    sys = dX';
 
% Return system outputs

elseif flag==3   %%%%  	mdloutputs
    
for i = 1:4
    if X(i) < 0
        X(i) = 0;
    else
        X(i) = X(i);
    end
end

for i = 1:4
    Y(i) = X(i);
end

% Calculo para el flujo de metano Qch4 (L/d) 
    mu2 = umax2_1*X(4)/(ks2_1+X(4)+(X(4)^2)/ki2_1);   % Acidos grasos volátiles - Haldane (1/d)
    kk9 = RT/Pt;       % m3/mol
    qm = k4*mu2*X(2);  % mmol CH4/d.L
	qch4 = V*qm*F;	   % mmol CH4/d
    Qch4 = kk9*qch4;   % Flujo de metano (L/d) 

 Y(5) = max(0,Qch4);   %eliminar value negative
 Y(6)= U(1);
   sys = Y;
    
% Return sizes	
elseif flag==0 %%% mdlinitializesizes

    %[ContStates  DiscStates  Outputs  Inputs  Feedthrough(0/1)  SampleTime(0/1)]
     sys = [4  0  6  3  0  1];  

     % Condiciones iniciales
X0(1) = 1.3;  % Bimasa 1
X0(2) = 1+2;    % Biomasa 2
X0(3) = 2.8+1;  % Sustrato 1
X0(4) = 3+2;    % Sustrato 2

% Continuous system 
else     
	sys=[];
end     