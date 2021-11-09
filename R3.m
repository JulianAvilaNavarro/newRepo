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

% R3: Reactor 3. Reactor final. 

function [sys, X0]=R3(t,X,U,flag);

% ODE Implementation Model No4. 
% With hydrolisis of particulate

global k1 k2 k3 k4 umax1_3 ks1_3 umax2_3 ks2_3 ki2_3 V Pt RT F

% Return state derivatives (Calculation of derivatives)
if flag==1

% Variables de entrada
X1in = U(1);    % Biomasa acidogénica (gCOD X1/L)
X2in = U(2);    % Biomasa metanogénica (gCOD X2/L)
S1in = U(3);    % Materia orgánica a a entrada (gCOD S1/L)
S2in = U(4);    % Acidos grasos volátiles (mmol/L)
Qch4_R1 = U(5); % Flujo de metano del reactor anterior (L/d)
Qin = U(6);     % Tasa de flujo a la entrada (L/d)


% Tasas de crecimiento específico 
mu1 = umax1_3*X(3)/(ks1_3+X(3));                % Materia orgánica DQO - Monod (1/d)
mu2 = umax2_3*X(4)/(ks2_3+X(4)+(X(4)^2)/ki2_3);   % Acidos grasos volátiles - Haldane (1/d)

r(1) = mu1*X(1);   % (gCOD_X1/L*d)
r(2) = mu2*X(2);   % (gCOD_X2/L*d)

% Tasa de Dilución 
D = Qin/V;      % (1/d)

%#################################################
% Ecuaciones diferenciales Reactor de alimentación
%#################################################

dX(1) = D*(X1in - X(1)) + mu1*X(1)*F;           % Biomass 1 (gCOD_X1/L)
dX(2) = D*(X2in - X(2)) + mu2*X(2)*F;           % Biomass 2 (gCOD_X2/L)
dX(3) = D*(S1in - X(3)) - k1*mu1*X(1)*F;        % Sustrato 1 - DQO (gCOD/L)
dX(4) = D*(S2in - X(4)) + k2*mu1*X(1)*F - k3*mu2*X(2)*F;  % Sustrato 2 - VFA (mmol/L)

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
mu2 = umax2_3*X(4)/(ks2_3+X(4)+(X(4)^2)/ki2_3);   % Acidos grasos volátiles - Haldane (1/d)
    kk9 = RT/Pt;       % m3/mol
    qm = k4*mu2*X(2);  % mmol CH4/d.L
	qch4 = V*qm*F;	   % mmol CH4/d
    Qch4 = kk9*qch4;   % Flujo de metano (L/d) 

 Y(5) = max(0,Qch4) + U(5);   %eliminar value negative
 Y(6)= U(6);
   sys = Y;
    
% Return sizes	
elseif flag==0 %%% mdlinitializesizes

    %[ContStates  DiscStates  Outputs  Inputs  Feedthrough(0/1)  SampleTime(0/1)]
     sys = [4  0  6  6  0  1];  

X0(1) = 4.2;         % Biomasa 1
X0(2) = 1.9;         % Biomasa 2
X0(3) = 1;           % Sustrato 1
X0(4) = 0.5;         % Sustrato 2


% Continuous system 
else     
	sys=[];
end     