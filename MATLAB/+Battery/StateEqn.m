function XNew = StateEqn(parameters,t,X,U,N,dt)
% StateEqn   Compute the new states of the battery model
%
%   XNew = StateEqn(parameters,t,X,U,N,dt) computes the new states of the
%   battery model given the parameters strcucture, the current time, the
%   current states, inputs, process noise, and the sampling time.
%
%   Copyright (c)Â 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Extract states


% 
% Tb = X(1,:);
% Vo = X(2,:);
% Vsn = X(3,:);
% Vsp = X(4,:);
% qnB = X(5,:);
% qnS = X(6,:);
% qpB = X(7,:);
% qpS = X(8,:);

XNew = OriginalModel(parameters,t,X,U,N,dt);





function XNew = OriginalModel(parameters,t,X,U,N,dt)
% Extract inputs
P = U(1,:);

% Thermal ECM paramters
   
[Tb,Vo,Vsn,Vsp,qnB,qnS,qpB,qpS,xnS,xpS,Ven,Vep,V] = Battery.CurrentState(parameters,t,X,U);

% Constraints

CnBulk = qnB./parameters.VolB;
CnSurface = qnS./parameters.VolS;
xSn = qnS./parameters.qSMax;
CpSurface = qpS./parameters.VolS;
% xnS = qnS./parameters.qSMax;
xSp = qpS./parameters.qBMax;
% xpS = qpS./parameters.qSMax;
CpBulk = qpB./parameters.VolB;
qdotDiffusionBSn = (CnBulk-CnSurface)./parameters.tDiffusion;
qnBdot = - qdotDiffusionBSn;
Jn0 = parameters.kn.*(1-xSn).^parameters.alpha.*(xSn).^parameters.alpha;
qdotDiffusionBSp = (CpBulk-CpSurface)./parameters.tDiffusion;
Jp0 = parameters.kp.*(1-xSp).^parameters.alpha.*(xSp).^parameters.alpha;
qpBdot = - qdotDiffusionBSp;
i = P; % P./V; Assuming constant current
qpSdot = i + qdotDiffusionBSp;
Jn = i./parameters.Sn;
VoNominal = i.*parameters.Ro;
Jp = i./parameters.Sp;
qnSdot = qdotDiffusionBSn - i;
VsnNominal = parameters.R.*Tb./parameters.F./parameters.alpha.*asinh(Jn./(2.*Jn0));
Vodot = (VoNominal-Vo)./parameters.to;
VspNominal = parameters.R.*Tb./parameters.F./parameters.alpha.*asinh(Jp./(2.*Jp0));
Vsndot = (VsnNominal-Vsn)./parameters.tsn;
Vspdot = (VspNominal-Vsp)./parameters.tsp;
U = Vep - Ven; 
voltage_eta = U - V;

% etakin = (Vsn) + (Vsp);
% Ptherm = i.*U;
% Pkin = i.*etakin;
% Pohm = i.*VoNominal;
h = 5;
mC = 37.04; % 0.045kg � 823J/kg/C = 37.04 kg/m2/(K-s^2)
tau = 50; % 0.045kg � 823J/kg/C/8.84cm2/5W/K-m2 = 8379 s
Rshort = 1;
% Ptherm + 
% Tbdot = (Ptherm + Pohm + Pkin + h.*(parameters.x0.Tb - Tb))./(7.5e6*0.35); %- i.*T.*deltaU./deltaT;
% Tbdot = ((voltage_eta.*i./mC) + (i.^2.*Rshort./mC) + ((parameters.x0.Tb - Tb)./tau)); %- i.*T.*deltaU./deltaT; %Newman
Tbdot = ((voltage_eta.*i./mC) + ((parameters.x0.Tb - Tb)./tau)); %- i.*T.*deltaU./deltaT; %Newman
% Update state
XNew = zeros(size(X));
XNew(1,:) = Tb + Tbdot*dt;
XNew(2,:) = Vo + Vodot*dt;
XNew(3,:) = Vsn + Vsndot*dt;
XNew(4,:) = Vsp + Vspdot*dt;
XNew(5,:) = qnB + qnBdot*dt;
XNew(6,:) = qnS + qnSdot*dt;
XNew(7,:) = qpB + qpBdot*dt;
XNew(8,:) = qpS + qpSdot*dt;

% Add process noise
XNew = XNew + dt*N;

function XNew = LiPoModel(parameters,t,X,U,N,dt)

error("Not implemented");






% Pout = i.*eta;
% Tbdot = (Pout + h.*(parameters.x0.Tb - Tb))./(830.*80); %- i.*T.*deltaU./deltaT;
% XNew = zeros(size(X));
% XNew(1,:) = Tb + Tbdot*dt;
% XNew(2,:) = Vo + Vodot*dt;
% XNew(3,:) = Vsn + Vsndot*dt;
% XNew(4,:) = Vsp + Vspdot*dt;
% XNew(5,:) = qnB + qnBdot*dt;
% XNew(6,:) = qnS + qnSdot*dt;
% XNew(7,:) = qpB + qpBdot*dt;
% XNew(8,:) = qpS + qpSdot*dt;

% 
% Tbold =  X(1,:);
% Tbnew =  XNew(1,:);
% deltaT = Tbnew - Tbold;
% [Vennew,Vepnew] = Unew(parameters,XNew);
% E = Vepnew - Vennew;
% deltaU = E - E0;
% 
% deltaU_deltaT = deltaU./deltaT;
% if isinf(deltaU_deltaT)
%     deltaU_deltaT = 0;
% end
% 
% % Ptherm + 
% Tbdot = (Pout + h.*(parameters.x0.Tb - Tb))./(1) - i.*Tb.*deltaU_deltaT;
% % Update state



function [Ven,Vep] = Unew(parameters,X)
    Tb = X(1,:);
    qnS = X(6,:);
    qpS = X(8,:);
    xnS = qnS./parameters.qSMax;
    xpS = qpS./parameters.qSMax;
    Ven = Ven_func(parameters,xnS,Tb);
    Vep = Vep_func(parameters,xpS,Tb);
    U = Vep - Ven; 
   
    
    