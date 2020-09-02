function B = ThresholdEqn(parameters,t,X,U)
% ThresholdEqn   Compute whether the battery has reached end-of-discharge
%
%   B = ThresholdEqn(parameters,t,X,U) computes whether the battery model
%   for the given time, states, and inputs, has crossed the voltage
%   threshold defining end-of-discharge. 
%
%   Copyright (c) 2016 United States Government as represented by the
%   Administrator of the National Aeronautics and Space Administration.
%   No copyright is claimed in the United States under Title 17, U.S.
%   Code. All Other Rights Reserved.

% Extract states
Tb = X(1,:);
Vo = X(2,:);
Vsn = X(3,:);
Vsp = X(4,:);
qnB = X(5,:);
qnS = X(6,:);
qpB = X(7,:);
qpS = X(8,:);

% Extract inputs
P = U(1,:);

% Constraints
xpS = qpS./parameters.qSMax;
Vep3 = parameters.Ap3.*((2.*xpS-1).^(3+1) - (2.*xpS.*3.*(1-xpS))./(2.*xpS-1).^(1-3))./parameters.F;
Vep8 = parameters.Ap8.*((2.*xpS-1).^(8+1) - (2.*xpS.*8.*(1-xpS))./(2.*xpS-1).^(1-8))./parameters.F;
Vep6 = parameters.Ap6.*((2.*xpS-1).^(6+1) - (2.*xpS.*6.*(1-xpS))./(2.*xpS-1).^(1-6))./parameters.F;
Vep5 = parameters.Ap5.*((2.*xpS-1).^(5+1) - (2.*xpS.*5.*(1-xpS))./(2.*xpS-1).^(1-5))./parameters.F;
Vep10 = parameters.Ap10.*((2.*xpS-1).^(10+1) - (2.*xpS.*10.*(1-xpS))./(2.*xpS-1).^(1-10))./parameters.F;
Vep9 = parameters.Ap9.*((2.*xpS-1).^(9+1) - (2.*xpS.*9.*(1-xpS))./(2.*xpS-1).^(1-9))./parameters.F;
Vep12 = parameters.Ap12.*((2.*xpS-1).^(12+1) - (2.*xpS.*12.*(1-xpS))./(2.*xpS-1).^(1-12))./parameters.F;
Vep4 = parameters.Ap4.*((2.*xpS-1).^(4+1) - (2.*xpS.*4.*(1-xpS))./(2.*xpS-1).^(1-4))./parameters.F;
Vep11 = parameters.Ap11.*((2.*xpS-1).^(11+1) - (2.*xpS.*11.*(1-xpS))./(2.*xpS-1).^(1-11))./parameters.F;
Vep2 = parameters.Ap2.*((2.*xpS-1).^(2+1) - (2.*xpS.*2.*(1-xpS))./(2.*xpS-1).^(1-2))./parameters.F;
Vep7 = parameters.Ap7.*((2.*xpS-1).^(7+1) - (2.*xpS.*7.*(1-xpS))./(2.*xpS-1).^(1-7))./parameters.F;
Vep0 = parameters.Ap0.*((2.*xpS-1).^(0+1))./parameters.F;
Vep1 = parameters.Ap1.*((2.*xpS-1).^(1+1) - (2.*xpS.*1.*(1-xpS))./(2.*xpS-1).^(1-1))./parameters.F;
xnS = qnS./parameters.qSMax;
Ven5 = parameters.An5.*((2.*xnS-1).^(5+1) - (2.*xnS.*5.*(1-xnS))./(2.*xnS-1).^(1-5))./parameters.F;
Ven1 = parameters.An1.*((2.*xnS-1).^(1+1) - (2.*xnS.*1.*(1-xnS))./(2.*xnS-1).^(1-1))./parameters.F;
Ven10 = parameters.An10.*((2.*xnS-1).^(10+1) - (2.*xnS.*10.*(1-xnS))./(2.*xnS-1).^(1-10))./parameters.F;
Ven7 = parameters.An7.*((2.*xnS-1).^(7+1) - (2.*xnS.*7.*(1-xnS))./(2.*xnS-1).^(1-7))./parameters.F;
Ven2 = parameters.An2.*((2.*xnS-1).^(2+1) - (2.*xnS.*2.*(1-xnS))./(2.*xnS-1).^(1-2))./parameters.F;
Ven8 = parameters.An8.*((2.*xnS-1).^(8+1) - (2.*xnS.*8.*(1-xnS))./(2.*xnS-1).^(1-8))./parameters.F;
Ven4 = parameters.An4.*((2.*xnS-1).^(4+1) - (2.*xnS.*4.*(1-xnS))./(2.*xnS-1).^(1-4))./parameters.F;
Ven3 = parameters.An3.*((2.*xnS-1).^(3+1) - (2.*xnS.*3.*(1-xnS))./(2.*xnS-1).^(1-3))./parameters.F;
Vep = parameters.U0p + parameters.R.*Tb./parameters.F.*log((1-xpS)./xpS) + Vep0 + Vep1 + Vep2 + Vep3 + Vep4 + Vep5 + Vep6 + Vep7 + Vep8 + Vep9 + Vep10 + Vep11 + Vep12;
Ven0 = parameters.An0.*((2.*xnS-1).^(0+1))./parameters.F;
Ven11 = parameters.An11.*((2.*xnS-1).^(11+1) - (2.*xnS.*11.*(1-xnS))./(2.*xnS-1).^(1-11))./parameters.F;
Ven12 = parameters.An12.*((2.*xnS-1).^(12+1) - (2.*xnS.*12.*(1-xnS))./(2.*xnS-1).^(1-12))./parameters.F;
Ven6 = parameters.An6.*((2.*xnS-1).^(6+1) - (2.*xnS.*6.*(1-xnS))./(2.*xnS-1).^(1-6))./parameters.F;
Ven9 = parameters.An9.*((2.*xnS-1).^(9+1) - (2.*xnS.*9.*(1-xnS))./(2.*xnS-1).^(1-9))./parameters.F;
Ven = parameters.U0n + parameters.R.*Tb./parameters.F.*log((1-xnS)./xnS) + Ven0 + Ven1 + Ven2 + Ven3 + Ven4 + Ven5 + Ven6 + Ven7 + Ven8 + Ven9 + Ven10 + Ven11 + Ven12;
V = Vep - Ven - Vo - Vsn - Vsp;
Vm = V;

% Return true if voltage is less than the voltage threshold
B = Vm < parameters.VEOD;
