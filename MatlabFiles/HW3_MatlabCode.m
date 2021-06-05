% This part is used to plot bode diagram of buck converter
% with non-idealities (ESR of capacitor).

% Parameters
L = 10e-6; % H - Inductor
C = 10e-6; % F - Capacitor
Rc = 10e-3; % Ohm - ESR of capacitor
Vin = 5; % V - Input voltage
Vo = 3.3; % V - Output voltage
Io = 3; % A - Output current
RL = Vo/Io; % Ohm - Load resistance
fsw = 200e3; % Hz - Switching freq
Vref = 1.2; % V - Ref voltage
Vosc = 1.8; % V - Oscillator peak voltage for creating duty cycle

num_nonideal = (Vin*RL).*[C*Rc 1]; % Numerator
den_nonideal = [L*C*(RL+Rc) L+RL*C*Rc RL]; % Denominator
G_nonideal = tf(num_nonideal,den_nonideal);

num_ideal = Vin.*[1];
den_ideal = [L*C (L/RL) 1];
G_ideal = tf(num_ideal,den_ideal);

figure(1);
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
opts.Title.String = 'Bode Plot of Buck Converter with and without ESR of Capacitor';
opts.Title.FontSize = 16;
%opts.XLabel.Limit = [1e]
h = bodeplot(G_nonideal,G_ideal,opts);
legend('Non-Ideal','Ideal','Fontsize',16)
set(findall(gcf,'Type','line'),'LineWidth',3)
% legend('Non-Ideal','Ideal')
% figure(1);
% bode(G_nonideal)
% hold on
% bode(G_ideal)
% set(findall(gcf,'Type','line'),'LineWidth',3)
% title('Bode Plot of Buck Converter with and without ESR of Capacitor')
% legend('Non-Ideal','Ideal')
% grid minor;
%% This part is used to find poles and zeros of TF

syms s;

num_eqn = Vin*RL*(C*Rc*s+1);
solve(num_eqn,s);
zero1 = double(ans)/(2*pi);

den_eqn = Vosc*(L*C*(s^2)*(RL+Rc)+s*(L+RL*C*Rc)+RL);
solve(den_eqn,s);
pole1 = abs(double(ans(1)))/(2*pi);
pole2 = abs(double(ans(2)))/(2*pi);

% den_eqn_ideal = L*C*(s^2)+s*(L/RL)+1;
% solve(den_eqn_ideal,s);
% pole1_ideal = abs(double(ans(1)));
% pole2_ideal = abs(double(ans(2)));
%% This part is used to find compensator poles&zeros
fo = (1/5)*fsw;
Fp3 = fsw/2;
teta = pi*(70/180);
Fz2 = fo*sqrt((1-sin(teta))/(1+sin(teta)));
Fp2 = fo*sqrt((1+sin(teta))/(1-sin(teta)));
Fz1 = 0.5*Fz2;


%% This part is used for Type III-B compensator

% Parameters
Cf3 = 2.2e-9; % F
Rf3 = 320; % Ohm
%Rf3 = 1/(2*pi*Cf3*Fp2);
Rf1 = 10e3; % Ohm
%Rf1 = (1/(2*pi*Cf3*Fz2))-Rf3;
Rf2 = 5.69e3; % Ohm
%Rf2 = (Vref/(Vo-Vref))*Rf1;
Rc1 = 4.12e3; % Ohm
%Rc1 = (2*pi*fo*L*C*Vosc)/(Vin*Cf3);
Cc1 = 12e-9; % F
%Cc1 = 1/(2*pi*Rc1*Fz1);
Cc2 = 390e-12; % F
%Cc2 = 1/(2*pi*Rc1*Fp3);

s = tf('s');
num_comp = -(1+Rc1*Cc1*s)*(1+s*Cf3*(Rf1+Rf3));
den_comp = s*Rf1*Cc1*(Rc1*Cc2*s+1)*(1+s*Rf3*Cf3);
H = num_comp/den_comp;

%% This part is used to combine transfer function
% i.e. Getting loop gain
k = Vo/Vref;
%T = H*G_nonideal*(1/k);
T = H*G_nonideal;
figure(2);
bode(G_nonideal)
hold on
bode(T)
%hold on
%bode(H);
set(findall(gcf,'Type','line'),'LineWidth',2)
legend('Non-Ideal','Loop Gain','H')
grid minor;
