pkg load control;

%%%%%% fizioloske konstante %%%%%%
SG_const = 0.014;
C_const = 0.1725;
ka = 6*10^(-6);
kb = 0.01;

p = [SG_const, C_const, ka, kb];
%%%%%% fizioloske konstante %%%%%%

I0 = 15;
x20 = I0 * p(3)/p(4);
x1 = p(2) / (p(1) + x20);

pol_1 = p(1) + x20;

s = tf("s");
G = -x1*p(3)/((s+p(4))*(s+pol_1)); % Funkcija prenosa

%%%%%% PI Regulator %%%%%%
Tao = 150;
K = dcgain(G);
%K = -x1*p(3)/(p(4)*pol_1);
T1 = abs(1/pol_1);
T2 = abs(1/p(4));

Ti = min(T1, T2);

Kp = T1/(abs(K)*(Tao+T2));
Ki = 3*Kp/Ti;

G_r = -(Kp*s + Ki)/s;
W = series(G_r, G);

nyquist(W)
[d, Phi, Wpi, Wpf] = margin(W);
set(gca, 'FontSize', 24, 'LineWidth', 1.5);
xlim([-1.5, 0.5]);
ylim([-1, 1]);
xlabel("Realna osa"); ylabel("Imaginarna osa");
title("Nikvistova kriva");
disp(["Pretek pojacanja: ", num2str(d)]);
disp(["Pretek faze: ", num2str(Phi)]);
