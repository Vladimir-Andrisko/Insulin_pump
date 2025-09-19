pkg load control;

graphics_toolkit("fltk");

%%%%%% fizioloske konstante %%%%%%
SG_const = 0.014;
C_const = 0.1725;
ka = 6*10^(-6);
kb = 0.01;
%%%%%% fizioloske konstante %%%%%%

t = 1:0.1:500;
zakasnjenje = 60;
delta_I = 5;

%%%%%% Nelinearni %%%%%%
I0 = 15;
X0 = I0 * ka/kb;
G0 = C_const / (SG_const + X0);
x0 = [G0, X0];


I = @(t) 15 + delta_I*(t >= zakasnjenje); % BOLUS doza insulina

d_stanja = @(t, x) [
    -(SG_const + x(2))*x(1) + C_const;
    ka*I(t) - kb*x(2);
];

[tl, x_nelin] = ode45(d_stanja, t, x0);

plot(tl, x_nelin(:, 1), 'b-', 'LineWidth', 2);
%%%%%% Nelinearni %%%%%%

%%%%%% Linearni %%%%%%
u0 = I0;
x20 = u0 * ka/kb;
x10 = C_const / (SG_const + x20);

A = [-(SG_const + x20) -x10; 0 -kb];
B = [0; ka];
C = [1, 0];
D = 0;

sys = ss(A, B, C, D);
G = tf(sys);

u = zeros(size(t));
u(t>=zakasnjenje) = delta_I;

g_lin = lsim(G, u, t) + x10;
hold on;
plot(t, g_lin, 'r--', 'LineWidth', 1.5);
set(gca, 'FontSize', 24, 'LineWidth', 1.5);
xlim([0, 500]);
ylim([6, 8]);
xlabel("Vreme [min]");
ylabel("Glukoza [mmol/L]");
legend("Nelinearan model", "Linearan model");
grid on;
title("Model glukoze u krvi");

%%%%%% Linearni %%%%%%


