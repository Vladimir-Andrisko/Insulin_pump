pkg load control;

graphics_toolkit("fltk");

%%%%%% fizioloske konstante %%%%%%
SG_const = 0.014;
C_const = 0.1725;
ka = 6*10^(-6);
kb = 0.01;

P = [SG_const, C_const, ka, kb];
%%%%%% fizioloske konstante %%%%%%


% funkcija koja vraca izracunate vrendnosti na osnovu parametara
function [g_abs, I_abs, e_dev] = sistem(I0, r, p, t);
  % radna tacka
  x20 = I0 * p(3)/p(4);
  x1 = p(2) / (p(1) + x20);

  pol_1 = p(1) + x20;

  s = tf("s");
  G = -x1*p(3)/((s+p(4))*(s+pol_1)); % Funkcija prenosa

  %%%%%% PI Regulator %%%%%%
  Tao = 30;
  K = dcgain(G)
  %K = -x1*p(3)/(p(4)*pol_1);
  T1 = abs(1/pol_1);
  T2 = abs(1/p(4));

  Ti = min(T1, T2);

  Kp = T1/(abs(K)*(Tao+T2))
  Ki = 1.1*Kp/Ti

  G_r = -(Kp*s + Ki)/s;
  G_rg = series(G_r, G);

  G_sp = feedback(G_rg, 1); % Funkcija spregnutog prenosa
  u = (r-x1) * ones(size(t));

  g_lin = lsim(G_sp, u, t);
  g_abs = g_lin + x1;

  e_dev = u - g_lin;
  [u_dev, _] = lsim(G_rg, e_dev, t);
  I_abs = -u_dev+I0;

end

%%%%%% Simulacija %%%%%%

h = 1;
t = (0:1:3600*h)';
ref = 5.5;      % zeljena vrednost

ref_vector = ref * ones(size(t));

[g_abs1, I_abs1, e_dev1] = sistem(16, ref, P, t);
[g_abs2, I_abs2, e_dev2] = sistem(15, ref, P, t);
[g_abs3, I_abs3, e_dev3] = sistem(14, ref, P, t);
[g_abs4, I_abs4, e_dev4] = sistem(0, ref, P, t);

%%%%%% Plot %%%%%%
lw = 1.5;

figure(1);
plot(t, g_abs1, 'b-', 'LineWidth', lw)
hold on;
plot(t, g_abs2, 'r-', 'LineWidth', lw);
plot(t, g_abs3, 'g-', 'LineWidth', lw);
plot(t, g_abs4, 'm-', 'LineWidth', lw);
plot(t, ref_vector, 'k-.', 'LineWidth', 1);
plot(t, 4*ones(size(t)),'c--', 'LineWidth', 1.5);
plot(t, 6*ones(size(t)),'c--', 'LineWidth', 1.5);

set(gca, 'FontSize', 24, 'LineWidth', lw);
xlim([0, 1000]);
ylim([2, 13]);
xlabel("Vreme [min]");
ylabel("Glukoza [mmol/L]");
title("Glukoza u krvi sa PI regulatorom");
legend("I_0 = 16", "I_0 = 15", "I_0 = 14", "I_0 = 0" ,"r = 5.5", "opseg (4, 6)", 'location', 'best');
grid on;

figure(2);
plot(t, I_abs1, 'b-', 'LineWidth', lw)
hold on;
plot(t, I_abs2, 'r-', 'LineWidth', lw);
plot(t, I_abs3, 'g-', 'LineWidth', lw);
plot(t, I_abs4, 'm-', 'LineWidth', lw);

set(gca, 'FontSize', 24, 'LineWidth', lw);
xlim([0, 1000]);
ylim([-1, 20]);
xlabel("Vreme [min]");
ylabel("Insulin I(t) [mU/L]");
title("Kontrolni signal (insulin)");
legend("I_0 = 16", "I_0 = 15", "I_0 = 14", "I_0 = 0", 'location', 'best');
grid on;

figure(3);
plot(t, e_dev1, 'b-', 'LineWidth', lw)
hold on;
plot(t, e_dev2, 'r-', 'LineWidth', lw);
plot(t, e_dev3, 'g-', 'LineWidth', lw);
plot(t, e_dev4, 'm-', 'LineWidth', lw);
plot(t, zeros(size(t)), 'k-.', 'LineWidth', 1);

set(gca, 'FontSize', 24, 'LineWidth', lw);
xlim([0, 1000]);
ylim([-8, 4]);
xlabel("Vreme [min]");
ylabel("e(t)");
title("Greska (r(t)-y_m(t))");
legend("I_0 = 16", "I_0 = 15", "I_0 = 14", "I_0 = 0", "0", 'location', 'best');
grid on;


