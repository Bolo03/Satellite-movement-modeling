clc; clear; close all
load constants.mat
%% Cerinta 1
epsilon0 = 1e4;
k = 5;
%% Cerinta 2
% realizata in Simulink
%% Cerinta 3
g = @(t, r, dr) -((G*Mp)/norm(r)^2) .* (r/norm(r)) - ...
       (3/2 * J2 * G * Mp * Rp^2 / norm(r)^5) .* ...
       [(r(1) - 5 * r(1) * r(3)^2)/norm(r)^2;
        (r(2) - 5 * r(2) * r(3)^2)/norm(r)^2;
        (r(3) - 5 * r(3) * r(3)^2)/norm(r)^2] + ...
        omegap^2 .* [r(1); r(2); 0] + ...
        2*omegap .* [dr(2); -dr(1); 0];

f = @(t, x) [x(4:6); g(t, x(1:3), x(4:6)); x(10:12); g(t, x(7:9), x(10:12)); norm(x(1:3) - x(7:9))/norm(x(7:9))];
%% Cerinta 4
% pas de integrare si Tmax
Tmax = 1800;
h = 1;

% variabile de start
x = [r1_0; dr1_0; r2_0; dr2_0; epsilon0];
t = 0;

num_steps = ceil(Tmax / h);
X = zeros(size(x, 1), num_steps);
T = zeros(1, num_steps);

% Metoda Runge-Kutta
i = 1;
while t < Tmax
    k1 = f(t, x);
    k2 = f(t + h/2, x + h*k1/2);
    k3 = f(t + h/2, x + h*k2/2);
    k4 = f(t + h, x + h*k3);

    X(:, i) = x;
    T(i) = t;

    x = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    t = t + h;

    i = i + 1;
end

r1 = X(1:3, :);
r2 = X(7:9, :);

figure;
hold on; grid on;
plot(r1(1, :), "LineWidth", 1);
plot(r1(2, :), "LineWidth", 1);
plot(r1(3, :), "LineWidth", 1);
xlabel("$t$",  "Interpreter","latex");
ylabel("$r_1$", "Interpreter","latex");
legend("$r_{1x}$", "$r_{1y}$", "$r_{1z}$", "Interpreter","latex");
title("Evolutia $r_{1}$ pe axele de coordonate", "Interpreter", "latex", "FontSize", 12);

figure
hold on; grid on;
plot(r2(1, :), "LineWidth", 1);
plot(r2(2, :), "LineWidth", 1);
plot(r2(3, :), "LineWidth", 1);
xlabel("$t$",  "Interpreter","latex");
ylabel("$r_2$", "Interpreter","latex");
legend("$r_{2x}$", "$r_{2y}$", "$r_{2z}$", "Interpreter","latex");
title("Evoluia $r_{2}$ pe axele de coordonate", "Interpreter", "latex", "FontSize", 12);
%% Cerinta 5
% Runge-Kutta
figure;
view(3);
hold on; grid on;
plot3(r1(1, :), r1(2, :), r1(3, :), "b", "LineWidth", 1);
plot3(r2(1, :), r2(2, :), r2(3, :), "k", "LineWidth", 1);
legend("$r_1$", "$r_2$", "Interpreter","latex", "Location", "northeast");
title("Orbita satelitilor: Runge-Kutta", "Interpreter", "latex", "FontSize", 12);
xlabel("x");
ylabel("y");
zlabel("z");

% Simulink
mdl = 'satelit_mdl';
% intrarea
t = 0:h:Tmax;
u = k * 1e-3 * double(t>=0);
in_u = timeseries(u, t);

% iesirea
load_system(mdl);
set_param(mdl, "StopTime", num2str(Tmax));
out = sim(mdl);

slx_r1 = squeeze(out.r1.Data);
slx_r2 = squeeze(out.r2.Data);

figure;
view(3);
hold on; grid on;
plot3(slx_r1(1, :), slx_r1(2, :), slx_r1(3, :), "b", "LineWidth", 1);
plot3(slx_r2(1, :), slx_r2(2, :), slx_r2(3, :), "k", "LineWidth", 1);
legend("$r_1$", "$r_2$", "Interpreter","latex", "Location", "northeast");
title("Orbita satelitilor: Simulink", "Interpreter", "latex", "FontSize", 12);
xlabel("x");
ylabel("y");
zlabel("z");


% Interpolare rezultat Simulink, pentru un plot mai fin
% Acest lucru este vizibil la Tmax mare
num_points = size(slx_r1, 2);
t_original = linspace(0, 1, num_points);
t_fine = linspace(0, 1, num_points * 10);

slx_r1_smooth = interp1(t_original, slx_r1', t_fine, 'spline')';
slx_r2_smooth = interp1(t_original, slx_r2', t_fine, 'spline')';

figure;
view(3);
hold on; grid on;
plot3(slx_r1_smooth(1, :), slx_r1_smooth(2, :), slx_r1_smooth(3, :), "b", "LineWidth", 1);
plot3(slx_r2_smooth(1, :), slx_r2_smooth(2, :), slx_r2_smooth(3, :), "k", "LineWidth", 1);
legend("$r_1$", "$r_2$", "Interpreter","latex", "Location", "northeast");
title("Orbita satelitilor: Simulink dupa interpolare", "Interpreter", "latex", "FontSize", 12);
xlabel("x");
ylabel("y");
zlabel("z");

%% Cerinta 6
y_rk = X(13, :);

y_slx = squeeze(out.y.Data);

% interpolare
num_points = size(y_rk, 2);
t_original = linspace(0, 1, size(y_slx, 1));
t_fine = linspace(0, 1, num_points);
y_slx_smooth =  interp1(t_original, y_slx', t_fine, 'spline');

figure;
hold on;
plot(y_rk, "LineWidth",2);
plot(y_slx_smooth, "LineWidth",2);
xlabel("$t$","Interpreter","latex");
ylabel("$y$", "Interpreter","latex");
legend("$y^{RK}$", "$y^{Slx}$", "Interpreter","latex");
title("Iesirea sistemului", "Interpreter", "latex", "FontSize", 12);

e = vecnorm(y_rk - y_slx_smooth, 2, 1);

figure;
plot(e, 'r-',"LineWidth",2);
xlabel("$t$","Interpreter","latex");
ylabel("$\epsilon$", "Interpreter","latex");
title("Eroarea de integrare", "Interpreter", "latex", "FontSize", 12);

%% Cerinta 7
t = 0:h:Tmax;
k_star = [5, 10, 50, 100, 250, 500]; % pentru u_star
epsilon_star = zeros(1, numel(k_star));
u_star = zeros(1, numel(k_star));

for i = 1:numel(k_star)
    u = k_star(i) * 1e-3 * double(t>=0);
    in_u = timeseries(u, t);

    out = sim(mdl);
    epsilon_out = squeeze(out.y.Data);

    epsilon_star(i) = epsilon_out(end);
    u_star(i) = u(end);
end

figure;
plot(u_star, epsilon_star, 'r*', "LineWidth", 2);
xlabel("$u^*$", "Interpreter", "latex");
ylabel("$\varepsilon^*$", "Interpreter", "latex");
title("Caracteristica $\varepsilon^*(u^*)$", "Interpreter", "latex", "FontSize", 12);

%% Cerinta 8
p = polyfit(u_star, epsilon_star, 1);

epsilon_aprox = polyval(p, u_star);

figure;
hold on;
plot(u_star, epsilon_star, 'r*', "LineWidth", 2);
plot(u_star, epsilon_aprox, 'b--', "LineWidth", 1);
xlabel("$u^*$", "Interpreter", "latex");
ylabel("$\varepsilon^*$", "Interpreter", "latex");
title("Caracteristica $\varepsilon^*(u^*)$", "Interpreter", "latex", "FontSize", 12);
%% Cerinta 9
load valori.mat % pentru a reseta valoarea lui r2_0 inainte de simulari
% folosim intrarea aleasa la cerinta 1
t = 0:h:Tmax;
u = k * 1e-3 * double(t>=0);
in_u = timeseries(u, t);

% N(0, 0.1)
mu = 0; % media
var = 0.1; % varianta

r2_0_copy = r2_0;

nr_iter = 100;
sample_size = Tmax;
y = zeros(nr_iter, sample_size);

figure;
hold on;
for i = 1:nr_iter
    alpha = sqrt(var) .* randn(1, 1) + mu;
    r2_0 = (1 + alpha) .* r2_0_copy;

    out = sim(mdl);

    y_slx = squeeze(out.y.Data)';

    % interpolare pentru a corespunde cu suportul de timp
    t_original = linspace(0, 1, size(y_slx, 2));
    t_fine = linspace(0, 1, sample_size);
    y_slx_smooth =  interp1(t_original, y_slx', t_fine, 'spline');

    y(i, :) = y_slx_smooth;

    plot(y(i, :), ':', 'HandleVisibility','off');    
end

y_mean = mean(y, 1);
plot(y_mean, 'r', 'LineWidth', 2);
legend("$\bar{y}$", "Interpreter","latex");
xlabel("$t$","Interpreter","latex");
ylabel("$\varepsilon$", "Interpreter","latex");
title("Incetitudine multiplicativa", "Interpreter", "latex", "FontSize", 12);

% Analiza probabilitatii distantei relative dintre sateliti
% pragul este definit ca ultima valoare din y_mean
epsilon_prag = y_mean(end);

% cautam simularile care depasesc pragul
depasesc_prag = max(y, [], 2) > epsilon_prag;

probabilitate = sum(depasesc_prag) / nr_iter;

fprintf('Probabilitatea ca ε(t) să depășească y_{mean}(end) = %.2f este %.2f%%\n', ...
        epsilon_prag, probabilitate * 100);

% Reprezentare grafica a simularilor care depasesc pragul
figure;
hold on;
for i = 1:nr_iter
    if depasesc_prag(i)
        plot(y(i, :), 'r:', 'HandleVisibility','off'); % simulari care depasesc pragul
    else
        plot(y(i, :), 'b:', 'HandleVisibility','off'); % simulari care nu depasesc pragul
    end
end

plot(y_mean, 'r-', 'LineWidth', 2, 'DisplayName', '$\overline{y}$');
legend('Interpreter', 'latex', 'FontSize', 12);
xlabel("$t$", "Interpreter", "latex");
ylabel("$\varepsilon$", "Interpreter", "latex");
title("Analiza depasirii pragului $\bar{y}(\mathrm{end})$", "Interpreter", "latex", "FontSize", 12);
hold off;

%% Cerinta 10
load valori.mat % pentru a reseta valoarea lui r2_0 inainte de simulari
% folosim intrarea aleasa la cerinta 1
t = 0:h:Tmax;
u = k * 1e-3 * double(t>=0);
in_u = timeseries(u, t);

% N(0, 5)
mu = 0; % media
var = 5; % varianta

r2_0_copy = r2_0;

nr_iter = 100;
sample_size = Tmax;
y = zeros(nr_iter, sample_size);

figure;
hold on;
for i = 1:nr_iter
    alpha = sqrt(var) .* randn(numel(r2_0), 1) + mu;
    r2_0 = alpha + r2_0_copy;

    out = sim(mdl);

    y_slx = squeeze(out.y.Data)';

    % interpolate pentru outputuri mai mici decat sample_size
    num_points = size(y_slx, 2);
    t_original = linspace(0, 1, size(y_slx, 2));
    t_fine = linspace(0, 1, sample_size);
    y_slx_smooth =  interp1(t_original, y_slx', t_fine, 'spline');

    y(i, :) = y_slx_smooth;

    plot(y(i, :), ':', 'HandleVisibility','off');    
end

y_mean = mean(y, 1);
plot(y_mean, 'r', 'LineWidth', 2);
xlim([1008.786 1008.792]);
legend("$\bar{y}$", "Interpreter","latex");
xlabel("$t$","Interpreter","latex");
ylabel("$\varepsilon$", "Interpreter","latex");
title("Incetitudine aditiva", "Interpreter", "latex", "FontSize", 12);

% Analiza probabilitatii distantei relative dintre sateliti
% pragul este definit ca ultima valoare din y_mean
epsilon_prag = y_mean(end);

% cautam simularile care depasesc pragul
depasesc_prag = max(y, [], 2) > epsilon_prag;

probabilitate = sum(depasesc_prag) / nr_iter;

fprintf('Probabilitatea ca ε(t) să depășească y_{mean}(end) = %.2f este %.2f%%\n', ...
        epsilon_prag, probabilitate * 100);

% Reprezentare grafica a simularilor care depasesc pragul
figure;
hold on;
for i = 1:nr_iter
    if depasesc_prag(i)
        plot(y(i, :), 'r:', 'HandleVisibility','off'); % simulari care depasesc pragul
    else
        plot(y(i, :), 'b:', 'HandleVisibility','off'); % simulari care nu depasesc pragul
    end
end

plot(y_mean, 'r-', 'LineWidth', 2, 'DisplayName', '$\overline{y}$');
xlim([1008.782 1008.795]);
legend('Interpreter', 'latex', 'FontSize', 12);
xlabel("$t$", "Interpreter", "latex");
ylabel("$\varepsilon$", "Interpreter", "latex");
title("Analiza depasirii pragului $\bar{y}(\mathrm{end})$", "Interpreter", "latex", "FontSize", 12);
hold off;
%% Cerinta 11
t = 0:1:Tmax; t = t';
sample_size = Tmax;

u = randn(numel(t), 1); % apartine N(0,1)
in_u = timeseries(u, t);

out = sim(mdl);

y_slx = squeeze(out.y.Data)';

% interpolare pentru outputuri mai mici decat sample_size
num_points = size(y_slx, 2);
t_original = linspace(0, 1, size(y_slx, 2));
t_fine = linspace(0, 1, sample_size);
y_slx_smooth =  interp1(t_original, y_slx', t_fine, 'spline');

second_dy = diff(y_slx_smooth, 2);

figure;
hold on;
plot(t, u, 'b');
plot(second_dy, 'r', "LineWidth", 2);
ylim([-0.005 0.005]);
xlabel("$t$", "Interpreter", "latex");
legend("$u$", "$y$", "Interpreter", "latex");
title("Semnalele de intrare si iesire", "Interpreter", "latex", "FontSize", 12);

second_dy_mean = mean(second_dy);
second_dy_std = std(second_dy);