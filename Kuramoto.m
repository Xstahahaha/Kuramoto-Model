clear all
format long;
N = 1000; 
tspan = [0 50];
tstep = 0.05;
t = tspan(1):tstep:tspan(2);
theta_mean = zeros(N,tspan(2)/tstep+1);
R_mean = zeros(1,tspan(2)/tstep+1);
loops = 1;
for p = 1:loops
%omega = randn(N,1); 
omega = normrnd(0.3,0.05,N,1); 
K = 0.291920622; % Acoplamiento global

theta0 = rand(N,1) * 2 * pi;


theta = zeros(N,length(t));
theta(:,1) = theta0;

for i = 2:length(t)
    k1 = kuramoto_ode(t(i-1), theta(:,i-1), omega, K);
    k2 = kuramoto_ode(t(i-1) + tstep/2, theta(:,i-1) + tstep/2 * k1, omega, K);
    k3 = kuramoto_ode(t(i-1) + tstep/2, theta(:,i-1) + tstep/2 * k2, omega, K);
    k4 = kuramoto_ode(t(i-1) + tstep, theta(:,i-1) + tstep * k3, omega, K);
    
    theta(:,i) = theta(:,i-1) + (tstep/6) * (k1 + 2*k2 + 2*k3 + k4);
end
theta_mean = theta_mean + theta; 
R = abs(mean(exp(1i * theta), 1));
R2 = sum(exp(1i * theta), 1)/N;
R_mean = R_mean + R ; 
end
theta_mean = theta_mean/loops; 
R_mean = R_mean/loops; 

omega_normalized = (omega-min(omega))/(max(omega)-min(omega));
cmap = jet(256); % Obtener el mapa de colores 'jet'
omega_colors = interp1(linspace(0, 1, 256), cmap, omega_normalized); % Mapear los valores normalizados al mapa de colores 'jet'
figure(4), clf(4);
tiledlayout(1,2);
nexttile;
polarscatter(theta(:,1), 1,50,omega_colors,'filled'); % Trazar el primer conjunto de puntos
hold on;
polarplot(angle(R2(1,1)), abs(R2(1,1)), 'x'); % Trazar el segundo conjunto de puntos
hold off;
title('Configuración Inicial');
  rticks(0:0.2:1); % Sin tics en el radio
thetaticks(0:60:360); % Sin tics en el ángulo
thetaticklabels({'0', '\pi/3', '2\pi/3', '\pi', '4\pi3', '5\pi/3', '2\pi'}); % Etiquetas de los ángulos en radianes
% Segundo mosaico
nexttile;
polarscatter(theta(:,end), 1,50,omega_colors,'filled'); % Trazar el primer conjunto de puntos
hold on;
polarplot(angle(R2(1,end)), abs(R2(1,end)), 'x'); % Trazar el segundo conjunto de puntos
hold off;
title('Tiempo Final');
  rticks(0:0.2:1); % Sin tics en el radio
thetaticks(0:60:360); % Sin tics en el ángulo
thetaticklabels({'0', '\pi/3', '2\pi/3', '\pi', '4\pi3', '5\pi/3', '2\pi'}); % Etiquetas de los ángulos en radianes

function dtheta_dt = kuramoto_ode(t, theta, omega, K)
    N = length(theta);
    dtheta_dt = zeros(N,1);
    for i = 1:N
        sum_term = 0;
        for j = 1:N
            sum_term = sum_term + sin(theta(j) - theta(i));
        end
        dtheta_dt(i) = omega(i) + (K/size(theta,1)) * sum_term;
    end
end
