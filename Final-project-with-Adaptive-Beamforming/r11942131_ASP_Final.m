clear;
clc;

data = load("ASP_Final_Data.mat");
X=data.matX;
L=length(X(1,:));
N=length(X(:,1));
theta_s_noisy=data.theta_s_noisy;
theta_i_noisy=data.theta_i_noisy;
%L=length(theta_i_noisy);
t=1:L;

delta = 0.01;
lambda = 0.5;

theta_s_hat = Denoise_theta(delta, lambda, L, theta_s_noisy);
theta_i_hat = Denoise_theta(delta, lambda, L, theta_i_noisy);

figure();
plot(t,theta_s_hat);
hold on;
plot(t,theta_i_hat);
axis([0 2000 -10 20]);
title('DOAs $\hat{\theta}$','Interpreter','latex','fontweight','bold');
xlabel('Time index t');
ylabel('DOA (degree)','fontweight','bold');
yticks([-10 -5 0 10 20]);
yticklabels({'-10','-5','0','10','20'});
legend('DOAs $\hat{\theta}_s(t)$','DOAs $\hat{\theta}_i(t)$','Interpreter','latex');

%% s_t_hat_MVDR plot

y_MVDR_hat = MVDR_beamformer(X, N, L, theta_s_hat);
    
figure();
subplot(2,1,1)
plot(t,real(y_MVDR_hat));
title('Estimated source signal $\hat{s}(t)$ through MVDR beamformer','Interpreter','latex','fontweight','bold');
xlabel('Time index t');
ylabel('real part of $\hat{s}(t)$','Interpreter','latex','fontweight','bold');
grid on

subplot(2,1,2)
plot(t,imag(y_MVDR_hat));
title('Estimated source signal $\hat{s}(t)$ through MVDR beamformer','Interpreter','latex','fontweight','bold');
xlabel('Time index t');
ylabel('imag part of $\hat{s}(t)$','Interpreter','latex','fontweight','bold');
grid on;

%% s_t_hat_LCMV plot

y_LCMV_hat = LCMV_beamformer(X, N, L, theta_s_hat, theta_i_hat);
    
figure();
subplot(2,1,1)
plot(t,real(y_LCMV_hat));
title('Estimated source signal $\hat{s}(t)$ through LCMV beamformer','Interpreter','latex','fontweight','bold');
xlabel('Time index t');
ylabel('real part of $\hat{s}(t)$','Interpreter','latex','fontweight','bold');
grid on

subplot(2,1,2)
plot(t,imag(y_LCMV_hat));
title('Estimated source signal $\hat{s}(t)$ through LCMV beamformer','Interpreter','latex','fontweight','bold');
xlabel('Time index t');
ylabel('imag part of $\hat{s}(t)$','Interpreter','latex','fontweight','bold');
grid on

%% s_t_hat_Designed plot

s_t_hat = Designed_beamformer(X, N, L, theta_s_hat, theta_i_hat, delta, lambda);

figure();
subplot(2,1,1)
plot(t,real(s_t_hat));
title('Estimated source signal $\hat{s}(t)$ through my designed beamformer','Interpreter','latex','fontweight','bold');
xlabel('Time index t');
ylabel('real part of $\hat{s}(t)$','Interpreter','latex','fontweight','bold');
grid on

subplot(2,1,2)
plot(t,imag(s_t_hat));
title('Estimated source signal $\hat{s}(t)$ through my designed beamformer','Interpreter','latex','fontweight','bold');
xlabel('Time index t');
ylabel('imag part of $\hat{s}(t)$','Interpreter','latex','fontweight','bold');
grid on

%% Output

save('r11942131_馬振洋_ASPFinal_PerformanceEvaluation.mat','theta_s_hat','theta_i_hat','s_t_hat');

%% Denoise algorithm
function theta_hat = Denoise_theta(delta, lambda, L, theta_noisy)

    P = 1/delta*ones(L);
    k=[];
    q=linspace(0,2000,2000);
    t_peak = [];
    theta_noisy_peak = [];
    t_dip = [];
    theta_noisy_dip = [];
    for t=1:L
        k(t) = ((1/lambda)*P(t)) / (1+(1/lambda)*P(t));
        if t==1
            t_peak = [t_peak q(t)];
            theta_noisy_peak = [theta_noisy_peak theta_noisy(t)];
            t_dip = [t_dip q(t)];
            theta_noisy_dip = [theta_noisy_dip theta_noisy(t)];
        end
        if t~=1 && t~=L
            if (theta_noisy(t)>theta_noisy(t-1)) && (theta_noisy(t)>theta_noisy(t+1))
                t_peak = [t_peak q(t)];
                theta_noisy_peak = [theta_noisy_peak theta_noisy(t)/(k(t)^2)];
            elseif (theta_noisy(t)<theta_noisy(t-1)) && (theta_noisy(t)<theta_noisy(t+1))
                t_dip = [t_dip q(t)];
                theta_noisy_dip = [theta_noisy_dip theta_noisy(t)*(k(t)^2)];
            end
        end
        if t==L
            t_peak = [t_peak q(t)];
            theta_noisy_peak = [theta_noisy_peak theta_noisy(t)];
            t_dip = [t_dip q(t)];
            theta_noisy_dip = [theta_noisy_dip theta_noisy(t)];
            break;
        end
        P(t+1)=((1/lambda)*P(t)-(1/lambda)*k(t)*P(t))/2 + P(t)*k(t);
    end
    theta_noisy_1_peak = spline(t_peak,theta_noisy_peak,q);      %connect local peaks
    theta_noisy_1_dip = spline(t_dip,theta_noisy_dip,q);         %connect local dips
    theta_hat = (theta_noisy_1_peak + theta_noisy_1_dip)/2;      %compute mean

end


%% MVDR beamformer
function y_MVDR_hat = MVDR_beamformer(X, K, M, theta_s_hat)

    
    R_hat=zeros(K);
    for m=1:M
        R_hat=R_hat+(1/K)*((X(:,m))*(X(:,m)'));
    end

    a_theta_s=zeros(K,M);
    for m=1:M
        for k=1:K
        a_theta_s(k,m)= exp(1j*(k-1)*pi*sin(theta_s_hat(m)*pi/180));
        end
    end
    
    w_MVDR=zeros(K,M);
    for m=1:M
        w_MVDR(:,m)=(R_hat\a_theta_s(:,m))/((a_theta_s(:,m)')/R_hat*(a_theta_s(:,m)));
    end

    y_MVDR_hat=zeros(1,M);
    for m=1:M    
        y_MVDR_hat(m)=(w_MVDR(:,m)')*(X(:,m));
    end

end


%% LCMV beamformer
function y_LCMV_hat = LCMV_beamformer(X, K, M, theta_s_hat, theta_i_hat)

    
    g=[1;10^-8];
    R_hat=zeros(K);
    for m=1:M
        R_hat=R_hat+(1/K)*((X(:,m))*(X(:,m)'));
    end

    a_theta_s=zeros(K,M);
    a_theta_i=zeros(K,M);
    for m=1:M
        for k=1:K
        a_theta_s(k,m)= exp(1j*(k-1)*pi*sin(theta_s_hat(m)*pi/180));
        a_theta_i(k,m)= exp(1j*(k-1)*pi*sin(theta_i_hat(m)*pi/180));
        end
        C=[a_theta_s(:,m) a_theta_i(:,m)];
    end
    
    w_LCMV=zeros(K,M);
    for m=1:M
        w_LCMV(:,m)=(R_hat\C)/((C')/R_hat*(C))*g;
    end

    y_LCMV_hat=zeros(1,M);
    for m=1:M    
        y_LCMV_hat(m)=(w_LCMV(:,m)')*(X(:,m));
    end

end

%% Designed beamformer
function s_t_hat = Designed_beamformer(X, N, L, theta_s_hat, theta_i_hat, delta, lambda)


    g=[1;0];
    R_hat=zeros(N);
    for t=1:L
        R_hat=R_hat+(1/N)*((X(:,t))*(X(:,t)'))/(norm(X(:,t))*norm(X(:,t)'))^2;
    end
    
    R_hat=R_hat + delta*lambda^L*eye(N);

    a_theta_s=zeros(N,L);
    a_theta_i=zeros(N,L);
    for t=1:L
        for n=1:N
            a_theta_s(n,t)= exp(1j*(n-1)*pi*sin(theta_s_hat(t)*pi/180));
            a_theta_i(n,t)= exp(1j*(n-1)*pi*sin(theta_i_hat(t)*pi/180));
        end
        C=[a_theta_s(:,t) a_theta_i(:,t)];
    end

    w_designed=zeros(N,L);
    for t=1:L
        w_designed(:,t)=(R_hat\C)/((C')/R_hat*(C))*g;
    end

    s_t_hat=zeros(1,L);
    for t=1:L    
        s_t_hat(t)=(w_designed(:,t)')*(X(:,t));
    end

end



