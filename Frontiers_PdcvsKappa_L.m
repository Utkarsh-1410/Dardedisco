% Shubhi Singhal
% Modified by Shobha Sundar Ram
% 1/12/2021
% Range Limited Cell 
% Frontiers Journal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pd for the scenario when Kappa is varying
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BISTATIC RADAR PROBLEM SPACE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Area working on extends from -X_max_in_m to X_max_in_m in x-axis 
X_max_in_m = 100;
% Area working on extends from -Y_max_in_m to Y_max_in_m in y-axis 
Y_max_in_m = 100;
% No. of iterations
i_max = 10000;
% threshold SCNR 
threshold = 1;
% kappa (in m)
Kappas = linspace(10,50,20);
% Kappas = 50;
sizeKappas = size(Kappas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BISTATIC RADAR PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline length
L = [0 5 10]; 
Tx = -L/2;
Rx = +L/2;
Radar_posx=[Tx,Rx];
Radar_posy=[0,0];
% scatter(Radar_posx,Radar_posy,'ko');hold on;
% set(gca,'Xlim',[-100 100],'Ylim',[-100 100]);
% transmitted power (in Watts)
P_tx = 5;
% Aperture Efficiency (no units)
A_0 = 1;
% Fractional time 
epsilon = 0.1;

% Total time (seconds)
T = 1e-3;
% Beam time
Tbeam = 5e-6;
% B0 constant
B0 = (T/(2*pi*Tbeam));
% Beamwidth (radians)
Delta_theta = 0*epsilon;
for q = 1:length(epsilon)
    Delta_theta(q) = 1/(B0*epsilon); 
end
Delta_theta_deg = Delta_theta*180/pi;

% Gain constant
G0 = 1;
% transmitted frquency (in Hz) 
fc = 60e9;
% speed of light (in m/s)
c = 3*10^8;
% wavelength (m)
wavelength = c/fc;
% Path loss coefficient
H0 = wavelength^2/((4*pi)^3);
% angle with the origin
theta_min = 0;
theta_max = 2*pi;
% Temperature (in kappa_m)
Ts = 300;
% Boltzmann's constant
kB = 1.38064852*10^(-23);
% bandwidth (in Hz)
BW = 2e7;
tau = 1/BW;
% Noise (in Watts)
noise = kB*Ts*BW;
% Range resolution (m)
delr = c/(2*BW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLUTTER PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average clutter RCS (in m^2)
sigma_c_avg = 1;
% clutter density (1/m^2)
rho = 0.001;
% Number of clutter scatterers
Nc_avg = rho*200*200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TARGET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average target RCS (in m^2)
sigma_m_avg = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO SIMULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pd = zeros(length(L),length(Kappas));
Throughput = 0*Pd;

for q = 1:length(L)
    % Probability of detection
    for k = 1:sizeKappas(2)
        kappa_m = Kappas(k);
        % n is the number of iterations in which target is detected
        n = 0;
        % throughput constant
        tp_const = (2*pi*kappa_m) - ((3*pi*L(q)*L(q))/(8*kappa_m));
        delrk= delr/sqrt(1-(L(q)^2/(4*kappa_m^2)));
        tp_const = tp_const/delrk;
        tp_const = tp_const*(1-epsilon);
        
        for i = 1:i_max
            % selection of polar angle
            theta = rand*(theta_max-theta_min);
            % distance of target from the origin
            r = ((-(L(q)^2/2-L(q)^2*cos(theta))+((L(q)^2/2-L(q)^2*cos(theta))^2-4*1*(L(q)^4/16-kappa_m^4))^0.5)/(2*1))^0.5;
            % distance from target to the transmitter
            R_tt = ((r^2+L(q)^2/4)+r*L(q)*cos(theta))^0.5;
            % distance from target to the receiver
            R_rt = ((r^2+L(q)^2/4)-r*L(q)*cos(theta))^0.5;
            kappa_ch = sqrt(R_tt*R_rt);
            sigma_m = exprnd(sigma_m_avg);
            
            % received power
            S = (P_tx*H0*sigma_m*B0*epsilon)/(kappa_ch^4);
            % beta
            beta = acos((R_tt^2+R_rt^2-L(q)^2)/(2*R_tt*R_rt));
            % theta_t 
            theta_t = acos((R_tt^2+L(q)^2-R_rt^2)/(2*R_tt*L(q)));
            % theta_r - is this right
            theta_r = pi-theta_t-beta;

            x0 = r*cos(theta);
            y0 = r*sin(theta);
    %             scatter(x0,y0,'bx'); hold on;
            if y0<0
                theta_t = pi-theta_t;
            end

            % slopes of the two lines containing the resolution cell on
            % the transmitter side
            slope_t1 = tan(theta_t - Delta_theta/2);
            slope_t2 = tan(theta_t + Delta_theta/2);

            C = 0;
            Nc = poissrnd(Nc_avg);
            n1=0;

            %pointx = [Tx-100 Tx+100];
            %pointy = [slope_t1*(-100) slope_t1*(100)];
            %line(pointx,pointy,'Color','red','LineStyle','--'); hold on;
            %pointx = [Tx-100 Tx+100];
            %pointy = [slope_t2*(-100) slope_t2*(100)];
            %line(pointx,pointy,'Color','red','LineStyle','--'); hold on;

            for j = 1:Nc
                % location of the clutters
                x = -100 + 200*rand;
                y = -100 + 200*rand;
                rc = sqrt(x^2+y^2);
                % difference of the 'slope made by the line joining the
                % transmitter to the clutter' and the 'slopes of the two 
                % lines containing the resolution cell on
                % the transmitter side' 
                m1 = (y-0)/(x+L(q)/2)-slope_t1;
                m2 = (y-0)/(x+L(q)/2)-slope_t2;
                flag_beamwidth=1;
                % checking whether the clutter lies within the resolution
                % cell or not
                if m1*m2<0 
                    flag_beamwidth = 0;
                    if ((slope_t1*slope_t2)>0)
                        if (y0*y)<0
                            flag_beamwidth = 1;
                        end
                    elseif (slope_t1*slope_t2)<0
                        if ((x0-Tx(q))*(x-Tx(q)))<0
                            flag_beamwidth = 1;
                        end
                    end
                end
                flag_range = 1;
                if abs(rc-r) < delr/sqrt(1-(L(q)^2/(4*Kappas(k)^2)))
                    flag_range=0;
                end
                flag_beamwidth=flag_beamwidth+flag_range;
                if flag_beamwidth==0
                    % distance from clutter to the transmitter
                    R_tc = ((Tx(q)-x)^2+(0-y)^2)^0.5;
                    % distance from clutter to the receiver
                    R_rc = ((Rx(q)-x)^2+(0-y)^2)^0.5;
                    sigma_c = exprnd(sigma_c_avg);           
                    % received power
                    C = C + (P_tx*H0*sigma_c*B0*epsilon)/(R_tc^2*R_rc^2);          
                    n1=n1+1;
    %                     scatter(x,y,'rx');hold on;
    %                 else
    %                     scatter(x,y,'gx'); hold on;
                end
                
            end % for number of clutter scatterers
            SCNR = S/(C+noise);
            if SCNR>=threshold
                n = n + 1;
            end
        end % for i Monte Carlo iteration
        % Probability of detection
        
        Pd(q,k) = n/i_max;
        Throughput(q,k) = Pd(q,k)*tp_const;
    end % for k (kappa)
end
%     hold on
figure(1);
subplot(1,2,1);
plot(Kappas,Pd(1,:),'bv-','LineWidth',2,'MarkerSize',10);hold on;
plot(Kappas,Pd(2,:),'bo-','LineWidth',2,'MarkerSize',10);hold on;
plot(Kappas,Pd(3,:),'b*-','LineWidth',2,'MarkerSize',10);hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOCHASTIC GEOMETRY FORMULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pdc = zeros(length(L),length(Kappas));
SNR = Pdc;

for q = 1:length(L)
    
    for i = 1:length(Kappas)
        if Kappas(i)< (0.5*L(q))
            disp('Not co-site region');
        end
        kappa_m = Kappas(i);
        SNR_Nr = threshold*noise*(kappa_m^4);
        SNR_Dr = P_tx*H0*sigma_m_avg*B0*epsilon;
        SNR_const = SNR_Nr/SNR_Dr;
    
        SNR(i) = SNR_const;
        
        SCR_Nr = threshold*rho*c*tau*(kappa_m^2)*sigma_c_avg;
        SCR_Dr = B0*epsilon*(sigma_m_avg+(threshold*sigma_c_avg)*(kappa_m+sqrt(kappa_m^2-L(q)^2)));
        SCR_const = SCR_Nr/SCR_Dr;
        SCR(i) = SCR_const;
        
        Pdc(q,i) = exp(-(SNR(i))-SCR(i));
    
    end
end


plot(Kappas,Pdc(1,:),'rv-','LineWidth',2,'MarkerSize',10); hold on;
plot(Kappas,Pdc(2,:),'ro-','LineWidth',2,'MarkerSize',10); hold on;
plot(Kappas,Pdc(3,:),'r*-','LineWidth',2,'MarkerSize',10); hold on;
grid on;
set(gca,'Fontsize',16);
set(gca,'Ylim',[0 1]);
xlabel('Bistatic range $\kappa$ (m)','interpreter','latex','Fontsize',20,'FontWeight','bold');
ylabel('$P_{DC}^{Bi}$','interpreter','latex','Fontsize',20,'FontWeight','bold');
% ylabel('Bistatic radar detection coverage probability','Fontsize',14,'FontWeight','bold');
legend({'L = 0m (M.C.)','L = 5m (M.C.)','L = 10m (M.C)','L = 0m (SG)','L = 5m (SG)','L = 10m (S.G)'},'Location','southwest','NumColumns',2);

% saveas(gcf,'PdvsKappa_L','fig');
% saveas(gcf,'PdcvsKappa_L','png');

% figure(2)
subplot(1,2,2);
plot(Kappas,Throughput(1,:),'rv-','LineWidth',2,'MarkerSize',10);hold on;
plot(Kappas,Throughput(2,:),'bo-','LineWidth',2,'MarkerSize',10);hold on;
plot(Kappas,Throughput(3,:),'k*-','LineWidth',2,'MarkerSize',10);hold on;
grid on;
set(gca,'Fontsize',16);
xlabel('Bistatic range $\kappa$ (m)','interpreter','latex','Fontsize',20,'FontWeight','bold');
ylabel('Network Throughput $\Upsilon$ (bps)','interpreter','latex','Fontsize',20,'FontWeight','bold');
% ylabel('Bistatic radar detection coverage probability','Fontsize',14,'FontWeight','bold');
legend({'L = 0m','L = 5m','L = 10m'},'Location','southwest','NumColumns',2);

saveas(gcf,'ResultsvsKappa_L','fig');
saveas(gcf,'ResultsvsKappa_L','png');