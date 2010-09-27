% CPTEMDFunction(y, Ts, DesiredArcLength, f_max)

clear
close all
clc

% define the test signal
% t_max = 1;
% Fs = 5e3;
% Ts = 1/Fs;
% NSamples = t_max*Fs;
% t = linspace(0,t_max,NSamples);
% f = 30;
% y = cos(2*pi*f*t);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Fs = 1e3;                               % samples/second
% Ts = 1/Fs;                              % sample period (radians)
% Duration = 3;                         % seconds
% NSamples = Duration*Fs;
% t{1} = linspace(0,Duration,NSamples);
% 
% f1 = 7;                    Theta1 = 2*pi*f1*t{1};                % frequency Hz
% f2 = 13;                  Theta2 = 2*pi*f2*t{1};
% f_max = 1*f2;
% 
% alpha = 1;
% A1 = 1/(f1^alpha);              % power falls off at 1/(f^2) and amplitude falls away at 1/f
% A2 = 1/(f2^alpha);   
% 
% x1 = A1*cos(Theta1);
% x2 = A2*cos(Theta2);
% 
% y = x1 + x2;

% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~

% load CSGridDepth200608CS33_19900_Chs_1_48_200s

load HGSP7_98s_seizure
Fs = 4069.010498046875;     % Hz
Ts = 1/Fs;

% pick channel
x = detrend(Data(1:end-round(Fs*2),end) - Data(1:end-round(2*Fs),end-1));
% x = detrend(Data(:,end-1));
NSamples = length(x);
t{1} = 0:Ts:(length(x)-1)/Fs;        % seconds

clear data

% use a median filter to give the data a first clean
x = medfilt1(x,20);

% low pass filter file at 40 Hz
Fc = [2.5 95];                                % Hz
f_max = 200;
Wc = Fc/(Fs/2);                     % normalised digital frequency
[b a] = butter(3,Wc);
y = filtfilt(b,a,x)';

% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~



%parameters for the cpt
psi = 3*pi/2;

% initialise emd
foundarc = true;
res{1} = y;
m=1;
InitialArcSamples = 10;                                   % Number of samples in the first try to find the 'DesiredArcLength'
init_b_f = InitialArcSamples*ones(1,NSamples);
zeta = 1;
% run emd
while foundarc
    
%     figure
%     plot(init_b_f)
%     drawnow
    
    [foundarc phi{m} r{m} t{m} firstindex lastindex  ArcPoints TangentPoints] = CPTfunction(res{m}, t{m}, Ts, psi, f_max, init_b_f,zeta);
    zeta = zeta+2;
    temp = res{m};
    res{m} = temp(firstindex:lastindex);            % now res is the same length as t
    
    if foundarc
        
        % find the detrending line
        cos_phi = cos(phi{m});
        p1 = findpeaks(cos_phi);        % indexes for maxima
        p2 = findpeaks(-cos_phi);       % indexes for minima
    %     p3 = sort([p1 p2]);
        
        if (length(p1) < 2) || (length(p2) < 2)
            foundarc  = false;
        else
            if p1(1) > p2(1);
                min_index = p1(1);
            else
                min_index = p2(1);
            end

            if p1(end) > p2(end)
                max_index = p2(end);
            else
                max_index = p1(end);
            end
            spline_index = min_index:max_index;

            if length(spline_index) < 2
                foundarc  = false;
            else
                
                
                maxima = res{m};
                maxima = maxima(p1);
                minima = res{m};
                minima = minima(p2);

                s1 = spline(p1,maxima,1:length(res{m}));
                s2 = spline(p2,minima,1:length(res{m}));

                m=m+1;

                res{m} = (s1(spline_index)+s2(spline_index))/2;        % this is the residue or the detrending line
                NSamples = length(res{m});
                t_temp = t{m-1};
                t{m} = t_temp(spline_index);                    % t is shrunk according to the spline indexes
                
%                 temp = res{m-1};
%                 temp = temp(spline_index);
%                 r_hat{m-1} = abs(temp-res{m});
%                 temp = phi{m-1};
%                 temp = temp(spline_index);
%                 phi_hat{m-1} = temp;
%                 temp = t{m-1};
%                 t_hat{m-1} = temp(spline_index);
                
                % here we cut down the edges of the arcs the define the bit
                % we fit the circle to and the tangent
                temp = [ArcPoints(:,spline_index); TangentPoints(:,spline_index)];
                init_b_f  = min(temp,[],1);
                
                % ~~~~~~~~~~~
                figure
                plot(res{m-1})
                hold on
                plot(p1,maxima,'g.')
                plot(p2,minima,'r.')
                plot(spline_index,s1(spline_index),'k')
                plot(spline_index,s2(spline_index),'k')        
                plot(spline_index,res{m},'c')
                hold off
                drawnow
                %  ~~~~~~~~~~~~~~
                
%                 figure
%                 plot(phi{m-1})
%                 ylabel(['\phi_' num2str(m) '(n)'])
                % ~~~~~~~~~~~~~~~
            end
        end
    end
end

% subplotmax = 3;
% for n=1:subplotmax
%     
%     subplot(subplotmax,1,n)
%     plot(t_hat{n},r_hat{n}.*cos(phi_hat{n})),hold on
%     plot(t{n},res{n},'r')
%     hold off
% end

figure
subplotmax = 3;
for n=1:subplotmax
    
    subplot(subplotmax,1,n)
    plot(t{n},r{n}.*cos(phi{n})),hold on
    plot(t{n},res{n},'r')
    hold off
end