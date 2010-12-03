% CPTEMDFunction(y, Ts, DesiredArcLength, f_max)

clear
close all
clc

% define the test signal
% ~~~~~~~~~~~~~
% t_max = 1;
% Fs = 5e3;
% Ts = 1/Fs;
% NSamples = t_max*Fs;
% t = linspace(0,t_max,NSamples);
% f = 30;
% y = cos(2*pi*f*t);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fs = 1e3;                               % samples/second
Ts = 1/Fs;                              % sample period (radians)
Duration = 1;                         % seconds
NSamples = Duration*Fs;
t = linspace(0,Duration,NSamples);

f1 = 7;                    Theta1 = 2*pi*f1*t;                % frequency Hz
f2 = 17;                  Theta2 = 2*pi*f2*t;
f_max = 1*f2;

alpha = 1;
A1 = 1/(f1^alpha);              % power falls off at 1/(f^2) and amplitude falls away at 1/f
A2 = 1/(f2^alpha);   

x1 = A1*cos(Theta1);
x2 = A2*cos(Theta2);

y = x1 + x2;

% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~
% 
% % load CSGridDepth200608CS33_19900_Chs_1_48_200s
% 
% load HGSP7_98s_seizure
% Fs = 4069.010498046875;     % Hz
% Ts = 1/Fs;
% 
% % pick channel
% x = detrend(Data(1:end-round(Fs*2),end) - Data(1:end-round(2*Fs),end-1));
% % x = detrend(Data(:,end-1));
% NSamples = length(x);
% t{1} = 0:Ts:(length(x)-1)/Fs;        % seconds
% 
% clear data
% 
% % use a median filter to give the data a first clean
% x = medfilt1(x,20);
% 
% % low pass filter file at 40 Hz
% Fc = [2 65];                                % Hz
% f_max = 200;
% Wc = Fc/(Fs/2);                     % normalised digital frequency
% [b a] = butter(2,Wc);
% y = filtfilt(b,a,x)';

% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~

% parameters for the cpt
psi = 2*pi/2;

% initialise emd
foundarc = true;
res{1} = y;
m=1;
InitialArcSamples = 20;                                   % Number of samples in the first try to find the 'DesiredArcLength'
init_b_f = InitialArcSamples*ones(1,NSamples);
zeta = 1;
start_index_offset = 0;
end_index_offset = 0;
% run emd
while foundarc
    
    tic 
    [x0 Hx0 foundarc phi{m} phi_unwrapped r{m} firstindex lastindex ArcPoints TangentPoints] ...
        = CPTfunction(res{m}, Ts, psi, f_max, init_b_f, zeta, start_index_offset, end_index_offset);
    toc
    zeta = 1.3*zeta;        % this changes the step size of b and f when finding the arc and tangent

    if foundarc
        
        % find the detrending line
        cos_phi = cos(phi{m});
        p1 = findpeaks(cos_phi);        % indexes for maxima
        p2 = findpeaks(-cos_phi);       % indexes for minima
        
        if (length(p1) < 2) || (length(p2) < 2)
            foundarc  = false;
        else
            if p1(1) > p2(1);           % check which extrema comes first
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

            % these values make it so the estimation
            start_index_offset = min_index;
            end_index_offset = NSamples-max_index;
            
            if length(spline_index) < 2
                foundarc  = false;
            else
                
                % find the points in the signal that correspond to the
                % maxima of x
                maxima = res{m};
                maxima = maxima(p1);
                
                % now for the minima
                minima = res{m};
                minima = minima(p2);

                s1 = spline(p1,maxima,1:length(res{m}));
                s2 = spline(p2,minima,1:length(res{m}));

                m=m+1;
                
                emperical_offset = (s1+s2)/2; 

                % now set the samples to zero that are before the extrema
                emperical_offset(1:min_index-1) = 0;
                emperical_offset(max_index+1:end) = 0;
                                
                res{m} = emperical_offset;        % this is the residue or the detrending line
                NSamples = length(res{m});

                % here we cut down the edges of the arcs the define the bit
                % we fit the circle to and the tangent
%                 temp = [ArcPoints(:,spline_index); TangentPoints(:,spline_index)];
%                 init_b_f  = min(temp,[],1);
                
                % ~~~~~~~~~~~
                figure
                plot(res{m-1})
                hold on
                plot(p1,maxima,'g.')
                plot(p2,minima,'r.')
                plot(s1,'k')
                plot(s2,'k')        
                plot(res{m},'c')
                hold off
                drawnow
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
    plot(t,r{n}.*cos(phi{n})),hold on
    plot(t,res{n},'r')
    hold off
end