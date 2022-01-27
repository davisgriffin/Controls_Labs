%% ECE 5390 - Practicum 1 - Least Squares Fitting
%  G.Davis
%  01/25/2022

clc; clear; close all;

%% Data Setup

global y
global t

lab_data = [0 0
 0.5000 0.3268
 1.0000 0.5913
 1.5000 0.7521
 2.0000 0.8496
 2.5000 0.9088
 3.0000 0.9447
 3.5000 0.9664
 4.0000 0.9796
 4.5000 0.9877
 5.0000 0.9925
 5.5000 0.9955
 6.0000 0.9972
 6.5000 0.9983
 7.0000 0.9990
 7.5000 0.9994
 8.0000 0.9996
 8.5000 0.9998
 9.0000 0.9999
 9.5000 0.9999
 10.0000 0.9999];

t = lab_data(:,1);
y = lab_data(:,2);

%% Practicum Examples
% TF = 1/(s^2 + alpha*s + 2)
% alpha = 2;
% sys = tf([0 0 1], [1 alpha 2]);
% [y,t] = step(sys,10); % get the step response of the system for 10 seconds
% figure(1)
% plot(t,y); grid
% 
% % create reference vector to compute errors
% % play with blue curve (y) til it matches red curve (yref)
% yref = 1.4*y;
% figure(2)
% plot(t,y,t,yref); grid
% 
% errordiff = y-yref;
% LSE = (errordiff' * errordiff);
% 
% for i=1:20
%     
% end

%% First Order System

% educate our guesses based on the data.
[~,idx] = min(abs(y-0.632));
Tau = y(idx);

alpha = 1/Tau;
k = alpha*y(4*idx);

% create initial reference
LSE = getLSE(alpha, k);
prevLSE = LSE;

%-------------------
alpha_pos = alpha; alpha_neg = alpha;
aLSE_pos = LSE; aLSE_neg = LSE;
aMag_pos = 0.2; aMag_neg = 0.2;

k_pos = k; k_neg = k;
kLSE_pos = LSE; kLSE_neg = LSE;
kMag_pos = 0.2; kMag_neg = 0.2;

mag_order = 0.001;
%-------------------

for i=1:1000
    if LSE == 0
        break
    end
    
    if mod(i,2)
        % work on alpha       
        alpha_pos = alpha + aMag_pos;
        alpha_neg = alpha - aMag_neg;
        
        aLSE_pos = getLSE(alpha_pos, k);
        aLSE_neg = getLSE(alpha_neg, k);
        
        if aLSE_pos > prevLSE && aLSE_neg > prevLSE
            aMag_pos = aMag_pos - mag_order;
            aMag_neg = aMag_neg - mag_order;
        else
            if aLSE_pos > aLSE_neg
                aMag_pos = aMag_pos - mag_order;
                alpha = alpha_neg;
            else
                aMag_neg = aMag_neg - mag_order;
                alpha = alpha_pos;
            end
        end
    else
        % work on k
        k_pos = k + kMag_pos;
        k_neg = k - kMag_neg;
        
        kLSE_pos = getLSE(alpha, k_pos);
        kLSE_neg = getLSE(alpha, k_neg);
        
        if kLSE_pos > prevLSE && kLSE_neg > prevLSE
            kMag_pos = kMag_pos - mag_order;
            kMag_neg = kMag_neg - mag_order;
        else
            if kLSE_pos > kLSE_neg
                kMag_pos = kMag_pos - mag_order;
                k = k_neg;
            else
                kMag_neg = kMag_neg - mag_order;
                k = k_pos;
            end
        end
    end
    prevLSE = LSE;
    LSE = getLSE(alpha, k);
end

sys = tf([0 k], [1 alpha]);
[y_m,~] = step(sys,1:length(t));
figure, plot(t,y,t,y_m), grid
xlim([min(t)-1 max(t)+1])
ylim([min(y)-std(y) max(y)+std(y)])

%% Second Order System

function LSE = getLSE(alpha, k)
    global y
    global t
    sys = tf([0 k], [1 alpha]);
    [y_m,~] = step(sys,1:length(t));
    errordiff = y - y_m;
    LSE = errordiff' * errordiff;
end
