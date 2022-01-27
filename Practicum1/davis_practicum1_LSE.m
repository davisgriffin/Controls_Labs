%% ECE 5390 - Practicum 1 - Least Squares Fitting
%  G.Davis
%  01/25/2022

clc; clear; close all;

%% Data Setup

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
magn = 0.1; % magnitude by which to change
operator = 0; % 0 add 1 subtract
state = 0;

% create initial reference
LSE = getLSE(alpha, k, y, t);
prevLSE = ones(1,3); % two-bit memory of error
tempLSE = 1

for i=1:100

    if LSE == 0
        break
    end
     % add logic that if it's only greater than the most previous one to
     % step back and decrease magnitude.
     % if greater than both reset magnitude (*2) and restore with new
     % operator?
    while state ~= 2
         if LSE > tempLSE
            alpha = alterParam(alpha, ~operator, magn);
            magn = magn/2;
            alpha = alterParam(alpha, operator, magn);
            tempLSE = getLSE(alpha, k, y, t);
            state = state + 1;
         else
             break
         end
    end
    if LSE > tempLSE && state == 3
        
    end
    state = 0;
     
%     if LSE > prevLSE(3)
%         alpha = alterParam(alpha, ~operator, magn);
%         magn = magn/2;
%     end
    if LSE > prevLSE(2) && LSE > prevLSE(3)
        operator = ~operator;
%         magn = magn*2;
        % reset to previous
        alpha = alterParam(alpha, operator, magn*2);
        alpha = alterParam(alpha, operator, magn*2);
%         magn = magn/2; % likely need new logic for this line
    end
    
    alpha = alterParam(alpha,operator,magn);
    
    prevLSE = [prevLSE(2) prevLSE(3) LSE];
    tempLSE = LSE;
    LSE = getLSE(alpha, k, y, t);
end

figure, plot(t,y,t,y_m), grid
xlim([min(t)-1 max(t)+1])
ylim([min(y)-std(y) max(y)+std(y)])

function LSE = getLSE(alpha, k, y, t)
    sys = tf([0 k], [1 alpha]);
    [y_m,~] = step(sys,1:length(t));
    errordiff = y - y_m;
    LSE = errordiff' * errordiff;
end

function a = alterParam(x, operator, magn)
    if operator
        a = x - magn;
    else
        a = x + magn;
    end
end
