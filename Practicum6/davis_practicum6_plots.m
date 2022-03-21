%% ECE 5390 - Practicum 1 - Modeling from Experimental Data
%  G.Davis
%  01/25/2022

clc; clear; close all;

%% Useful Encoding

% Arbitray V2 = 5V, V1 = 0V chosen.

theta = 0:0.01:2*pi;
V2 = 5; V1 = 0;
Vo = 2*(V2-V1)*(theta-pi)/(2*pi);
Vo(1:floor(numel(theta)/2)) = 0;

figure(1)
plot(theta,Vo, 'LineWidth', 1.5); grid on
title('Useful Degrees of Rotation for Angular Encoding')
xlabel('Degrees of Rotation (rad)')
ylabel('Voltage')
text(0.25, 2.5, "Arbitrary V2 = 5V")

minCx = floor(min(theta)/pi);
maxCx = ceil(max(theta)/pi);
ticks = piRationM(minCx:0.25:maxCx);
xticks((minCx:0.25:maxCx) * pi)
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', ticks)

%% Vo Normalized

theta = -pi/2:0.01:pi/2;
V2 = 10; Kp = (V2-V1)/pi;
Vo = Kp*(theta + pi/2) - (V2-V1)/2;

figure(2)
plot(theta,Vo, 'r', 'LineWidth', 1.5); grid on
title('$V_o$ for $V_2 = 10 V$', 'interpreter', 'latex')
xlabel('Potentiometer Wiper Position (rad)')
ylabel('Voltage')

minCx = floor(min(theta)/pi);
maxCx = ceil(max(theta)/pi);
ticks = piRationM(minCx:0.25:maxCx);
xticks((minCx:0.25:maxCx) * pi)
set(gca, 'TickLabelInterpreter', 'latex', 'XTickLabel', ticks)

%% Negligible Inertia and Dampening

Jm = 1; Dm = 1;
tfcn = tf([0 0 1],[Jm Dm 1]);
figure(3)
step(tfcn)

% Now considering spring friction
K = 1;
tfcn = tf([0 0 K], [Jm Dm K]);
figure(4)
step(tfcn)

%% Functions

function str = piRation(x)
    [num, den] = rat(x);
    if x == 0
        str = '$0$';
    elseif x == 1
        str = '$\pi$';
    elseif den == 1
        str = '$' + string(num) + '\pi$';
    else
        if num < 0
            str = '$-\frac{' + string(abs(num)) + '}{' + string(den) + '}\pi$';
        else
            str = '$\frac{' + string(num) + '}{' + string(den) + '}\pi$';
        end
    end
end

function Y = piRationM(X)
    Y = strings(size(X));
    for i = 1:numel(X)
        Y(i) = piRation(X(i));
    end
end
