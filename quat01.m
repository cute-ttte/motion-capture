% 椅子に座った状態で片脚の膝の曲げ伸ばしを行った。
clear;
close all;
clc;

% Excel 読み込み
T = readtable('1link_KIT_08.xlsx');


% 膝クォータニオン
x1 = T{9:end, 3};
y1 = T{9:end, 4};
z1 = T{9:end, 5};
w1 = T{9:end, 6};

% 脛クォータニオン
x2 = T{9:end, 22};
y2 = T{9:end, 23};
z2 = T{9:end, 24};
w2 = T{9:end, 25};

% データ数
N = length(x1);           

% 配列確保
yaw   = zeros(N,1);
pitch = zeros(N,1);
roll  = zeros(N,1);


for i = 1:N
    qx1 = x1(i);
    qy1 = y1(i);
    qz1 = z1(i);
    qw1 = w1(i);

    qx2 = x2(i);
    qy2 = y2(i);
    qz2 = z2(i);
    qw2 = w2(i);

    % 正規化
    norm1 = sqrt(qx1^2 + qy1^2 + qz1^2 + qw1^2);
    qx1 = qx1 / norm1;
    qy1 = qy1 / norm1;
    qz1 = qz1 / norm1;
    qw1 = qw1 / norm1;
    
    norm2 = sqrt(qx2^2 + qy2^2 + qz2^2 + qw2^2);
    qx2 = qx2 / norm2;
    qy2 = qy2 / norm2;
    qz2 = qz2 / norm2;
    qw2 = qw2 / norm2;

    % 回転行列 膝
    R1 = [ 1 - 2*(qy1^2 + qz1^2),      2*(qx1*qy1 - qz1*qw1),        2*(qx1*qz1 + qy1*qw1);
          2*(qx1*qy1 + qz1*qw1),          1 - 2*(qx1^2 + qz1^2),    2*(qy1*qz1 - qx1*qw1);
          2*(qx1*qz1 - qy1*qw1),          2*(qy1*qz1 + qx1*qw1),        1 - 2*(qx1^2 + qy1^2)];

    % 回転行列 脛
    R2 = [ 1 - 2*(qy2^2 + qz2^2),      2*(qx2*qy2 - qz2*qw2),        2*(qx2*qz2 + qy2*qw2);
          2*(qx2*qy2 + qz2*qw2),          1 - 2*(qx2^2 + qz2^2),    2*(qy2*qz2 - qx2*qw2);
          2*(qx2*qz2 - qy2*qw2),          2*(qy2*qz2 + qx2*qw2),        1 - 2*(qx2^2 + qy2^2)];

    % 大腿から見た脛回転行列
    R0 = R1.' * R2;

    %  オイラー角変換
    yaw(i)   = atan2(R0(2,1), R0(1,1));
    pitch(i) = asin(-R0(3,1));
    roll(i)  = atan2(R0(3,2), R0(3,3));
    


end


% degに変換
roll_deg  = rad2deg(roll);
pitch_deg = rad2deg(pitch) - 90;
yaw_deg   = rad2deg(yaw);

% 時間データ
t = T{9:end, 2};

% グラフ出力
figure % yaw
    plot(t,yaw_deg , 'linestyle', '-');
    yline(-90, '--');
    hold on

    xlabel('t[s]', 'Interpreter', 'latex');
    ylabel('yaw[deg]', 'Interpreter', 'latex');

figure % pitch
    plot(t,pitch_deg ,'g', 'linestyle', '-');
    yline(-90, '--');
    hold on

    xlabel('t[s]', 'Interpreter', 'latex');
    ylabel('q[deg]', 'Interpreter', 'latex');


figure % roll
    plot(t,roll_deg ,'m', 'linestyle', '-');
    yline(-90, '--');
    hold on

    xlabel('t[s]', 'Interpreter', 'latex');
    ylabel('roll[deg]', 'Interpreter', 'latex');

% 初期角度
pitch_deg(1)