% 椅子に座った状態で片脚の膝の曲げ伸ばしを行った。
clear;
close all;
clc;

%% パス追加
addpath('KIT実験');
addpath('富大山内研');

%% Excel 読み込み
T = readtable('1link_KIT_08.xlsx');

%% 脛の剛体マーカ重心位置
x = T{9:end,26};
y = T{9:end,27};
z = T{9:end,28};

p = [x y z];

%% 脛クォータニオン
xq = T{9:end,22};
yq = T{9:end,23};
zq = T{9:end,24};
wq = T{9:end,25};

q = [wq xq yq zq];

%% 時間データ
t = T{9:end,2};

%% データ数
N = length(x);

%% 回転行列へ変換
R = zeros(3,3,N);

for k = 1:N
    
    % MATLAB quaternion形式 [w x y z]
    quat = quaternion(q(k,:));
    
    % 回転行列
    R(:,:,k) = rotmat(quat,'frame');
    
end

%% 相対変換計算
dR = zeros(3,3,N-1);

for k = 1:N-1

    R1 = R(:,:,1);      % 初期姿勢を基準
    R2 = R(:,:,k+1);

    % 相対回転 R(theta)
    dR(:,:,k) = R2 * R1';

end

%% 膝関節位置推定
pwk = zeros(N-1,3);

p0 = p(1,:)';   % 初期脛位置

I = eye(3);

for k = 1:N-1

    pt = p(k+1,:)';

    Rtheta = dR(:,:,k);

    % (I-R)p = pt - R*p0
    A = I - Rtheta;
    b = pt - Rtheta*p0;

    % 特異回避
    if rank(A) < 3
        pwk(k,:) = [NaN NaN NaN];
    else
        p_est = A\b;
        pwk(k,:) = p_est';
    end

end

%% 平均値を回転中心として採用
pwk_mean = mean(pwk,'omitnan');

disp('推定膝位置')
disp(pwk_mean)

%% 膝位置推移の描画

t_plot = t(2:end);   % pwkはN-1個なので時間も合わせる

figure;

subplot(3,1,1)
plot(t_plot,pwk(:,1))
hold on
yline(pwk_mean(1),'--')
xlabel('Time [s]')
ylabel('X [mm]')
title('Estimated Knee Position X')
grid on

subplot(3,1,2)
plot(t_plot,pwk(:,2))
hold on
yline(pwk_mean(2),'--')
xlabel('Time [s]')
ylabel('Y [mm]')
title('Estimated Knee Position Y')
grid on

subplot(3,1,3)
plot(t_plot,pwk(:,3))
hold on
yline(pwk_mean(3),'--')
xlabel('Time [s]')
ylabel('Z [mm]')
title('Estimated Knee Position Z')
grid on

%% 3次元表示

figure

plot3(pwk(:,1),pwk(:,2),pwk(:,3))
hold on

plot3(pwk_mean(1),...
      pwk_mean(2),...
      pwk_mean(3),...
      'o','MarkerSize',10)

xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')

title('Estimated Knee Center Trajectory')
grid on
axis equal