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
% ΔT_k = T_(k+1) * inv(T_k)

dR = zeros(3,3,N-1);
dt = zeros(N-1,3);

for k = 1:N-1
    
    R1 = R(:,:,k);
    R2 = R(:,:,k+1);
    
    p1 = p(k,:)';
    p2 = p(k+1,:)';
    
    % 相対回転
    dR(:,:,k) = R2 * R1';
    
    % 相対並進
    dt(k,:) = (p2 - dR(:,:,k)*p1)';
    
end

%% 各フレームの瞬間回転軸方向
u_all = zeros(N-1,3);

for k = 1:N-1
    
    A = dR(:,:,k);
    
    % 固有値・固有ベクトル
    [V,D] = eig(A);
    
    eigvals = diag(D);
    
    % 固有値1に最も近いもの
    [~,idx] = min(abs(eigvals - 1));
    
    u = real(V(:,idx));
    
    % 正規化
    u = u / norm(u);
    
    u_all(k,:) = u';
    
end

%% 軸方向平均
u_mean = mean(u_all,1)';
u_mean = u_mean / norm(u_mean);

disp('推定回転軸方向');
disp(u_mean);

%% 軸上点 c を最小二乗推定
%
% (I - dR_k)c = dt_k
%

Aall = [];
ball = [];

for k = 1:N-1
    
    A = eye(3) - dR(:,:,k);
    b = dt(k,:)';
    
    Aall = [Aall; A];
    ball = [ball; b];
    
end

% 最小二乗解
c = pinv(Aall) * ball;

disp('推定膝関節中心');
disp(c);

%% 各フレームの膝角度推定
theta = zeros(N,1);

R0 = R(:,:,1);

for k = 1:N
    
    Rrel = R(:,:,k) * R0';
    
    tr = trace(Rrel);
    
    val = (tr - 1)/2;
    
    % 数値誤差対策
    val = max(min(val,1),-1);
    
    theta(k) = acos(val);
    
end

% degree変換
theta_deg = rad2deg(theta);

%% 可視化
figure;

plot(t, theta_deg,'LineWidth',2);

xlabel('Time [s]');
ylabel('Knee Flexion Angle [deg]');
title('Estimated Knee Flexion Angle');

grid on;

%% 3D表示
figure;
hold on;
grid on;
axis equal;

% 脛重心軌跡
plot3(p(:,1),p(:,2),p(:,3),'b');

% 推定膝関節中心
scatter3(c(1),c(2),c(3),100,'r','filled');

% 回転軸描画
L = 300;

p1 = c - L*u_mean;
p2 = c + L*u_mean;

plot3([p1(1) p2(1)], ...
      [p1(2) p2(2)], ...
      [p1(3) p2(3)], ...
      'k','LineWidth',3);

xlabel('X');
ylabel('Y');
zlabel('Z');

title('Estimated Knee Joint Axis');

legend('Tibia trajectory','Knee center','Joint axis');
view(3);