% 椅子に座った状態で片脚の膝の曲げ伸ばし
clear;
close all;
clc;

%% パス追加
addpath('KIT実験');
addpath('富大山内研');

%% Excel読み込み
T = readtable('1link_KIT_08.xlsx');

%% 脛部重心位置
x = T{9:end,26};
y = T{9:end,27};
z = T{9:end,28};

p = [x y z];

%% クォータニオン
xq = T{9:end,22};
yq = T{9:end,23};
zq = T{9:end,24};
wq = T{9:end,25};

% MATLAB形式 [w x y z]
q = [wq xq yq zq];

%% 時間
t = T{9:end,2};

%% データ数
N = length(t);

%% 回転行列へ変換
R = zeros(3,3,N);

for k = 1:N

    quat = quaternion(q(k,:));

    % 回転行列
    R(:,:,k) = rotmat(quat,'point');

end

%% 膝関節位置（回転中心）推定

% 初期姿勢
pwf_0 = p(1,:)';
Rwf_0 = R(:,:,1);

A = [];
B = [];

I = eye(3);

for k = 2:N

    %% 現在時刻
    pwf_t = p(k,:)';
    Rwf_t = R(:,:,k);

    %% 相対回転行列
    R_theta = Rwf_t * Rwf_0';

    %% 回転角算出
    theta = acos( ...
        max(min((trace(R_theta)-1)/2,1),-1));

    %% 小さい回転を除去
    if theta < deg2rad(5)
        continue
    end

    %% 最小二乗行列作成
    A_k = I - R_theta;

    B_k = pwf_t - R_theta*pwf_0;

    A = [A;A_k];
    B = [B;B_k];

end

%% ランク確認
fprintf('rank(A)=%d\n',rank(A));

%% 最小ノルム最小二乗解
pwk = lsqminnorm(A,B);

fprintf('\n推定膝関節位置\n')
fprintf('X = %.3f mm\n',pwk(1))
fprintf('Y = %.3f mm\n',pwk(2))
fprintf('Z = %.3f mm\n',pwk(3))

%% ==========================
%% 可視化
%% ==========================

figure
hold on
grid on
axis equal
view(3)

%% 脛部重心軌跡
plot3( ...
    p(:,1),...
    p(:,2),...
    p(:,3),...
    'LineWidth',2)

%% 推定膝位置
scatter3( ...
    pwk(1),...
    pwk(2),...
    pwk(3),...
    200,...
    'filled')

%% 初期位置
scatter3( ...
    p(1,1),...
    p(1,2),...
    p(1,3),...
    100,...
    'filled')

xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')

title('Knee joint center estimation')

legend( ...
    'Shank trajectory',...
    'Estimated knee',...
    'Initial position')

rotate3d on