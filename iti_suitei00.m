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

% 推定値保存
pwk = zeros(3,N);

for k = 2:N

    % 現在姿勢
    pwf_t = p(k,:)';
    Rwf_t = R(:,:,k);

    % 相対回転
    Rf = Rwf_t * Rwf_0';

    % 左辺
    A = eye(3) - Rf;

    % 右辺
    B = pwf_t - Rf*pwf_0;

    % 単純解
    pwk(:,k) = A\B;

end

%% x,y,zを別々表示

figure

subplot(3,1,1)
plot(t,pwk(1,:))
ylabel('x [m]')
grid on

subplot(3,1,2)
plot(t,pwk(2,:))
ylabel('y [m]')
grid on

subplot(3,1,3)
plot(t,pwk(3,:))
ylabel('z [m]')
xlabel('Time [s]')
grid on
