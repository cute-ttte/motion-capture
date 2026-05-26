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
    R(:,:,k) = rotmat(quat,'point');
end

%% ==========================
%% 1. 膝関節位置 pwk の推定
%% ==========================
pwf_0 = p(1,:)';
Rwf_0 = R(:,:,1);

%% 膝関節姿勢 Rwk 推定

Rwk = zeros(3,3,N);

for k = 2:N

    % 現在姿勢
    Rwf_t = R(:,:,k);

    % 相対回転
    Rf = Rwf_t*Rwf_0';

    %% 回転軸ベクトル抽出

    v = [
        Rf(3,2)-Rf(2,3)
        Rf(1,3)-Rf(3,1)
        Rf(2,1)-Rf(1,2)
        ];

    % 回転軸方向
    ky = v/norm(v);

    %% グラムシュミット

    a = [1;0;0];

    % x軸
    kx = a - (a'*ky)*ky;
    kx = kx/norm(kx);

    % z軸
    kz = cross(kx,ky);
    kz = kz/norm(kz);

    %% 回転行列

    Rwk(:,:,k) = [kx ky kz];

end

%% 膝関節角度

theta = zeros(N,1);

for k=2:N

    Rkf = Rwk(:,:,k)'*R(:,:,k);

    theta(k) = atan2(Rkf(1,3),Rkf(1,1));

end

theta = rad2deg(theta);

figure
plot(t,theta)
xlabel('Time [s]')
ylabel('Knee angle [deg]')
grid on