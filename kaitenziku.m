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

% 各時刻の回転角を保存する配列（後で最大屈曲時を探すため）
theta_array = zeros(N, 1);

for k = 2:N
    %% 現在時刻
    pwf_t = p(k,:)';
    Rwf_t = R(:,:,k);

    %% 相対回転行列 (世界座標系基準)
    R_theta = Rwf_t * Rwf_0';

    %% 回転角算出
    theta = acos( ...
        max(min((trace(R_theta)-1)/2,1),-1));
    
    theta_array(k) = theta;

    %% 小さい回転を除去 (ノイズ対策)
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
% 1自由度の純粋な回転の場合、回転軸(線)が一意に決まるためランクは理想的には2になります
fprintf('rank(A) = %d (※理論上、純粋な1軸回転なら2になります)\n', rank(A));

%% 最小ノルム最小二乗解
% ランク落ち行列に対して、原点に最も近い回転軸上の点を計算します
pwk = lsqminnorm(A,B);

fprintf('\n推定膝関節位置\n')
fprintf('X = %.3f mm\n',pwk(1))
fprintf('Y = %.3f mm\n',pwk(2))
fprintf('Z = %.3f mm\n',pwk(3))

%% ==========================
%% 回転行列 Rwk の推定
%% ==========================

%% 最も膝が曲がった状態（回転角が最大）のフレームを抽出
[~, max_idx] = max(theta_array);

%% 相対回転行列（初期姿勢から最大屈曲時まで）
R_f = R(:,:,max_idx) * Rwf_0';

%% 回転軸ベクトル (固有値分解)
[V,D] = eig(R_f);
eigval = diag(D);
[~,idx] = min(abs(eigval-1));
nw = real(V(:,idx));

%% 固有ベクトルの向き(符号)の修正
% eig()で得られるベクトルは正負が不定なため、反対称行列成分を用いて回転方向と一致させる
v_skew = [R_f(3,2)-R_f(2,3); R_f(1,3)-R_f(3,1); R_f(2,1)-R_f(1,2)];
if dot(nw, v_skew) < 0
    nw = -nw;
end

%% 正規化
ey = nw/norm(nw);

%% 仮の基準軸
a = [0;0;1];

%% 平行回避
if abs(dot(a,ey)) > 0.9
    a = [1;0;0];
end

%% x軸 (グラム・シュミットの直交化)
ex = a - (a'*ey)*ey;
ex = ex/norm(ex);

%% z軸
ez = cross(ex,ey);
ez = ez/norm(ez);

% ※ex = cross(ey,ez) は数学的に元のexに戻るだけなので削除しました

%% 膝座標系回転行列
Rwk = [ex ey ez];

disp('推定回転行列 Rwk')
disp(Rwk)

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

figure
hold on
grid on
axis equal
view(3)

plot3(p(:,1),p(:,2),p(:,3),...
    'LineWidth',2)

scatter3(pwk(1),pwk(2),pwk(3),...
    150,...
    'filled')

L = 100;

quiver3( ...
    pwk(1),pwk(2),pwk(3),...
    L*ex(1),L*ex(2),L*ex(3), 'r', 'LineWidth', 2)

quiver3( ...
    pwk(1),pwk(2),pwk(3),...
    L*ey(1),L*ey(2),L*ey(3), 'g', 'LineWidth', 2)

quiver3( ...
    pwk(1),pwk(2),pwk(3),...
    L*ez(1),L*ez(2),L*ez(3), 'b', 'LineWidth', 2)

xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')

legend( ...
    'Shank trajectory',...
    'Knee center',...
    'x-axis',...
    'y-axis',...
    'z-axis')

title('Estimated knee coordinate')
rotate3d on
%% ==========================
%% 膝関節角度の時間変化グラフ（絶対変位）
%% ==========================
figure
plot(t, rad2deg(theta_array), 'LineWidth', 2, 'Color', '#0072BD')
grid on
xlabel('Time [s]')
ylabel('Angle displacement [deg]')
title('Knee Joint Angle over Time (Magnitude)')