% 椅子に座った状態で片脚の膝の曲げ伸ばし
clear;
close all;
clc;

%% パス追加
addpath('KIT実験');
addpath('富大山内研');

%% Excel読み込み
T = readtable('1link_KIT_08.xlsx');

%% 脛剛体マーカ重心位置
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
    R(:,:,k) = rotmat(quat,'frame');

end

%% 膝関節位置（回転軸）の推定
pwf_0 = p(1,:)';
Rwf_0 = R(:,:,1);

A = [];
B = [];

I = eye(3);

for k = 1:N

    % 現在時刻
    pwf_t = p(k,:)';
    Rwf_t = R(:,:,k);

    %% 相対回転行列
    R_theta = Rwf_t * Rwf_0';

    %% 回転角計算
    val = (trace(R_theta)-1)/2;

    % 数値誤差対策
    val = max(min(val,1),-1);

    theta = acos(val);

    %% 小さい回転を除外
    if theta < deg2rad(5)
        continue;
    end

    %% 最小二乗行列作成
    A = [A;
         R_theta-I];

    B = [B;
         R_theta*pwf_0-pwf_t];

end

%% 最小二乗解
pwk = A\B;

%% 結果表示
disp('====================================');

fprintf('推定された膝関節(回転軸)の位置 [m]\n');

fprintf('x = %.5f\n',pwk(1));
fprintf('y = %.5f\n',pwk(2));
fprintf('z = %.5f\n',pwk(3));

disp('====================================');

%% 時系列表示
figure('Name','Position Time Series','Color','w');

axis_labels = {'X Position [m]','Y Position [m]','Z Position [m]'};

for i=1:3

    subplot(3,1,i)

    plot(t,p(:,i),'b','LineWidth',1.5)
    hold on

    plot([t(1),t(end)],...
         [pwk(i),pwk(i)],...
         'r--','LineWidth',2)

    grid on

    xlabel('Time [s]')
    ylabel(axis_labels{i})

    legend('Foot Marker',...
           'Estimated Knee',...
           'Location','best')

end

%% ======3次元表示======

figure('Name',...
       'Estimated Knee Position',...
       'Color','w')

% 脛部マーカ点群
plot3(p(:,1),...
      p(:,2),...
      p(:,3),...
      'b.',...
      'MarkerSize',10)

hold on

% 推定膝位置
plot3(pwk(1),...
      pwk(2),...
      pwk(3),...
      'r.',...
      'MarkerSize',40)

grid on
axis equal

xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

title('Estimated Knee Position and Foot Marker Trajectory')

legend('Foot Marker',...
       'Estimated Knee')

view(3)

rotate3d on

%% ======膝からの距離確認======

r = vecnorm((p-pwk'),2,2);

figure('Name',...
       'Distance from Knee',...
       'Color','w')

plot(t,r,'LineWidth',1.5)

grid on

xlabel('Time [s]')
ylabel('Distance [m]')

title('Distance from Estimated Knee')