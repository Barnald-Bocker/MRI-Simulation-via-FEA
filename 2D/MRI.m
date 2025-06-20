clearvars;

%-----------------------------
% Estableciendo geometrias
L = 0.4;
H = 0.4;
S1 = [3;4;-L/2;L/2;L/2;-L/2;-H/2;-H/2;H/2;H/2];
S2 = [3;4;-L/2+1/2;L/2+1/2;L/2+1/2;-L/2+1/2;-H/2;-H/2;H/2;H/2];
S3 = [3;4;-L/2;L/2;L/2;-L/2;-H/2+1/2;-H/2+1/2;H/2+1/2;H/2+1/2];
S4 = [3;4;-L/2+1/2;L/2+1/2;L/2+1/2;-L/2+1/2;-H/2+1/2;-H/2+1/2;H/2+1/2;H/2+1/2];
gd = [S1,S2,S3,S4];
ns = char('S1','S2','S3','S4');
ns = ns';
sf = 'S1+S2+S3+S4';
[dl,bt] = decsg(gd,sf,ns);

% Graficar el dominio
%-----------------------------------
% figure(1)
% pdegplot(dl, "EdgeLabels","on","FaceLabels","on")

% Crear modelo (clase) ecuaciones diferenciales de matlab
model = createpde(3); % 3 equations for Mx, My, Mz

% Se ingresa la geometría creada al modelo
geometryFromEdges(model,dl);



% Definir función fcoeffunction correctamente
function f = fcoeffunction1(location, state)

    b = 2e-2;
    % Parámetros de las ecuaciones de Bloch
    omega0 = 42.58;
    T2 = 0.1;
    T1 = 0.92;
    
    N = 3; % número de ecuaciones
    nr = numel(location.x); % número de puntos
    f = zeros(N, nr);
    
    if nargin < 2 || isempty(state)
        f = NaN(N, nr); % Return NaN when state is not provided
        return;
    end
    
    % state.u contiene la solución: [Mx; My; Mz] (3 filas x nr columnas)
    Mx = state.u(1,:);
    My = state.u(2,:);
    Mz = state.u(3,:);
    
    % Ecuaciones de Bloch
    f(1,:) = omega0 * My - (1/T2) * Mx;    % dMx/dt
    f(2,:) = -omega0 * Mx - (1/T2) * My;   % dMy/dt
    f(3,:) = (1/T1) * (1*exp((-(location.x).^2 - (location.y).^2)/b) - Mz);             % dMz/dt 
end

function f = fcoeffunction2(location, state)

    b = 2e-2;
    % Parámetros de las ecuaciones de Bloch
    omega0 = 42.58;
    T2 = 3.7;
    T1 = 3.7;
    
    N = 3; % número de ecuaciones
    nr = numel(location.x); % número de puntos
    f = zeros(N, nr);
    
    if nargin < 2 || isempty(state)
        f = NaN(N, nr); % Return NaN when state is not provided
        return;
    end
    
    % state.u contiene la solución: [Mx; My; Mz] (3 filas x nr columnas)
    Mx = state.u(1,:);
    My = state.u(2,:);
    Mz = state.u(3,:);
    
    % Ecuaciones de Bloch
    f(1,:) = omega0 * My - (1/T2) * Mx;    % dMx/dt
    f(2,:) = -omega0 * Mx - (1/T2) * My;   % dMy/dt
    f(3,:) = (1/T1) * (1*exp((-(location.x-1/2).^2 - (location.y).^2)/b) - Mz);             % dMz/dt 
end

function f = fcoeffunction3(location, state)

    b = 2e-2;
    % Parámetros de las ecuaciones de Bloch
    omega0 = 42.58;
    T2 = 0.05;
    T1 = 0.9;
    
    N = 3; % número de ecuaciones
    nr = numel(location.x); % número de puntos
    f = zeros(N, nr);
    
    if nargin < 2 || isempty(state)
        f = NaN(N, nr); % Return NaN when state is not provided
        return;
    end
    
    % state.u contiene la solución: [Mx; My; Mz] (3 filas x nr columnas)
    Mx = state.u(1,:);
    My = state.u(2,:);
    Mz = state.u(3,:);
    
    % Ecuaciones de Bloch
    f(1,:) = omega0 * My - (1/T2) * Mx;    % dMx/dt
    f(2,:) = -omega0 * Mx - (1/T2) * My;   % dMy/dt
    f(3,:) = (1/T1) * (1*exp((-(location.x).^2 - (location.y-1/2).^2)/b) - Mz);             % dMz/dt 
end

function f = fcoeffunction4(location, state)

    b = 2e-2;
    % Parámetros de las ecuaciones de Bloch
    omega0 = 42.58;
    T2 = 0.08;
    T1 = 0.25;
    
    N = 3; % número de ecuaciones
    nr = numel(location.x); % número de puntos
    f = zeros(N, nr);
    
    if nargin < 2 || isempty(state)
        f = NaN(N, nr); % Return NaN when state is not provided
        return;
    end
    
    % state.u contiene la solución: [Mx; My; Mz] (3 filas x nr columnas)
    Mx = state.u(1,:);
    My = state.u(2,:);
    Mz = state.u(3,:);
    
    % Ecuaciones de Bloch
    f(1,:) = omega0 * My - (1/T2) * Mx;    % dMx/dt
    f(2,:) = -omega0 * Mx - (1/T2) * My;   % dMy/dt
    f(3,:) = (1/T1) * (1*exp((-(location.x-1/2).^2 - (location.y-1/2).^2)/b) - Mz);             % dMz/dt 
end
% Establecer coeficientes PDE
specifyCoefficients(model, 'm',0, 'd',1, 'c',0, 'a',0, 'f',@fcoeffunction1,'Face',1);
specifyCoefficients(model, 'm',0, 'd',1, 'c',0, 'a',0, 'f',@fcoeffunction2,'Face',2);
specifyCoefficients(model, 'm',0, 'd',1, 'c',0, 'a',0, 'f',@fcoeffunction3,'Face',3);
specifyCoefficients(model, 'm',0, 'd',1, 'c',0, 'a',0, 'f',@fcoeffunction4,'Face',4);
% Condiciones iniciales (3 componentes)
a = -0.05;
b = 3e-2;
u01 = @(location) [ ...
    ones(1, numel(location.x)); ...  % Mx initial = 0
    ones(1, numel(location.x)); ...  % My initial = 0
    1+1*exp((-(location.x).^2 - (location.y).^2)/b) ... % Mz inicial = 1 + perturbación
];

u02 = @(location) [ ...
    ones(1, numel(location.x)); ...  % Mx initial = 0
    ones(1, numel(location.x)); ...  % My initial = 0
    1+1*exp((-(location.x-1/2).^2 - (location.y).^2)/b) ... % Mz inicial = 1 + perturbación
];

u03 = @(location) [ ...
    ones(1, numel(location.x)); ...  % Mx initial = 0
    ones(1, numel(location.x)); ...  % My initial = 0
    1+1*exp((-(location.x).^2 - (location.y-1/2).^2)/b) ... % Mz inicial = 1 + perturbación
];

u04 = @(location) [ ...
    ones(1, numel(location.x)); ...  % Mx initial = 0
    ones(1, numel(location.x)); ...  % My initial = 0
    1+1*exp((-(location.x-1/2).^2 - (location.y-1/2).^2)/b) ... % Mz inicial = 1 + perturbación
];

setInitialConditions(model, u01,Face=1);
setInitialConditions(model, u02,Face=2);
setInitialConditions(model, u03,Face=3);
setInitialConditions(model, u04,Face=4);

% Generar malla y resolver
generateMesh(model, 'Hmax', 0.05);
t = 0:0.001:5;
R = solvepde(model, t);
% Extract solution data
u = R.NodalSolution;  % [numNodes x 3 x numTimeSteps]
t = R.SolutionTimes;  % Time points

%% 1. Animation of Mz (sin cambios - muestra todo el dominio)
titles = {'Cara 1 (Materia Gris)', 'Cara 2 (Agua)', 'Cara 3 (Músculo)', 'Cara 4 (Tejido Adiposo)'};
centers = [0, 0;   % Centro de S1 (Face 1)
           1/2, 0;   % Centro de S2 (Face 2)
           0, 1/2;   % Centro de S3 (Face 3)
           1/2, 1/2];  % Centro de S4 (Face 4)

positions = [0, +0.18;   % Centro de S1 (Face 1)
           1/2, +0.18;   % Centro de S2 (Face 2)
           0, 1/2+0.18;   % Centro de S3 (Face 3)
           1/2, 1/2+0.18];  % Centro de S4 (Face 4)
tissueNames = {'Materia Gris', 'Agua', 'Músculo', 'Tejido Adiposo'};
figure(2);
p = model.Mesh.Nodes; 
pdeplot(model, 'XYData', u(:,3,1), 'ColorMap', 'gray', 'Contour', 'off');
title(sprintf('Mz Distribution at t = %.3f', t(1)));
clim([min(u(:,3,:), [], 'all'), max(u(:,3,:), [], 'all')]);
colorbar;
ax = gca;
patchHandle = findobj(ax, 'Type', 'patch');

textHandles = gobjects(4,1);
for i = 1:4
    textHandles(i) = text(positions(i,1), positions(i,2), tissueNames{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'BackgroundColor','white',...
        'EdgeColor', 'black', ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'Margin', 1);
end

for k = 1:10:length(t)
    if ~isempty(patchHandle)
        set(patchHandle, 'CData', u(:,3,k));
    end
    title(sprintf('Distribución de Mz en t = %.3f', t(k)));
    pause(0.01)
    exportgraphics(gcf,'evolución_mag.gif','Append',true)
    drawnow;
end

%% 2. Time evolution at the center of each square

centerIdx = zeros(4,1);
p = model.Mesh.Nodes';
for i = 1:4
    dist = sum((p - centers(i,:)).^2, 2);
    [~, idx] = min(dist);
    centerIdx(i) = idx;
end

% Extract solution at centers
centerSol = zeros(length(t), 3, 4); % [time, component, square]
for i = 1:4
    centerSol(:,:,i) = squeeze(u(centerIdx(i), :, :))';
end

figure(3);

for i = 1:4
    subplot(2,2,i);
    plot(t, centerSol(:,1,i), 'r', t, centerSol(:,2,i), 'g', t, centerSol(:,3,i), 'b');
    title(sprintf('Evolución de %s', titles{i}));
    xlabel('Tiempo');
    ylabel('Magnetización');
    legend('Mx', 'My', 'Mz');
    grid on;
end
saveas(gcf, 'Evolucion_Centros.png');