clearvars;

%-----------------------------
% Estableciendo geometrias
L = 0.4;
H = 0.4;
S1 = [3;4;-L/2;L/2;L/2;-L/2;-H/2;-H/2;H/2;H/2];
gd = S1;
ns = char('S1');
ns = ns';
sf = 'S1';
[dl,bt] = decsg(gd,sf,ns);

% Graficar el dominio
%-----------------------------------
figure(1)
pdegplot(dl, "EdgeLabels","on","FaceLabels","on")

% Crear modelo (clase) ecuaciones diferenciales de matlab
model = createpde(3); % 3 equations for Mx, My, Mz

% Se ingresa la geometría creada al modelo
geometryFromEdges(model,dl);

% Parámetros de las ecuaciones de Bloch
omega0 = 0.1;
T2 = 0.1;
T1 = 0.4;

% Definir función fcoeffunction correctamente
function f = fcoeffunction(location, state)

    b = 2e-2;
    omega0 = 42.58;
    T2 = 0.1;
    T1 = 0.4;
    
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

% Establecer coeficientes PDE
specifyCoefficients(model, 'm',0, 'd',1, 'c',0, 'a',0, 'f',@fcoeffunction);
% Condiciones iniciales (3 componentes)
a = -0.05;
b = 3e-2;
u0 = @(location) [ ...
    ones(1, numel(location.x)); ...  % Mx initial = 0
    ones(1, numel(location.x)); ...  % My initial = 0
    1+1*exp((-(location.x).^2 - (location.y).^2)/b) ... % Mz inicial = 1 + perturbación
];

setInitialConditions(model, u0);

% Generar malla y resolver
generateMesh(model, 'Hmax', 0.05);
t = 0:0.001:3;
R = solvepde(model, t);
% Extract solution data
u = R.NodalSolution;  % [numNodes x 3 x numTimeSteps]
t = R.SolutionTimes;  % Time points

%% 1. Animation of Mz spatial distribution (FIXED)
figure(2);
p = model.Mesh.Nodes; % Get node coordinates from the mesh
% First plot - get all graphics objects
pdeplot(model, 'XYData', u(:,3,1), 'ColorMap', 'gray', 'Contour', 'off');
title(sprintf('Mz Distribution at t = %.3f', t(1)));
clim([min(u(:,3,:), [], 'all'), max(u(:,3,:), [], 'all')]);
colorbar;

% Find the patch object that contains the face data
ax = gca;
patchHandle = findobj(ax, 'Type', 'patch');


for k = 1:10:length(t)/2
    % Update the face color data.
    if ~isempty(patchHandle)
        set(patchHandle, 'CData', u(:,3,k));
    end
    title(sprintf('Distribución de Mz en t = %.3f', t(k)));
    pause(0.01)
    exportgraphics(gcf,'evolución_mag.gif','Append',true)
    drawnow;
end

% close(vidObj); % Uncomment if recording video

%% 2. Time evolution at domain center
% Find node closest to (0,0)
[~, centerIdx] = min(sum(p.^2, 1));
centerSol = squeeze(u(centerIdx, :, :))'; % [time x components]

figure(3);
subplot(2,1,1);
plot(t, centerSol(:,1), 'r', t, centerSol(:,2), 'g', t, centerSol(:,3), 'b');
title('Componentes de la magnetización en el centro');
xlabel('Tiempo');
ylabel('Magnetización');
legend('Mx', 'My', 'Mz');
grid on;
savefig()
% Calculate magnitude and transverse magnetization
transverse = sqrt(centerSol(:,1).^2 + centerSol(:,2).^2);
magnitude = sqrt(transverse.^2 + centerSol(:,3).^2);

subplot(2,1,2);
plot(t, transverse, 'm', t, centerSol(:,3), 'b', t, magnitude, 'k--');
title('Magnitud de Magnetización en el centro');
xlabel('Tiempo');
ylabel('Magnitud');
legend('Transversal (M_{xy})', 'Longitudinal (M_z)', 'Magnitud Total');
grid on;
saveas(gcf, 'Evolucion_Magnetizacion_Centro.png');
%% 3. Time evolution of spatial average
avgM = squeeze(mean(u, 1))'; % [time x components]

figure(4);
plot(t, avgM(:,1), 'r', t, avgM(:,2), 'g', t, avgM(:,3), 'b');
title('Promedio Espacial de la Magnetización');
xlabel('Tiempo');
ylabel('Magnetización Promedio');
legend('Mx', 'My', 'Mz');
grid on;
saveas(gcf, 'Promedio_Espacial_Magnetizacion.png');