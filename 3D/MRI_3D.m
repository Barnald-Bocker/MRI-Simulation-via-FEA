clearvars;

%-----------------------------
% Estableciendo geometrias
% L = 0.4;
% H = 0.4;
% S1 = [3;4;-L/2;L/2;L/2;-L/2;-H/2;-H/2;H/2;H/2];
% gd = S1;
% ns = char('S1');
% ns = ns';
% sf = 'S1';
% [dl,bt] = decsg(gd,sf,ns);

% Graficar el dominio
%-----------------------------------
% figure(1)
% pdegplot(dl, "EdgeLabels","on","FaceLabels","on")

% Crear modelo (clase) ecuaciones diferenciales de matlab
model = createpde(3); % 3 equations for Mx, My, Mz
importGeometry(model,"ModeloPavo.stl")
generateMesh(model, 'Hmax', 5);
figure(1)
pdemesh(model)
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
    f(3,:) = (1/T1) * (1+1*exp((-(location.x).^2 - (location.y).^2)/b) - Mz);             % dMz/dt 
end

% Establecer coeficientes PDE
specifyCoefficients(model, 'm',0, 'd',1, 'c',0, 'a',0, 'f',@fcoeffunction);
% Condiciones iniciales (3 componentes)

b = 1e-2;
u0 = @(location) [ ...
    ones(1, numel(location.x)); ...  % Mx initial = 1
    ones(1, numel(location.x)); ...  % My initial = 1
    1+1*exp((-(location.x).^2 - (location.y).^2)/b) ... % Mz inicial = 1 + perturbación
];

setInitialConditions(model, u0);

% Generar malla y resolver
t = 0:0.001:3;
R = solvepde(model, t);
% Extract solution data
u = R.NodalSolution;  % [numNodes x 3 x numTimeSteps]
t = R.SolutionTimes;  % Time points

%% 1. Animation of Mz spatial distribution - FINAL CORRECTED VERSION
figure(2);

% Extract surface mesh
mesh = model.Mesh;
vertices = mesh.Nodes';        % Nx3 vertex coordinates
elements = mesh.Elements';     % Tetrahedral connectivity

% Find surface faces (boundary of tetrahedral mesh)
faces = meshBoundaryFaces(elements);

% Create patch object
patchHandle = patch('Faces', faces, ...
                   'Vertices', vertices, ...
                   'FaceColor', 'flat', ...
                   'EdgeColor', 'none', ...
                   'FaceVertexCData', u(:,3,1));
axis equal;
view(3);
title('Mz Distribution at t=0');
clim([min(u(:,3,:),[],'all'), max(u(:,3,:),[],'all')]);
colormap('gray');
colorbar;

% Set up animation parameters
numFrames = min(50, length(t));  % Max 50 frames
frameIndices = round(linspace(1, length(t), numFrames));

% Capture and write first frame
frame = getframe(gcf);
im = frame2im(frame);
[A, map] = rgb2ind(im, 256);
imwrite(A, map, 'evolución_mag_3D.gif', 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);

% Update and write subsequent frames
for idx = 2:numFrames
    k = frameIndices(idx);
    % Update patch data
    set(patchHandle, 'FaceVertexCData', u(:,3,k));
    title(sprintf('Mz Distribution at t = %.3f', t(k)));
    drawnow;
    
    % Capture frame and add to GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    imwrite(A, map, 'evolución_mag_3D.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    pause(0.1)
end

% Helper function to find boundary faces
function faces = meshBoundaryFaces(elements)
    % Find all triangular faces of tetrahedrons
    allFaces = [ 
        elements(:,[1 2 3]);
        elements(:,[1 2 4]);
        elements(:,[1 3 4]);
        elements(:,[2 3 4])
    ];
    
    % Sort each face to ensure consistent ordering
    allFaces = sort(allFaces, 2);
    
    % Find unique faces and their counts
    [uniqueFaces, ~, ic] = unique(allFaces, 'rows');
    faceCounts = accumarray(ic, 1);
    
    % Boundary faces are those that appear only once
    faces = uniqueFaces(faceCounts == 1, :);
end
%% 2. Time evolution at domain center
% Find node closest to (0,0)
p = model.Mesh.Nodes; % Get node coordinates from the mesh
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
saveas(gcf, 'Evolucion_Magnetizacion_Centro_3D.png');
%% 3. Time evolution of spatial average
avgM = squeeze(mean(u, 1))'; % [time x components]

figure(4);
plot(t, avgM(:,1), 'r', t, avgM(:,2), 'g', t, avgM(:,3), 'b');
title('Promedio Espacial de la Magnetización');
xlabel('Tiempo');
ylabel('Magnetización Promedio');
legend('Mx', 'My', 'Mz');
grid on;
saveas(gcf, 'Promedio_Espacial_Magnetizacion_3D.png');