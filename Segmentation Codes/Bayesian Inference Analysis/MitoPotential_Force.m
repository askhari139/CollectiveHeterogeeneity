% ----- INPUT: 
%
% filename = filename of Matlab workspace generated from Tissue Analyser
% segmentation (imported bonds.csv and cells.csv files). 
% Tissue Analyser is a segmentation plugin for FIJI (B. Aigouy)
% https://grr.gred-clermont.fr/labmirouse/software/WebPA/
%
load('cellbonddata'); load('muFile');
evolver=0;
SquareSize = 128; nRadius = 1.5;
pixScale = 1;
isDelFreeEdge = 0; showLLL= 1; showDist = 1;
isLeaderDist = 0; leaderPos = [1014 548]; % if isLeaderDist = 1, set the leader cell position
showCellForce = 1;
%% INPUT POTENTIAL IMAGES

[nFlname, nPathname] = uigetfile('*.tif','Select the file with TMRM staining');
nFlpathname = strcat(nPathname,nFlname);
nI = imread(nFlpathname);
if size(nI,3)==3
    nI = rgb2gray(nI);
    i = 1
end
nI = double(nI); 
[n_xLen n_yLen] = size(nI);
min_I = min(min(nI)); max_I = max(max(nI));
nIs = (nI-min_I)*255/(max_I-min_I);

%% TRANSLATE TISSUE ANALYZER GEOMETRY TO MATLAB

tic;
disp('Generating geometry from Tisse Analyzer data...')
TissueAnalyzerToMatlab;

%% REMOVE FREE EDGES FROM T-P CALCULATION
% Find edges of each cell (CellEdges)
Ne=length(E); Nc=length(C);
for i=1:Nc % for each cell
    clear e
    nv=length(C{i}); % how many vertices (and edges) in the cell
    for j=1:nv-1
        e(j,:)=[C{i}(j) C{i}(j+1)]; % each couple of adjacent vertices is an edge
        for k=1:Ne % find which edge
            if (E(k,1)== e(j,1) && E(k,2)== e(j,2)) || (E(k,1)== e(j,2) && E(k,2)== e(j,1))
                CellEdges{i}(j)=k;
            end
        end
    end
    e(nv,:)=[C{i}(nv) C{i}(1)];
    for k=1:Ne
        if (E(k,1)== e(nv,1) && E(k,2)== e(nv,2)) || (E(k,1)== e(nv,2) && E(k,2)== e(nv,1))
            CellEdges{i}(nv)=k; % List of edges (index in E) in Cell i
        end
    end
    CellEdges{i}(CellEdges{i}==0)=[]; 
    % some zeros can appear because of image boundaries, remove them !
    % (edges connecting two vertices on the image boundary 
    % are not in E, because tension is not calculated for these edges)
end

% Remove the edges associated with only one cell
cellpEdge = zeros(Ne,1);
for c=1:C_tot
    numE = length(CellEdges{c});
    if numE<10                   % cell free contributes to cell initially
        for k=1:numE
            cE = CellEdges{c}(k);
            cellpEdge(cE) = cellpEdge(cE) + 1;
        end
    end
end

edgeToDel = find(cellpEdge == 1);
verEdge = EDGES(edgeToDel,:); verEdgeN = length(verEdge)*2;
verEdgeList = reshape(verEdge,[verEdgeN,1]);
verEdge = unique(verEdgeList); verEdgeN = length(verEdge);
vEx = real(VERTICES(verEdge)); vEy = imag(VERTICES(verEdge));
vE = [vEx vEy];
% 
if isDelFreeEdge == 1
    E(edgeToDel,:) = []; E_tot = length(E);
end

%% COMPUTATION OF THE COEFFICIENTS OF MATRIX A 

[A,B,g]=subprog_Calcul_Matrix_A(V, E, C, V_ent, E_tot, C_tot,evolver);
disp(['Done (' num2str(round(toc,2)) ' seconds).'])

%% INVERSE PROBLEM

tic;
disp('Solving the inverse problem...')
[INFERENCE] = subprog_inverse_matrix_A(mu_bayes, A, B, g, E_tot);
disp(['Done (' num2str(round(toc,2)) ' seconds).'])

%% FIGURE / TENSIONS
% this should be customized depending on user's needs

figure; set(gcf,'color','w'); couleur = jet(101); hold on

tmin=1-3.5*std(INFERENCE.TENSIONS); % 3*std : good compromise for color limits
tmax=1+2.5*std(INFERENCE.TENSIONS);
Ti = (INFERENCE.TENSIONS-tmin)./(tmax-tmin);
Ti(Ti<0)=0;
Ti(Ti>1)=1;
Ti = round(Ti*100)+1; %scale to jet101
%color scaled from tmin to tmax
colormap jet(101)
cbr = colorbar;
cbr.Ticks=[0 1];
cbr.TickLabels={num2str(round(tmin,1)) num2str(round(tmax,1))};

for cpt = 1:E_tot,
    plot(V(E(cpt,:),1),V(E(cpt,:),2),...
        'Color',couleur(Ti(cpt),:),'LineWidth',2)
end
axis equal; axis tight; title('Bayesian Inference Tension Map')
axis off
set(gca,'Ydir','reverse');
%saveas(gcf,'tensions.fig');

%% FIGURE / PRESSURES
% this should be customized depending on user's needs

figure; set(gcf,'color','w'); couleur = cool(101);  hold on; 
% Change colormap for representation (jet) and computing (gray)

c_len = min(length(INFERENCE.PRESSURES),C_ent);
P_plot=INFERENCE.PRESSURES(1:c_len);
Pmax=3*std(P_plot); % 3*std : good compromise for color limits
Pmin=-Pmax;
Pi = (P_plot-Pmin)./(Pmax-Pmin);
Pi(Pi<0)=0;
Pi(Pi>1)=1;
Pi = round(Pi*100)+1; %scale to jet101
colormap(couleur)
cbr = colorbar;
cbr.Ticks=[0 1 2];
cbr.TickLabels={num2str(round(Pmin,2)) '0' num2str(round(Pmax,2))};

for e = 1:E_tot
    plot(V(E(e,:),1),V(E(e,:),2),'k','LineWidth',1)
end
for c = 1:c_len
    v_in_f = C{c};
    patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',couleur(Pi(c),:))
end

axis equal; axis tight; title('Bayesian Inference Pressure Map')
axis off
set(gca,'Ydir','reverse');
%saveas(gcf,'pressures.fig');

%% COMPUTE MEAN CELL DIAMETERES
cellArea = zeros(c_len,1);
for i=1:c_len
    CellArea(i)=polyarea(V(C{i},1),V(C{i},2));
end
meanDia = 2*((mean(CellArea)/pi)^0.5);

%% OPTION: COMPUTE STRESS TENSOR (BATCHELOR)

if (isDelFreeEdge == 0)&&(showCellForce == 0)
    tic;
    disp('Computing the stress tensor...')
    [N,S,P,BoxCenter]=batchelor(INFERENCE.TENSIONS,INFERENCE.PRESSURES,C,E,V, SquareSize);
    %[N,S,P,BoxCenter]=batchelor(ones(1,length(INFERENCE.TENSIONS)),zeros(1,length(INFERENCE.PRESSURES)),C,E,V, SquareSize);
    disp(['Done (' num2str(round(toc,2)) ' seconds).'])
%     saveas(gcf,'stress.fig');
    Bayesian_Inference.Batchelor.N=N;
    Bayesian_Inference.Batchelor.S=S;
    Bayesian_Inference.Batchelor.P=P;
    Bayesian_Inference.Batchelor.BoxCenter=BoxCenter;
    set(gca,'Ydir','reverse');
    
    lenN = length(N);
    sAvg = zeros(1,lenN); Xgrid = zeros(1,lenN); Ygrid = zeros(1,lenN);
    for i=1:lenN
        sAvg(i)=(N{i}(1,1) + N{i}(2,2))/2;
        Xgrid(i) = BoxCenter{i}(1); Ygrid(i) = BoxCenter{i}(2);
    end
    Xg = linspace(min(Xgrid), max(Xgrid), 100);
    Yg = linspace(min(Ygrid), max(Ygrid), 100);
    [xi, yi] = meshgrid(Xg, Yg);
    F = scatteredInterpolant(Xgrid',Ygrid',sAvg');
    zi = F(xi,yi);
    
    figure;
    surf(xi,yi,zi, 'EdgeAlpha', 0)
    colormap jet
    
elseif (isDelFreeEdge == 0)&&(showCellForce == 1)
    tic;
    disp('Computing the stress tensor at cell centroids...')
    [N,S,PI,theta1N,theta2N,lambdasN,CellCenter]=...
        cellbatchelor(INFERENCE.TENSIONS,INFERENCE.PRESSURES,C,E,V,C_ent,CellEdges,nRadius,n_xLen,n_yLen);
    disp(['Done (' num2str(round(toc,2)) ' seconds).'])
else
    disp('Not computing Batchelor stresses')
end

%% OUTPUT of Force Inference Calculations

Bayesian_Inference.Infered_Tensions=INFERENCE.TENSIONS;
Bayesian_Inference.Infered_Pressures=INFERENCE.PRESSURES;
Bayesian_Inference.Geometry.Cells=C;
Bayesian_Inference.Geometry.Vertices=V;
Bayesian_Inference.Geometry.Edges=E;
Bayesian_Inference.mu_bayes=mu_bayes;
Cells = C;

%% Mitochondrial Potential Landscape
figure; 
tic
disp('Calculating the Mitochondrial Potential Landscape ...')
hold on
pcolor = cool(255); colormap(pcolor);  
cbr = colorbar; cbr.Ticks=[0 1 2]; cbr.TickLabels={'0' '128' '255'};
amp_fac = 3;

for e = 1:E_tot,
    plot(V(E(e,:),1),V(E(e,:),2),'k','LineWidth',1)
end

mitoPot = zeros(c_len,1); 
for c = 1:c_len
    v_in_f = C{c};
    xV = V(v_in_f,1)'; yV = V(v_in_f,2)';
    BW = poly2mask(xV,yV,n_xLen,n_yLen);
    xg = 1:n_xLen; yg = 1:n_yLen; [X,Y] = meshgrid(xg,yg);
    
    currMitoPot = mean(mean(nIs.*BW))/mean(mean(BW));
    mitoPot(c) = currMitoPot;
end
maxPot = max(mitoPot); minPot = min(mitoPot);

for c=1:c_len
    v_in_f = C{c};
    cell_color = round(255*(mitoPot(c)-minPot)/(maxPot-minPot));
    cell_color = max(cell_color,1);
    patch(V(v_in_f, 1), V(v_in_f, 2), 1,...
        'FaceColor',pcolor(cell_color,:))
    hold on
end
hold off
disp(['Done (' num2str(round(toc,2)) ' seconds).'])
axis equal; axis tight; title('Mitochondrial Potential Landscape')
axis off
set(gca,'Ydir','reverse');

%% Plot Average Normal Stress
maxPI = max(PI(1:c_len)); minPI = min(PI(1:c_len));
figure; 
tic
disp('Plotting the Stress Landscape ...')
hold on
pcolor = summer(255); colormap(pcolor);  
cbr = colorbar; cbr.Ticks=[0 1 2]; cbr.TickLabels={'0' '128' '255'};
for e = 1:E_tot,
    plot(V(E(e,:),1),V(E(e,:),2),'k','LineWidth',1)
end

for c=1:c_len
    v_in_f = C{c};
    cell_color = round(255*(PI(c)-minPI)/(maxPI-minPI));
    cell_color = max(cell_color,1);
    patch(V(v_in_f, 1), V(v_in_f, 2), 1,...
        'FaceColor',pcolor(cell_color,:))
    hold on
end
hold off
disp(['Done (' num2str(round(toc,2)) ' seconds).'])
axis equal; axis tight; title('Normal Stress Landscape')
axis off
set(gca,'Ydir','reverse');


%% Plot stress tensor
% principal directions and amplitudes of the stress tensor N
% 
figure; set(gcf,'color','w'); hold on
for e = 1:Ne
    plot([V(E(e,1),1) V(E(e,2),1)],[V(E(e,1),2) V(E(e,2),2)],'color',[0.7 0.7 0.7]);
end
forceAniso = zeros(c_len,1); l1 = zeros(c_len,1); l2 = zeros(c_len,1);
for b=1:c_len
    lambdasN{b}=min(1000*lambdasN{b},100); % for plot
    xC = CellCenter{b}(1); yC = CellCenter{b}(2);
    l1(b) = lambdasN{b}(1); l2(b) = lambdasN{b}(2);
    % plot( [xC-0.5*lambdasN{b}(1)*cos(theta1N{b}) xC+0.5*lambdasN{b}(1)*cos(theta1N{b})] , [yC-0.5*lambdasN{b}(1)*sin(theta1N{b}) yC+0.5*lambdasN{b}(1)*sin(theta1N{b})],'linewidth',2,'color',[0.12 0.56 0.24])
    plot( [xC-0.5*lambdasN{b}(2)*cos(theta2N{b}) xC+0.5*lambdasN{b}(2)*cos(theta2N{b})] , [yC-0.5*lambdasN{b}(2)*sin(theta2N{b}) yC+0.5*lambdasN{b}(2)*sin(theta2N{b})],'linewidth',2,'color',[0.86 0.08 0.24])
    forceAniso(b) = (lambdasN{b}(2)-lambdasN{b}(1))/(lambdasN{b}(2)+lambdasN{b}(1));
end
axis equal; axis tight; title('Stress Tensor Map'); axis off
set(gca,'Ydir','reverse');
hold off


