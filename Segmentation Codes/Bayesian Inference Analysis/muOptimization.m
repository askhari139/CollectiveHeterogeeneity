% ----- INPUT: 
%
% filename = filename of Matlab workspace generated from Tissue Analyser
% segmentation (imported bonds.csv and cells.csv files). 
% Tissue Analyser is a segmentation plugin for FIJI (B. Aigouy)
% https://grr.gred-clermont.fr/labmirouse/software/WebPA/
%
load('cellbonddata');
evolver=0;
showLLL= 1;
isDelFreeEdge = 1;
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

%% MAXIMUM LIKEHOOD ESTIMATION (L)
% This is by far the longest step for large tissues, so the range of mu
% tested should be as small as possible...

tic;
disp('Estimating maximum likelihood...')
murange=0.4:0.1:2.5; 
% range of mu tested. mu is typically slightly smaller than 1
% look at logL(mu), it should have a maximum !
cpt=1;
for mu=murange  
    logL(cpt)=subprog_logLikelihood_Estimation(mu,A,B,g);
    if cpt>4 && logL(cpt)<logL(cpt-1) && logL(cpt-1)<logL(cpt-2) && logL(cpt-2)<logL(cpt-3)
        break %if a clear maximum has been found, break (to save time)
    end
    mu
    cpt=cpt+1;
end
disp(['Done (' num2str(round(toc,2)) ' seconds).'])
[~,I]=max(logL);
mu_bayes=murange(I);
mu_bayes

if showLLL == 1
    figure; set(gcf,'color','w');
    plot(murange(1:min(cpt,length(murange))),logL)
    ylabel('Log Likelihood')
    xlabel('Mu')
end
%saveas(gcf,'likelihood.fig');

save muFile mu_bayes