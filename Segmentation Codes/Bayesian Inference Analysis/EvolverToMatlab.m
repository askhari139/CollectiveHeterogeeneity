% Readapt data shaping from Evolver.

% -> Faces can be "entire" or "on border", i.e. directly surrounding an
% entire face. Other faces are removed
% + New ranking: entire cells first and border cells at the end
% -> Edges that are not belonging to entire cells or directly connected
% are removed + new ranking : edges belonging to entire cell first
% -> Only vertices that are referred in edges matric are kept
% + new ranking

% Input arguments
% VERTICES (matrix): 1st row: vertex index
%   2nd, 3rd rows: x, y coord.
% EDGES (matrix): 1st row: edge index
%   2nd, 3rd rows: vertices indices
% FACES (array of vectors): 1st vector elt: face (cell) index //
%   other vector elts: edge indices of the face. comment: The list
%   of edge indices is vectorial and counterCW (a sign "-"
%   indicates that the edge must be read from the second vertex to
%   the first.
% EVOLVER.TENSIONS (vector): Evolver tension of every edge
% EVOLVER.PRESSURES (vector): Evolver pressure of every face
% output arguments
% C (array of vectors): cells not connected to an entire cell
%   removed
% EDGES_2 (matrix):
% VERTICES_2 (matrix):
% V_ent : number of vertices shared by at least one entire cell
% C_ent : number of entire cells
% E_tot : number of bonds (entire + border)

EVOLVER.PRESSURES= -K1*(VOLUME2-VOLUME1);
EVOLVER.TENSIONS=TENSIONS;

if exist('PERIMETER')==1 % if elastic perimeter
    for c=1:length(FACES)
        e0 = abs(FACES{c}(2:end)); % Edges belonging to cell f
        for e=e0
            EVOLVER.TENSIONS(e) = EVOLVER.TENSIONS(e) + Kper*PERIMETER(c);
        end
    end
end

EE=zeros(1,length(EDGES));
EE2=zeros(1,length(EDGES));
Fent=zeros(1,length(FACES));
VV=zeros(1,length(VERTICES));

for cptF=1:length(FACES)
    e_in_f = abs(FACES{cptF}(2:end)); % Edges belonging to cell f
    EE(e_in_f) = EE(e_in_f)+1; % EE(e) = n ->  n cells share edge e
end

for cptF=1:length(FACES)
    e_in_f =  abs(FACES{cptF}(2:end));
    if sum(EE(e_in_f)) == 2 * (length(FACES{cptF})-1)
        % if all edges of the cell belong to 2 cells
        Fent(cptF) = 1; % F(f)=1 -> cell f is considered as entire
        EE2(e_in_f) = 1; % EE2(e)=1 -> edge e belongs to an entire cell
        VV(EDGES(e_in_f, 2:3)) = 2; % VV(v)=2 -> vertex v belongs to an
        % entire cell
    end
end

for cptF=1:length(FACES)
    if Fent(cptF)==0
        % if the cell is not entire
        e_in_f = abs(FACES{cptF}(2:end));
        % edges of e_in_f belong to entire cells ?
        EE = EE2(e_in_f);
        % edges of e_in_f belong to entire cells or directly connected
        % to a such edge ?
        EE = EE + EE([2:end 1]) + EE([end 1:end-1]);
        % List of edges fullfill this last condition :
        edgesEE = EDGES(e_in_f(EE>0),:);
        % VV(v)=2 -> vertex v belongs to an entire cells
        % VV(v)=1 -> connected to an entire cell by 1 edge(border cell)
        % VV(v)=0 -> otherwise
        VV(edgesEE(:,2:3)) = max(VV(edgesEE(:,2:3)),1);
    end
end

C_ent=sum(Fent==1); %% Number of entire cells
V_ent=sum(VV>1); % Nbr of vertices belonging to entire cells
V_tot=sum(VV>0); % Nbr of vertices belonging to entire or border cells
VERTICES(VV>1,4)=[1:V_ent]; % new ranking in row 4
VERTICES(VV==1,4)=[V_ent+1:V_tot]; % new ranking in row 4
% border vertices get greater rank

%% VERTICES
VERTICES_2 = zeros(V_tot,2);
for cptV = find(VERTICES(:,4)>0) % Loop on entire and border vertices
    v = VERTICES(cptV,4); % v = new ranking
    VERTICES_2(v,1:2) = VERTICES(cptV,2:3); % positions do not change
end

%% EDGES
e=0;
for cptE=1:length(EDGES)
    w = VERTICES(EDGES(cptE,2:3),4);
    if min(w)>0 && min(w)<=V_ent
        % if at least one vertex of the edges cptE belongs to an entire cell
        e = e+1; % new ranking
        EDGES(cptE,4) = e; % new ranking in row 4
        EDGES_2(e,1:2) = VERTICES(EDGES(cptE,2:3),4)';
        % new ranking for tension vector :
        EVOLVER2.TENSIONS(e) = EVOLVER.TENSIONS(cptE);
    end
end

%% CELLS (FACES)
f=0; C = cell(1,length(FACES));
[~,ord]=sort(Fent,'descend');
for cptF = ord % loop on all cells starting by entire cells
    e_rel = FACES{cptF}(2:end);
    e = EDGES(abs(FACES{cptF}(2:end)),4); % Edges belonging to cell f
    e_rel(e==0)=[];
    e(e==0)=[];
    v=[];
    for i=1:length(e)
        if e_rel(i)>0
            v=[v; EDGES_2(e(i),1); EDGES_2(e(i),2)];
        else
            v=[v; EDGES_2(e(i),2); EDGES_2(e(i),1)];
        end
    end
    [~, I]=unique(v,'first');
    v=v(sort(I));
    v=v(end:-1:1);
    %         v1 = EDGES_2(e,1); v2 = EDGES_2(e,2);
    %         v = [v1'.*(e_rel>0) + v2'.*(e_rel<0) ];
    %         if length(v1)>0 && length(v2)>0
    %             v = [v v1(end).*(e_rel(end)<0)+v2(end).*(e_rel(end)>0)];
    %         end
    f=f+1;
    C{f} = v;
    C{f}(C{f}==0)=[];
    if isempty(C{f})
        % cell f removed if the cell does not contain at least one vertex
        % belonging to an entire cell
        C(f)=[];
        f = f-1;
    else
        %         elseif f <= C_ent
        % new ranking for pressure vector
        EVOLVER2.PRESSURES(f) = EVOLVER.PRESSURES(cptF);
    end
end

V = VERTICES_2;
E = EDGES_2;
E_tot = length(E);
C_tot = length(C);

EVOLVER.PRESSURES=EVOLVER2.PRESSURES/mean(EVOLVER2.TENSIONS);
EVOLVER.TENSIONS=EVOLVER2.TENSIONS/mean(EVOLVER2.TENSIONS); % normalization (mean tension is not equal 1 because of perimeter elasticity)