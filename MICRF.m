function [out,w]=MICRF(nodefile,netfile,outputfile,pi0)
%%
% Inputs:
% nodefile: node score file with two columns: Gene and score, tab separated.
%           When this is the only input, genes should be listed with GeneName. (Required)
% netfile:  network speific input by users with three clolumns gene1 gene2
%           and betweenness; one line with one edge. Or a specific string to select
%           one network used in MICRF paper, which should be one of the following
%           strings: 'STRING', 'iRefIndex', 'HPRD', 'CORR', 'CoEXP', 'CoPrePPI'. 
%           Default is 'CoPrePPI'. (Optional)
% outputfile:   output file name; if not given, it will only return a cell 
%               variable. (Optional)
% pi0:  the prior fraction of non-risk genes, default is 0.94. (Optional)
%
% Outputs:
% out:  a cell variable the same with output file
% w:    the optimal parameters used in MICRF model

%% initial inputs
netstrs={'STRING', 'iRefIndex', 'HPRD', 'CORR', 'CoEXP', 'CoPrePPI'};
netpath='data/network/';
netfiles={'STRINGnetmap.txt', 'iRefIndexm.txt', 'StringNew_HPRD_mnet.txt', 'brainspan_net_cor.txt', 'brainspan_net_top5.txt', 'ComCo_PrePPI.txt'};
bepath='data/Network_betweenness/';
benesss={'Betweenness_edge_STRING.txt', 'Betweenness_edge_iRef.txt', 'Betweenness_edge_HPRD.txt', 'Betweenness_edge_corr1.txt', 'Betweenness_edge_coexp5.txt', 'Betweenness_edge_Co_PrePPI.txt'};
Ws=zeros(2,6);
Ws(1,:)=[6.8557,6.1026,7.3602,7.4168,3.9194,3.2318];
Ws(2,:)=[10.653,9.4211,9.0906,8.0247,13.117,12.861];
flag=0; % 0 read MICRF files and 1 read users' input files

if nargin == 1
    j=6;
    netfile=[netpath netfiles{j}];
    beness=[bepath benesss{j}];
    w0=Ws(:,j);
end

if nargin >= 2
    [~,j] = ismember(netfile,netstrs);
    if j > 0
        netfile=[netpath netfiles{j}];
        beness=[bepath benesss{j}];
        w0=Ws(:,j);
    elseif j==0
        flag = 1;
        w0=[mean(Ws(1,:));mean(Ws(2,:))];
    end
end

if nargin <= 2
    outputfile='';
end

if nargin < 4
    pi0=0.94;
end

%% read related files
addpath(genpath(pwd)) % use functions in UGM package
% read the network edges and betweenness
if flag == 1
    [node1,node2,be] = textread(netfile,'%s%s%f','delimiter','\t');
    enode1=node1;enode2=node2;
elseif flag == 0
    [node1,node2,~] = textread(netfile,'%s%s%f','delimiter','\t');
    [enode1,enode2,be] = textread(beness,'%s -- %s\t%f');
end
genes = union(node1,node2);
nStates = ones(1,length(genes)) * 2;
% read input node score
fcon =  fopen(nodefile,'r');
C = textscan(fcon,'%s%s','delimiter','\t');
fclose(fcon);

%% step 1: network adjacency matrix 
adj = step1_adj(node1,node2);
% Make structure that tracks edge information
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

%% step 2 and 3: compute initial edge and node features; MICRF training, decoding and inferring
[Xnode,Xedge,nodeMap,edgeMap]=step2_feature(genes,C,enode1,enode2,be,edgeStruct,pi0);
[Y,nps,w]=step3_MICRF(Xnode,Xedge,nodeMap,edgeMap,edgeStruct,w0);

%% step 4: output files
out=step4_output(genes,Y,nps,Xnode,outputfile);

end

function adj = step1_adj(node1,node2)
    genes = union(node1,node2);
    nNodes = length(genes);
    % Make Adjacency Matrix 
    adj = zeros(nNodes); % Symmetric {0,1} matrix containing edges
    [~,i] = ismember(node1,genes);
    [~,j] = ismember(node2,genes);
    sub1 = nNodes*(j-1) + i;
    sub2 = nNodes*(i-1) + j; 
    adj(sub1)=1;
    adj(sub2)=1;  
    adj(1:(nNodes+1):end)=0;
end

function [Xnode,Xedge,nodeMap,edgeMap]=step2_feature(genes,C,enode1,enode2,be,edgeStruct,pi0)
%% prepare to node feature
% TADA score
nodes =C{1};
score1=C{2};
% delete NA value in score1
subs=find(~strcmp(score1,'NA'));
nodes=nodes(subs);
score1=score1(subs);
score1=str2double(score1);
node_s = intersect(nodes,genes);

%% prepare to edge feature
F=log(be)/max(log(be)); % max: normal distritbution
edgesb1 = strcat(enode1,'_',enode2);
edgesb2 = strcat(enode2,'_',enode1);

%% initial node and edge features: % 1 risk state, 2 non-risk state
% risk genes prior information
Pri=zeros(2,1); 
Pri(1)=1-pi0;
Pri(2)=pi0;

nInstance = 1;
nNodes = length(genes);
nEdges = size(edgeStruct.edgeEnds,1);
nEdgeFeatures=4;

% Make node features
nNodeFeatures = 2;
Xnode = zeros(nInstance, nNodeFeatures, nNodes);
Xnode(1,1,:) = Pri(1)*Pri(2);
Xnode(1,2,:) = Pri(2)*Pri(1);
[~,j] = ismember(node_s,genes);
[~,i] = ismember(node_s,nodes);        
Xnode(1,1,j) = score1(i)*Pri(2);
Xnode(1,2,j) = (1 - score1(i))*Pri(1);
Xnode = Xnode / mean(Xnode(:));

% Make edge features
n1 = edgeStruct.edgeEnds(:,1);
n2 = edgeStruct.edgeEnds(:,2);
edges = strcat(cellstr(genes(n1)),'_',cellstr(genes(n2)));
[~,ind1] = ismember(edges,edgesb1);
[~,ind2] = ismember(edges,edgesb2);
sube = ind1 + ind2;
% Compute Edge Features (use node features from both nodes)
Xedge = zeros(nInstance,nEdgeFeatures,nEdges);
Xedge(1,1,:) = reshape(min(Xnode(1,1,n1), Xnode(1,1,n2)),[],1).*F(sube);
Xedge(1,2,:) = reshape(min(Xnode(1,1,n1), Xnode(1,2,n2)),[],1).*F(sube);
Xedge(1,3,:) = reshape(min(Xnode(1,2,n1), Xnode(1,1,n2)),[],1).*F(sube);
Xedge(1,4,:) = reshape(max(Xnode(1,2,n1), Xnode(1,2,n2)),[],1).*F(sube);
[nodeMap,edgeMap] = UGM_makeCRFmaps_ehq(Xnode,Xedge,edgeStruct);

end

function [Y,nps,w]=step3_MICRF(Xnode,Xedge,nodeMap,edgeMap,edgeStruct,w0)

%% training 
% initial parameters
f = 100000;
f0 = f;
w=w0;
flag = 0;
iter = 0;
maxiter=100;
inferFunc = @UGM_Infer_LBP; 

maxFunEvals = 20;
options = [];
options.maxFunEvals = maxFunEvals;
options.progTol = 1e-3;
options.optTol = 1e-3;
        
while (flag==0 && iter <= maxiter)
    % update potentials
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
    % decoding a CRF model
    Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    Y = int32(Y');

    % training a CRF model
    lambda = ones(size(w)); 
    regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);
    [w,f]= minFunc(regFunObj,w,options);
   
    % fix bugs for illegal direction
    if isnan(f)==1 
        %fprintf('%d Here\n',iter);
        w(1)=randi(20);
        w(2)=randi(20);
        [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
        Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
        Y = int32(Y');
        lambda = ones(size(w)); %lambda(2) = 10;
        regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc);
        [w,f]= minFunc(regFunObj,w,options);
    end
    
    if isnan(f)
        w = w0;
        break;
    end
    if norm(w-w0,1) <= 1e-5 && norm(f-f0,1) <= 1e-5
        flag = 1;
    else
        f0 = f;
        w0 = w;
    end
    
    iter = iter + 1;
    %fprintf('%d\n',iter);
end

%% decoding and inferring    
[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,1);
Y = UGM_Decode_LBP(nodePot,edgePot,edgeStruct); % decoding a CRF model
Y = int32(Y');
nps = UGM_Infer_LBP(nodePot,edgePot,edgeStruct); % infer conditional probability
if any(isnan(nps))
    % nps = UGM_Infer_TRBP(nodePot,edgePot,edgeStruct); % trying TRBP
    nps(isnan(nps)) = Xnode(isnan(nps)); 
end

end

function out=step4_output(genes,Y,nps,Xnode,outputfile)
% step 4: output files
nodescore = [nps(:,1),reshape(Xnode(1,1,:),[],1)];
% nps with size nNodes * nstate
[~,Ind]=sortrows(-1*nodescore); % 1 as risk, 2 as non-risk
ng = length(genes);
out=cell(ng,5);
for i = 1:ng
    out(i,:)={genes{Ind(i)},nps(Ind(i),1),nps(Ind(i),2),Xnode(1,1,Ind(i)),Y(Ind(i))};
end

if strcmp(outputfile,'') == 0
    fileID = fopen(outputfile,'w');
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\n','Gene','Risk-state','Non-risk-state','Node-feature','Label');
    for i = 1:length(genes)
        fprintf(fileID,'%s\t%12.8f\t%12.8f\t%12.8f\t%d\n',genes{Ind(i)},nps(Ind(i),1),nps(Ind(i),2),Xnode(1,1,Ind(i)),Y(Ind(i)));
    end
    fclose(fileID);
end

end

function [nodeMap,edgeMap] = UGM_makeCRFmaps_ehq(Xnode,Xedge,edgeStruct)
% Assumes that all nodes have the same number of states
nNodes = size(Xnode,3);
nEdges = edgeStruct.nEdges;
nStates = edgeStruct.nStates;
maxState = max(nStates);

UGM_assert(min(nStates)==maxState,'UGM_makeCRFMaps assumes that all nodes must have the same number of states');
nStates = nStates(1);
nNodeFeatures = size(Xnode,2);
nEdgeFeatures = size(Xedge,2);

nodeMap = zeros(nNodes,nStates,nNodeFeatures,'int32');
fNum=1;
nodeMap(:,1,1) = fNum;
nodeMap(:,1,2) = 0;
nodeMap(:,2,1) = 0;
nodeMap(:,2,2) = fNum;
nNodeParams = max(nodeMap(:));

fNum=fNum+1;
edgeMap = zeros(nStates,nStates,nEdges,nEdgeFeatures,'int32');
edgeMap(1,1,:,1) = fNum;
edgeMap(1,2,:,2) = fNum;
edgeMap(2,1,:,3) = fNum;
edgeMap(2,2,:,4) = fNum;

end