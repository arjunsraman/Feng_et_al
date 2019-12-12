%%%%Matlab script for Feng/Raman et al%%%%
%Date: 12/10/2019

%Needed workspace

%workspace.mat contains data pertaining to manuscript
%Structures included in workspace
    %BPM: Binary Phenotype Matrix data
    %colorMap: colormap matrix for plotting Fig. 2B
    %colorMap2: colormap matrix for plotting Fig. 4
    %fitness: COPRO SEQ fitness data with delineation of rows as organisms
    %nutrients: Data pertaining to metabolic naalytes
    %RNA_seq: Data pertaining to RNA-seq
    
%% Load necessary workspace

load Workspace_FengRaman.mat
    
%% Create Fig. 2B

%Running this section by itself will create Fig. 2B
%The excel file 'Fitness_Ev1_RNA_seq_forFig2A.xlsx' shows the linear
%correlation between projection along PC1 and bacterial fitness shown in
%Fig. 2C

%Step 1: Create mcSEED enrichment matrix (shown in Fig. 2A)

%Step 1a: mcSEED enrichment matrix for S2 organisms
    
fn = fieldnames(RNA_seq.S2);
counter = 0;
for k=1:numel(fn);
    if(isnumeric(RNA_seq.S2.(fn{k})));
        counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,21) = zeros(length(RNA_seq.S2.(fn{k})(:,1)),1); %Create column of zeros to separate values from means
        RNA_seq.S2.(fn{k})(:,22) = mean(RNA_seq.S2.(fn{k})(:,1:5)'); %mcSEED average of S2 only
        RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,22); %row 1 is S2 only
        %counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,23) = mean(RNA_seq.S2.(fn{k})(:,6:10)'); %mcSEED average of S1-->S2
        %RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,23); %row 1 is S1-->S2
        counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,24) = mean(RNA_seq.S2.(fn{k})(:,11:15)'); %mcSEED average of S2-->S1
        RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,24); %row 2 is S2-->S1
        counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,25) = mean(RNA_seq.S2.(fn{k})(:,16:20)'); %mcSEED average of S1 + S2
        RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,25); %row 3 is S1 + S2
    end;
end;    

%At the end of the above 'for-loop', RNA_seq.S2.S2_mat is a matrix of the
%mcSEED metabolic pathway content for only the S2 organisms


%Step 1b: mcSEED enrichment matrix for S1 organisms

%Create mcSEED enrichment matrices for each of the S1 organisms
clear fn
fn = fieldnames(RNA_seq.S1);
counter = 0;
for k=1:numel(fn);
    if(isnumeric(RNA_seq.S1.(fn{k})));
        RNA_seq.S1.(fn{k})(:,21) = zeros(length(RNA_seq.S1.(fn{k})(:,1)),1); %Create column of zeros
        RNA_seq.S1.(fn{k})(:,22) = mean(RNA_seq.S1.(fn{k})(:,1:5)'); %mcSEED average of S1 only
        RNA_seq.S1.(fn{k})(:,23) = mean(RNA_seq.S1.(fn{k})(:,6:10)'); %mcSEED average of S1-->S2
        RNA_seq.S1.(fn{k})(:,24) = mean(RNA_seq.S1.(fn{k})(:,11:15)'); %mcSEED average of S2-->S1
        RNA_seq.S1.(fn{k})(:,25) = mean(RNA_seq.S1.(fn{k})(:,16:20)'); %mcSEED average of S1 + S2
    end;
end;
tmp = RNA_seq.S1.E_coli(:,22:25);
tmp = tmp';
tmp2 = RNA_seq.S2.S2_mat;
RNA_seq.enrichment_matrices.Ecoli = vertcat(tmp,tmp2);
clear tmp 

tmp = RNA_seq.S1.E_cass(:,22:25);
tmp = tmp';
RNA_seq.enrichment_matrices.Ecass = vertcat(tmp,tmp2);
clear tmp

tmp = RNA_seq.S1.C_bolteae(:,22:25);
tmp = tmp';
RNA_seq.enrichment_matrices.Cbolteae = vertcat(tmp,tmp2);
clear tmp

tmp = RNA_seq.S1.C_inoc(:,22:25);
tmp = tmp';
RNA_seq.enrichment_matrices.Cinoc = vertcat(tmp,tmp2);
clear tmp

%Step 1c: Add pseudocount
clear fn
fn = fieldnames(RNA_seq.enrichment_matrices);
for k=1:numel(fn);
        RNA_seq.enrichment_matrices.(fn{k}) = RNA_seq.enrichment_matrices.(fn{k}) + 0.1; %Add pseudocount
        RNA_seq.enrichment_matrices.(fn{k}) = log(RNA_seq.enrichment_matrices.(fn{k})./RNA_seq.enrichment_matrices.(fn{k})(1,:)); %Compute enrichment with S1 only being the reference
end;

%Step 1d: Create enrichment matrix shown in right panel of Fig. 2A
tmp = RNA_seq.S1.E_coli(:,23:25)';
tmp1 = RNA_seq.S1.E_cass(:,23:25)';
tmp2 = RNA_seq.S1.C_inoc(:,23:25)';
tmp3 = RNA_seq.S1.C_bolteae(:,23:25)';
tmp4 = vertcat(tmp,tmp1,tmp2,tmp3);

tmp4 = vertcat(tmp4,RNA_seq.S2.S2_mat);
tmp4 = tmp4 + 0.1; %Add pseudocount
tmp4 = log10(tmp4./tmp4(28,:)); %'tmp4' is the enrichment matrix with P distasonis S1-->S2 as reference


%Step 2: Perform eigendecomposition of mcSEED enrichment matrix
[x,y] = eig(cov(tmp4')); x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);
%x are the projections along eigenvectors and y are the eigenvectors
figure; scatter3(x(:,1),x(:,2),x(:,3),80,colorMap,'filled'); %Fig. 2B

%% Create Fig. 3B

%PC1 contains the information regarding high fitness. Performing SVD on the
%tmp4 matrix and isolating the projections corresponding to the first
%singular value will tell us which mcSEED subsystems create the PC1 axis

[U,S,V] = svd(RNA_seq.mcSEED_space.matrix);
U = U(:,1); S = S(1,1); V = V(:,1); V1 = V'; %Right singular vectors contain necessary information
new_mat = U*S*V1; %New matrix only considers eigenvector 1
figure; histogram(V1,50); %Histogram shown in Fig. 3B
[a,b] = sortrows(V); 
mcSEED_pathways_V1_ordered = RNA_seq.pathway_modules(b); %Labels shown in Fig. 3B

%% Create Fig. 3C

%Create matrix of mcSEED_enrichment with top 9 and bottom 9 (top 10% and
%bottom 10% SVD projection)

RNA_seq.mcSEED_space.matrix_trm = RNA_seq.mcSEED_space.matrix(:,[RNA_seq.SVD_analysis.V_srtd_b(1:9),RNA_seq.SVD_analysis.V_srtd_b(end-8:end)]);

%Hierarchically cluster by strain, Fig. 3C, left panel
clear Z H T out
Z = linkage(RNA_seq.mcSEED_space.matrix_trm,'average');
[H,T,out] = dendrogram(Z,36); %Dendrogram shown in Fig. 3C
RNA_seq.mcSEED_space.matrix_trm_reordered = RNA_seq.mcSEED_space.matrix_trm(out,:);
figure; imagesc(RNA_seq.mcSEED_space.matrix_trm_reordered); colormap gray; axis equal

%% Create Fig. 4A
%E faecium is a microbe that doesn't survive the introduction of S2 as well
%as the others. Why? Analyze its RNA-seq data. However, the S1-->S2
%colonization condition only contained 0.2% and as a result had a low
%RNA-seq data generation. So only look at S2-->S1 and S1 + S2 as well as S1
%and compare to the other organisms that survive by normalizing to the S1
%alone condition. 

RNA_seq.S1.E_faecium(:,16) = zeros(length(RNA_seq.S1.E_faecium(:,1)),1); %Column of zeros
RNA_seq.S1.E_faecium(:,17) = mean(RNA_seq.S1.E_faecium(:,1:5)'); %Mean of S1 only
RNA_seq.S1.E_faecium(:,18) = mean(RNA_seq.S1.E_faecium(:,6:10)'); %Mean of S2-->S1
RNA_seq.S1.E_faecium(:,19) = mean(RNA_seq.S1.E_faecium(:,11:15)'); %Mean of S1 + S2

clear tmp tmp1 tmp2 tmp3
tmp = RNA_seq.S1.E_coli(:,22:25)';
tmp1 = RNA_seq.S1.E_cass(:,22:25)';
tmp2 = RNA_seq.S1.C_inoc(:,22:25)';
tmp3 = RNA_seq.S1.C_bolteae(:,22:25)';
tmp4 = RNA_seq.S1.E_faecium(:,17:19)';

tmp5 = vertcat(tmp, tmp1, tmp2, tmp3, tmp4);
tmp5 = tmp5(:,[RNA_seq.SVD_analysis.V_srtd_b(1:9),RNA_seq.SVD_analysis.V_srtd_b(end-8:end)]);
tmp5 = tmp5 +0.1; %Pseudocount
tmp5 = log10(tmp5./tmp5(17,:)); %Normalize
clear x y
[x y] = eig(cov(tmp5')); x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);

figure; scatter3(x(:,1),x(:,2),x(:,3),80,colorMap2,'filled'); %Fig. 4A

%% Create Fig. 4B-D

[U,S,V] = svd(RNA_seq.mcSEED_space_w_Efaecium.mat);
figure; histogram(V(:,1),10);
figure; histogram(V(:,2),10);
U_3 = U(:,3); S_3 = S(3,3); V_3 = V(:,3);
figure; histogram(V_3,20);

%% Create Fig. S2A
%First compare reproducibility between experiments
%GF v GF
nutrients.expt1.data.GF_mean_std(:,1) = mean(nutrients.expt1.data.GF');
nutrients.expt1.data.GF_mean_std(:,2) = std(nutrients.expt1.data.GF');

nutrients.expt1.data.S1_mean_std(:,1) = mean(nutrients.expt1.data.S1');
nutrients.expt1.data.S1_mean_std(:,2) = std(nutrients.expt1.data.S1');

nutrients.expt1.data.S2_mean_std(:,1) = mean(nutrients.expt1.data.S2');
nutrients.expt1.data.S2_mean_std(:,2) = std(nutrients.expt1.data.S2');

nutrients.expt1.data.S1S2_mean_std(:,1) = mean(nutrients.expt1.data.S1S2');
nutrients.expt1.data.S1S2_mean_std(:,2) = std(nutrients.expt1.data.S1S2');

nutrients.expt1.data.S2S1_mean_std(:,1) = mean(nutrients.expt1.data.S2S1');
nutrients.expt1.data.S2S1_mean_std(:,2) = std(nutrients.expt1.data.S2S1');

nutrients.expt1.data.S1pS2_mean_std(:,1) = mean(nutrients.expt1.data.S1pS2');
nutrients.expt1.data.S1pS2_mean_std(:,2) = std(nutrients.expt1.data.S1pS2');

nutrients.expt2.data.GF_mean_std(:,1) = mean(nutrients.expt2.data.GF');
nutrients.expt2.data.GF_mean_std(:,2) = std(nutrients.expt2.data.GF');

nutrients.expt2.data.S1_mean_std(:,1) = mean(nutrients.expt2.data.S1');
nutrients.expt2.data.S1_mean_std(:,2) = std(nutrients.expt2.data.S1');

nutrients.expt2.data.S2_mean_std(:,1) = mean(nutrients.expt2.data.S2');
nutrients.expt2.data.S2_mean_std(:,2) = std(nutrients.expt2.data.S2');

nutrients.expt2.data.S1S2_mean_std(:,1) = mean(nutrients.expt2.data.S1S2');
nutrients.expt2.data.S1S2_mean_std(:,2) = std(nutrients.expt2.data.S1S2');

%Second add pseudocount to data 

% nutrients.expt1.data.GF_mean_std(:,3) = nutrients.expt1.data.GF_mean_std(:,1) + 0.1;
% nutrients.expt1.data.S1_mean_std(:,3) = nutrients.expt1.data.S1_mean_std(:,1) + 0.1;
% nutrients.expt1.data.S2_mean_std(:,3) = nutrients.expt1.data.S2_mean_std(:,1) + 0.1;
% nutrients.expt1.data.S1S2_mean_std(:,3) = nutrients.expt1.data.S1S2_mean_std(:,1) + 0.1;

nutrients.expt1.data.GF_pseudo = nutrients.expt1.data.GF + 0.1;
nutrients.expt1.data.S1_pseudo = nutrients.expt1.data.S1 + 0.1;
nutrients.expt1.data.S2_pseudo = nutrients.expt1.data.S2 + 0.1;
nutrients.expt1.data.S1S2_pseudo = nutrients.expt1.data.S1S2 + 0.1;
nutrients.expt1.data.S2S1_pseudo = nutrients.expt1.data.S2S1 + 0.1;
nutrients.expt1.data.S1pS2_pseudo = nutrients.expt1.data.S1pS2 + 0.1;

nutrients.expt1.data.GF_pseudo_mean = mean(nutrients.expt1.data.GF_pseudo');
nutrients.expt1.data.GF_pseudo_mean = nutrients.expt1.data.GF_pseudo_mean';

nutrients.expt1.data.S1_norm = log10(nutrients.expt1.data.S1_pseudo./nutrients.expt1.data.GF_pseudo_mean);
nutrients.expt1.data.S2_norm = log10(nutrients.expt1.data.S2_pseudo./nutrients.expt1.data.GF_pseudo_mean);
nutrients.expt1.data.S1S2_norm = log10(nutrients.expt1.data.S1S2_pseudo./nutrients.expt1.data.GF_pseudo_mean);
nutrients.expt1.data.S2S1_norm = log10(nutrients.expt1.data.S2S1_pseudo./nutrients.expt1.data.GF_pseudo_mean);
nutrients.expt1.data.S1pS2_norm = log10(nutrients.expt1.data.S1pS2_pseudo./nutrients.expt1.data.GF_pseudo_mean);

%Concatenate horizontally for all conditions to make large 56x25 matrix
nutrients.expt1.data.total_mat = horzcat(nutrients.expt1.data.S1_norm,nutrients.expt1.data.S2_norm,nutrients.expt1.data.S1S2_norm,nutrients.expt1.data.S2S1_norm,nutrients.expt1.data.S1pS2_norm);

%Hierarchically cluster by nutrients
nutrients.expt1.data.HC_analysis.Z = linkage(nutrients.expt1.data.total_mat,'average');
[H,T,out] = dendrogram(nutrients.expt1.data.HC_analysis.Z,56);
nutrients.expt1.data.HC_analysis.rows_reordered = out;
nutrients.expt1.data.HC_analysis.rows_reordered_names = nutrients.expt1.data_row_names(out);



nutrients.expt1.data.HC_analysis.Y = linkage(nutrients.expt1.data.total_mat','average');

clear H T out
[H,T,out] = dendrogram(nutrients.expt1.data.HC_analysis.Y,25);
nutrients.expt1.data.HC_analysis.columns_reordered = out;
nutrients.expt1.data.HC_analysis.columns_reordered_names = nutrients.expt1.data_column_names(out);

%visualize HC matrix
nutrients.expt1.data.total_mat_reordered = nutrients.expt1.data.total_mat(nutrients.expt1.data.HC_analysis.rows_reordered,nutrients.expt1.data.HC_analysis.columns_reordered);
figure; imagesc(nutrients.expt1.data.total_mat_reordered); axis equal; colormap gray;

%% Create S2B
%Perform eigendecomposition on nutrient matrix shown in previous section
clear x y
[x,y] = eig(cov(nutrients.expt1.data.total_mat_reordered(1:end-9,:)));
x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);
figure; scatter(x(:,1),x(:,2)); %Fig. S2B
nutrients.expt1.data.PCA_analysis.x = x;
nutrients.expt1.data.PCA_analysis.y = y;
clear x y



