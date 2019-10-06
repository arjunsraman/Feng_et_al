%%%%%Matlab Script for Feng et al%%%
%Date: 09/30/2019

%%Needed workspace

%workspace.mat contains data pertaining to manuscript
%Structures included in workspace
    %BPM: Binary Phenotype Matrix data
    %fitness: COPRO SEQ fitness data with delineation of rows as organisms
    %nutrients: Data pertaining to metabolic naalytes
    %RNA_seq: Data pertaining to RNA-seq



%% RNA_seq analysis

%Creation of matrices for C. bolteae, C. inocuum, E.coli, E. casseliflavus

%Make the S2 matrix

%For all files in the RNA_seq.S2 structure, need to calculate the mean and
%standard deviation for the mcSEED pathway modules in increments of 5
%corresponding to the different colonization conditions

%Save the mean for a particular organism in a particular colonization
%condition as the row of the S2_RNAseq_count matrix

%Compute Mean and create matrix
fn = fieldnames(RNA_seq.S2);
counter = 0
for k=1:numel(fn);
    if(isnumeric(RNA_seq.S2.(fn{k})));
        counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,21) = zeros(length(RNA_seq.S2.(fn{k})(:,1)),1); %Create column of zeros to separate values from means
        RNA_seq.S2.(fn{k})(:,22) = mean(RNA_seq.S2.(fn{k})(:,1:5)'); %mcSEED average of S2 only
        %RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,22); %row 1 is S2 only
        %counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,23) = mean(RNA_seq.S2.(fn{k})(:,6:10)'); %mcSEED average of S1-->S2
        RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,23); %row 1 is S1-->S2
        counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,24) = mean(RNA_seq.S2.(fn{k})(:,11:15)'); %mcSEED average of S2-->S1
        RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,24); %row 2 is S2-->S1
        counter = counter + 1;
        RNA_seq.S2.(fn{k})(:,25) = mean(RNA_seq.S2.(fn{k})(:,16:20)'); %mcSEED average of S1 + S2
        RNA_seq.S2.S2_mat(counter,:) = RNA_seq.S2.(fn{k})(:,25); %row 3 is S1 + S2
    end;
end;

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

%Create mcSEED enrichment matrices for each of the S1 organisms
%Ecoli
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

%Add pseudocount to S1 matrices and compute enrichment matrices
clear fn
fn = fieldnames(RNA_seq.enrichment_matrices);
for k=1:numel(fn);
        RNA_seq.enrichment_matrices.(fn{k}) = RNA_seq.enrichment_matrices.(fn{k}) + 0.1; %Add pseudocount
        RNA_seq.enrichment_matrices.(fn{k}) = log(RNA_seq.enrichment_matrices.(fn{k})./RNA_seq.enrichment_matrices.(fn{k})(1,:)); %Compute enrichment with S1 only being the reference
end;

%Ecoli
clear x y
[x,y] = eig(cov(RNA_seq.enrichment_matrices.Ecoli'));
x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);

figure; scatter3(x(:,1),x(:,2),x(:,3));

%Ecass
clear x y
[x,y] = eig(cov(RNA_seq.enrichment_matrices.Ecass'));
x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);
figure; scatter3(x(:,1),x(:,2),x(:,3));
%Cinocuum
clear x y
[x,y] = eig(cov(RNA_seq.enrichment_matrices.Cinoc'));
x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);
figure; scatter3(x(:,1),x(:,2),x(:,3));

clear x y
[x,y] = eig(cov(RNA_seq.enrichment_matrices.Cbolteae'));
x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);
figure; scatter3(x(:,1),x(:,2),x(:,3));


%% Create a single enrichment matrix with P. distasonis S1 + S2 as the reference
tmp = RNA_seq.S1.E_coli(:,23:25)';
tmp1 = RNA_seq.S1.E_cass(:,23:25)';
tmp2 = RNA_seq.S1.C_inoc(:,23:25)';
tmp3 = RNA_seq.S1.C_bolteae(:,23:25)';
tmp4 = vertcat(tmp,tmp1,tmp2,tmp3);

tmp4 = vertcat(tmp4,RNA_seq.S2.S2_mat);
tmp4 = tmp4 + 0.1; %Add pseudocount
tmp4 = log10(tmp4./tmp4(28,:)); %P distasonis S1-->S2 as reference

%Perform eigendecomposition
[x,y] = eig(cov(tmp4')); x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);
figure; scatter3(x(:,1),x(:,2),x(:,3));


%Now plot fitness as a function of projection along Ev1
%Did this in excel. This shows convincingly that Ev1 contains information
%regarding the fitness of 

%Color the 3-D space above by S1-->S2, S2-->S1, S1 + S2
colorMap = [];
%Created in excel file (Fitness_Ev1_RNA_seq_forFig2.xlsx)

figure; scatter3(x(:,1),x(:,2),x(:,3),80,colorMap,'filled');

%save variables
RNA_seq.mcSEED_space.x = x;
RNA_seq.mcSEED_space.y = y;
RNA_seq.mcSEED_space.matrix = tmp4;

%% Extracting features that explain high fitness

%PC1 contains the information regarding high fitness. Performing SVD on the
%tmp4 matrix and isolating the projections corresponding to the first
%singular value will tell us which mcSEED subsystems create the PC1 axis

[U,S,V] = svd(RNA_seq.mcSEED_space.matrix);
U = U(:,1); S = S(1,1); V = V(:,1); V1 = V'; %Right singular vectors contain necessary information
new_mat = U*S*V1; %New matrix only considers eigenvector 1
figure; histogram(V1,50);
[a,b] = sortrows(V);
mcSEED_pathways_V1_ordered = RNA_seq.pathway_modules(b);
%Save variables
RNA_seq.SVD_analysis.U = U;
RNA_seq.SVD_analysis. V = V;
RNA_seq.SVD_analysis.V1 = V1;
RNA_seq.SVD_analysis.S = S;
RNA_seq.SVD_analysis.new_mat = new_mat;
RNA_seq.SVD_analysis.V_srtd_a = a;
RNA_seq.SVD_analysis.V_srtd_b = b;
RNA_seq.SVD_analysis.mcSEED_V_srtd = mcSEED_pathways_V1_ordered;
clear U V V1 S new_mat a b mcSEED_pathways_V1_ordered

%Create matrix of mcSEED_enrichment with top 9 and bottom 9 (top 10% and
%bottom 10% SVD projection)

RNA_seq.mcSEED_space.matrix_trm = RNA_seq.mcSEED_space.matrix(:,[RNA_seq.SVD_analysis.V_srtd_b(1:9),RNA_seq.SVD_analysis.V_srtd_b(end-8:end)]);
figure; imagesc(RNA_seq.mcSEED_space.matrix_trm); colormap gray; axis equal; 

%Reorder such that P. distasonis is on the bottom of the graph (did by
%hand)
figure; imagesc(RNA_seq.mcSEED_space.matrix_trm); colormap gray; axis equal; 


RNA_seq.mcSEED_space.matrix_row_names_reordered = RNA_seq.mcSEED_space.matrix_row_names(out);

%Hierarchically cluster by strain
clear Z H T out
Z = linkage(RNA_seq.mcSEED_space.matrix_trm,'average');
[H,T,out] = dendrogram(Z,36);
RNA_seq.mcSEED_space.matrix_trm_reordered = RNA_seq.mcSEED_space.matrix_trm(out,:);
figure; imagesc(RNA_seq.mcSEED_space.matrix_trm_reordered); colormap gray; axis equal



%% Comparison with E. faecium
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
%Normalize to E_faecium S1 only
tmp5 = tmp5 +0.1 %pseudocount
tmp5 = log10(tmp5./tmp5(17,:)); %Normalize
clear x y
[x,y] = eig(cov(tmp5')); x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);

figure; histogram(y,80);

figure; scatter3(x(:,1),x(:,2),x(:,3));
%So you can see that there's a separation between the transcriptomes of E.
%faecium and the others. What if you make the space only using the mcSEED
%pathway modules that separate the species by fitness in Fig. 2D?

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

figure; histogram(y,80);

figure; scatter3(x(:,1),x(:,2),x(:,3));

%Wow! This means that the S1 organisms that live are pre-adapted whereas E.
%faecium is actively adapting.!!!

%Need to color the space by S1 only, S1-->S2, S2-->S1, S1 + S2
%Do this in excel by hand. 
colorMap2 = [];

figure; scatter3(x(:,1),x(:,2),x(:,3),80,colorMap2,'filled');

%Save these variables
RNA_seq.mcSEED_space_w_Efaecium.x = x;
RNA_seq.mcSEED_space_w_Efaecium.y = y;
RNA_seq.mcSEED_space_w_Efaecium.mat = tmp5;

clear x y tmp tmp1 tmp2 tmp3 tmp4 tmp5


%Eigenvector 3 holds the adaptive information regarding E. faecium

[U,S,V] = svd(RNA_seq.mcSEED_space_w_Efaecium.mat);
%Isolate eigenvector 3
U_3 = U(:,3); S_3 = S(3,3); V_3 = V(:,3);

figure; histogram(V_3,20);

%Why can E. faecium mostly adapt through ev3 and not the other
%eigenvectors? One hypothesis is that it is missing the genetic ability to
%utilize the other pathways. 

%Q1: What is the histogram of projections along 1st right singular vector?
figure; histogram(V(:,1),10);
figure; histogram(V(:,2),10);

%Panels E-G are from excel file: mcSEED BPM subsystem panels E-G.xlsx
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.E_faecium_v_Ecoli_enrichment); colormap gray;



%Compare the BPM of E.coli and E-faecium for the features that project onto
%SV1

figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.E_faecium_v_Ecoli_BPM); colormap gray; axis equal


%E_faecium_v_ecoli_enrichment_column_names(1,:) is for enrichment and (2,:)
%is for BPM


%Now compare E_faecium with E_cass and Cinoc along SVD2

figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.E_faecium_v_Ecass_Cinoc_enrichment); colormap gray; axis equal
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.E_faecium_v_E_cass_Cinoc_BPM); colormap gray; axis equal; 


%Now compare E_faecium with C.bolteae
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.E_faecium_v_Cbolteae_enrichment); colormap gray; axis equal;

%% Comparison with E.faecium part II
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.mat_trm_ev1); colormap gray; axis equal
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.mat_trm_ev2); colormap gray; axis equal
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.mat_trm_ev3); colormap gray; axis equal

%Create genetic barrier matrices (made by hand in excel: 'Genetic barrier
%analysis.xlsx' in the /Main_Figures/Fig_3 folder

figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.genetic_barrier.ev1_mat); colormap gray; axis equal
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.genetic_barrier.ev2_mat); colormap gray; axis equal
figure; imagesc(RNA_seq.mcSEED_space_w_Efaecium.genetic_barrier.ev3_mat); colormap gray; axis equal
%% NUTRIENT LANDSCAPE
%Row names are saved
%Columns are mice

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
nutrients.expt1.data_column_names = {};
nutrients.expt1.data.HC_analysis.columns_reordered_names = nutrients.expt1.data_column_names(out);

%visualize HC matrix
nutrients.expt1.data.total_mat_reordered = nutrients.expt1.data.total_mat(nutrients.expt1.data.HC_analysis.rows_reordered,nutrients.expt1.data.HC_analysis.columns_reordered);
figure; imagesc(nutrients.expt1.data.total_mat_reordered); axis equal; colormap gray;

%Perform PCA on this matrix 
[x,y] = eig(cov(nutrients.expt1.data.total_mat_reordered(1:end-9,:)));
x = fliplr(x); y = diag(y); y = flipud(y); y = y/sum(y);
figure; scatter(x(:,1),x(:,2));
nutrients.expt1.data.PCA_analysis.x = x;
nutrients.expt1.data.PCA_analysis.y = y;
clear x y


%Perform SVD on the matrix to isolate eigenvector 2
[U,S,V] = svd(nutrients.expt1.data.total_mat_reordered(1:end-9,:));
figure; histogram(U(:,2),20);
nutrients.expt1.data.SVD_analysis.U = U;
nutrients.expt1.data.SVD_analysis.V = V;
nutrients.expt1.data.SVD_analysis.S = S;
clear U S V


%The metabolites that project most along eigenvector 2 are
%N-acetylgalactosamine, Proline, and xylitol. Isolate these and show
%heatmap for figure

nutrients.expt1.data.SVD_analysis.SVD2_important_rows = [24,45,47];
nutrients.expt1.data.SVD_analysis.SVD2_mat = nutrients.expt1.data.total_mat_reordered(nutrients.expt1.data.SVD_analysis.SVD2_important_rows,:);

figure; imagesc(nutrients.expt1.data.SVD_analysis.SVD2_mat); colormap gray; axis equal;
