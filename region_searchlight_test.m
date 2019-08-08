% Language fMRI Project - Creating matrix of rhos for each voxel using
% searchlights
% Toolbox for RSA Paper - implementation of techniques
% 6/18/19

% All language areas, with 10 subjects and no bootstrap p values

% Print: each voxel being analyzed (in which subject)

addpath('/ncf/gershman/Lab/scripts/matlab/Language-fMRI/xlarge-spanish-2layer-rnn/avg');
addpath('/ncf/gershman/Lab/scripts/spm12');

fprintf('test 1');

clear rhos h p ci stats

%% Configuration

model = 'xlarge-english-to-spanish-model-pred-layer2-avg'; % Model you are testing against 

region_name = 'LAngG'; % Region you are analyzing 
% region_name = 'Hippocampus_L'; 
% region_name = 'Insula_L'; 
% region_name = 'LMidPostTemp';
% region_name = 'LPostTemp';
% region_name = 'LMidAntTemp';
% region_name = 'LIFG';
% region_name = 'LAntTemp';
% region_name = 'LIFGorb'; 
% region_name = 'LMFG'; 

% Some choices: 'LAngG', 'Hippocampus_L', 'Hippocampus_R', 'Insula_L', 'Insula_R',
% 'LMidPostTemp', 'LPostTemp', 'LMidAntTemp', 'LIFG', 'LAntTemp', 'LIFGorb', 'LAngG', 'LMFG'}
datadir = '/ncf/gershman/Lab/scripts/matlab/Language-fMRI/examplesGLM'; % Where the fMRI data lives on your computer 

%% Load model data into feature matrix and initializing variables
load(model); 

fprintf('test 2');

n = 240     % Number of sentences
b = 10    % Number of bootstrap samples

for i=1:n
    feature(i,:)=eval(['sentence', num2str(i), ';']); 
end

fprintf('test 3');


%% Get a brain area
% TODO make this better. Currently assumes your region_name is correct, you
% should probably check the printed idx 
labels = load([datadir '/subj1/examplesGLM.mat'], 'labels_aal', 'labels_langloc');
if isempty(find(strcmp(labels.labels_langloc, region_name)))
    isLangArea = false;
    useLabels = labels.labels_aal;
elseif isempty(find(strcmp(labels.labels_aal, region_name)))
    isLangArea = true;
    useLabels = labels.labels_langloc;
end 

idx = find(strcmp(useLabels, region_name)) % Map region name to number 

%% Loop for each voxel

% Compute model and neural RDMs

load([datadir '/subj1/examplesGLM.mat'], 'sentencesPresent', 'volmask', 'vollangloc', 'volaal', 'examples_sentences');

%{
for i = 1:b
    scrambled_model_RDMs{i} = pdist(feature(randperm(n), :), 'cosine');
end
%}

%cell array

%%

rhos = zeros(size(volmask));
pvalues = zeros(size(volmask));
bstrp_pvalues = zeros(size(volmask));
realmask = logical(ones(size(volmask)));
all_sentencesPresent = logical(ones(n,1));


%% creating realmask
for s= [1 2 4 5 7 8 9 10 11]
    fprintf('test 4.%i\n',s);
    % Create region mask 
    if isLangArea
        load([datadir '/subj' num2str(s) '/examplesGLM.mat'], 'vollangloc','sentencesPresent','examples_sentences','volmask');
        realmask = realmask & (logical(vollangloc)); 
    else
        load([datadir '/subj' num2str(s) '/examplesGLM.mat'], 'volaal','sentencesPresent','examples_sentences','volmask');
        realmask = realmask & (volaal == idx);
    end
    all_sentencesPresent = all_sentencesPresent & sentencesPresent;
    volmasks{s} = volmask;
    brain_activity{s} = examples_sentences;
end

model_RDM = pdist(feature(logical(all_sentencesPresent), :), 'cosine');

%bstrp_samples = zeros(b,1);

%% 
for i=1:size(volmask,1)
    for j = 1:size(volmask,2)
        for k = 1:size(volmask,3)
            for s=[1 2 4 5 7 8 9 10 11] % Subjects. All: [1 2 3 4 5 7 8 9 10 11]; working: [1 2 4 5 7 8 9 10 11]
                isrelevant = true;
                if ~realmask(i,j,k)
                    isrelevant = false;
                    continue;
                end
                fprintf('%i:%i,%i,%i\n', s,i,j,k)
                
                % TODO: for each voxel 
                % TODO: look at create_spherical_mask_helper to create a searchlight mask (currently just using the region_mask)


                sphere_coords = [];

                sphere_mask = zeros(size(realmask));

                x = i; y = j; z = k; r = 4; min_x = 1; max_x = size(volmask,1); min_y = 1; max_y = size(volmask,2); min_z = 1; max_z = size(volmask,3);

                for newx = floor(x - r) : ceil(x + r)
                    if newx < min_x || newx > max_x, continue; end
                    for newy = floor(y - r) : ceil(y + r)
                        if newy < min_y || newy > max_y, continue; end
                        for newz = floor(z - r) : ceil(z + r)
                            if newz < min_z || newz > max_z, continue; end
                            if ~realmask(newx, newy, newz), continue; end
                            if (x - newx)^2 + (y - newy)^2 + (z - newz)^2 > r^2, continue; end
                            sphere_mask(newx, newy, newz) = 1;
                           % mni = cor2mni([newx newy newz], Vmask.mat);
                           % sphere_coords = [sphere_coords; mni];
                            sphere_coords = [sphere_coords; newx newy newz];
                        end
                    end
                end

                sphere_mask = logical(sphere_mask);

                sphere_mask = sphere_mask(volmasks{s}); % Convert 3D -> 1D mask (volmask is the whole brain)
                badvoxels = any(isnan(brain_activity{s}), 1);
                sphere_mask = sphere_mask & (~badvoxels');

                % Compute neural RDM for searchlight 
                neural_RDM = pdist(brain_activity{s}(:, sphere_mask)); % average RDMs from subjects (2D matrix):

                %cell array (see above)
                if s == 1
                    allRDMs = neural_RDM;
                else
                    allRDMs = cat(3,allRDMs,neural_RDM);
                end

                %{
                generate bootstrap samples
                model_RDM = pdist(feature(randperm(n), :), 'cosine');
                model_RDM = pdist(features((permute numbers 1-n), sphere_mask)); % average RDMs from subjects (2D matrix):
                    generate b
                    compute spearman rank coefficient with avg neural RDM; store in
                    array
                    calculate percentile of actual rho
                (do FDR over all voxels)

                %}

                % Do RSA between this neural RDM and the model RDM, and save stats 
                % TODO: edit rhos format to support the searchlight regions! 
            end
            
            if isrelevant
                avg_RDM = mean(allRDMs, 3);
                [rho, p] = corr(model_RDM', avg_RDM', 'type', 'spearman'); % Need to be cols, so transpose' 
                rhos(i,j,k) = rho;
                pvalues(i,j,k) = p;
                
                %{
                for l = 1:b
                    [rho, p] = corr(scrambled_model_RDMs{l}', avg_RDM', 'type', 'spearman');
                    bstrp_samples(l) = rho;
                end
                bstrp_samples = sort(bstrp_samples);
                count = 0;
                for l = 1:b %use variable for b
                    count = l;
                    if rhos(i,j,k) > bstrp_samples(l)
                        break;
                    end
                end

                bstrp_pvalues(i,j,k) = 1 - count/b;
                
                %}
            end
            
        end
    end
end

fprintf('test 5\n');


realmask = logical(realmask);
rhos = realmask.* rhos;
pvalues = realmask.*pvalues;
%bstrp_pvalues = realmask.*bstrp_pvalues

%%
rhos_threshold = rhos;
rhos_threshold(pvalues > .00001) = NaN;

fprintf('test 6\n');
ccnl_view_mask(rhos_threshold)

%create mask where pvalue > 0.05

%%
save ('region_searchlight_test1.mat','-v7.3')
%{
%% FDR
fdr = mafdr(bstrp_pvalues);
for i=1:size(fdr,1)
    for j = 1:size(fdr,2)
        for k = 1:size(fdr,3)
            if fdr(i,j,k) > .05
                fdr(i,j,k) = NaN;
            end
        end
    end
end

%}



% order? FDR before bootstrapping? 

%resize rhos for subj3? 

%% Compute t-test for rhos
% TODO: make sure still computing the correct stats if you edit rhos above


%{
rhos(6) = [] % delete col for subject 6 (doesnt exist)
rhos = atanh(rhos) % Fisher transform 
[h, p, ci, stats] = ttest(rhos);
fprintf('average rho=%f, t(%d)=%f, p=%f\n', mean(rhos), stats.df, stats.tstat, p); 
newstring = sprintf('average rho=%f, t(%d)=%f, p=%f\n', mean(rhos), stats.df, stats.tstat, p); 

% TODO create null distribution for actual comparisons 


rhos = [.3, .7, ...]
r = atanh(rhos)
[h,p,ci,stat] = ttest(r)
t = stat.tstat
%}