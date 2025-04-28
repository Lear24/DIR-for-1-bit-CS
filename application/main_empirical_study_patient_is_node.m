% The program for empirical study.
% The data contains expriments for 15 participants with 3 sessions and 24
% trails in a 62-channel system.
%
% We consider the first 15 trails for each participants and each session
% There 3 method to ensure   e nodes.
%
% 1. 15 participants with one channels (15 nodes)
% 2. 15 participants and 3 sessions with one channels (45 nodes)
% 3. 15 participants, 3 sessions and 15 trails with one channels (675 nodes)
%
%


%% 15 participants with one channels (15 nodes)
clear
clc

m = 15; d = 200;T1 = 1;T2 = 2000;Node_i = 1;Trail_i = 1;

es_data = cell(m,1);
es_beta = zeros(m,d);
es_name = ["ww";"wsf";"wyw";"xyl";"ys";"zjy";"djc";"jl";"jj";"lqj";"ly";"mhw";"phl";"sxy";"wk";...
    ];
clearvars  -except n d T1 T2 Node_i Trail_i es*


%%
currentFolder = pwd;
File = dir(fullfile(currentFolder,'*.mat'));
% File = orderfields(File);
% es_file_name = File(3*6+1).name;

%%
m_i = 1;
es_file_name = File(3*(m_i-1)+1).name;
load(es_file_name);
str_var = [es_name(m_i)+"_eeg"+num2str(Trail_i)];
[channal_num, sample_freaqency] = size(eval(str_var));
es_data(m_i,1) = mat2cell(eval(str_var),channal_num, sample_freaqency);
clearvars  -except n d T1 T2 Node_i Trail_i File es*

for m_i = 2:15
    es_file_name = File(3*(m_i-1)+1).name;
    load(es_file_name);
    str_var = [es_name(m_i)+"_eeg"+num2str(Trail_i)];
    [channal_num, sample_freaqency] = size(eval(str_var));
    es_data(m_i,1) = mat2cell(eval(str_var),channal_num, sample_freaqency);
    clearvars  -except n d T1 T2 Node_i Trail_i File es*
    
    
end

%%
save('expirical_stady_node15.mat','es_data')

%% 45 participants with one channels (15 nodes)
clear
clc

m = 15; d = 200;T1 = 1;T2 = 2000;Node_i = 1;Trail_i = 1;

es_data = cell(m,3);
es_beta = zeros(m,d);
es_name = ["ww";"wsf";"wyw";"xyl";"ys";"zjy";"djc";"jl";"jj";"lqj";"ly";"mhw";"phl";"sxy";"wk";...
    ];
clearvars  -except n d T1 T2 Node_i Trail_i es*



currentFolder = pwd;
File = dir(fullfile(currentFolder,'*.mat'));
% File = orderfields(File);
% es_file_name = File(3*6+1).name;


m_i = 1;
for i = 1:3
    es_file_name = File(3*(m_i-1)+i).name;
    load(es_file_name);
    str_var = [es_name(m_i)+"_eeg"+num2str(Trail_i)];
    [channal_num, sample_freaqency] = size(eval(str_var));
    es_data(m_i,i) = mat2cell(eval(str_var),channal_num, sample_freaqency);
    clearvars  -except m_i n d T1 T2 Node_i Trail_i File es*
end

for m_i = 2:15
    
    for i = 1:3
    es_file_name = File(3*(m_i-1)+1).name;
    load(es_file_name);
    str_var = [es_name(m_i)+"_eeg"+num2str(Trail_i)];
    [channal_num, sample_freaqency] = size(eval(str_var));
    es_data(m_i,i) = mat2cell(eval(str_var),channal_num, sample_freaqency);
    clearvars  -except m_i n d T1 T2 Node_i Trail_i File es*
    end
    
end

%%
save('expirical_stady_node45.mat','es_data')

