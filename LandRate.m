%     LandRate toolbox: Eye movement analysis and landscape rating tool
%     Copyright (C) 2017 Vassilios Krassanakis 
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%     For further information, please email me at [krasanakisv@gmail.com], [krasvas@mail.ntua.gr], or [krasvas@teiath.gr] 

clear all
close all
format long g
clc
fprintf('LandRate toolbox: Eye movement analysis and landscape rating tool\n\n')
fprintf('LandRate toolbox is a post experimental tool appropriate for the analysis\nof eye movements data and the generation of a relative automatic report.\n')
fprintf('The generated report involves typical eye tracking metrics and\nvisualizations (scanpaths and heatmaps) and an analysis report of\nLandscape Rate Index (LRI) (Krassanakis et al., 2018).\n')
fprintf('The computation of fixation events is based on the the implementation of\nEyeMMV toolbox algortihm (Krassanakis et al., 2014).\n')
fprintf('LandRate toolbox is freely distributed under GPLv.3 license.\n\n')
fprintf('Please prepare your import files according to the examples given\nwith LandRate toolbox before starting import process.\n\n\n')

start_toolbox=input('Press ENTER key to start import process...\n');
tic
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input raw eye tracking data in xlsx file format
%file headers: ID, Passing_time_EyeA(sec), X_EyeA(tracker
%units),Y_EyeA(tracker units), Passing_time_EyeB(sec), X_EyeB(tracker units),Y_EyeB(tracker units), Subject ID, Stimulus ID, Dominant Eye)
%data=xlsread('test_data.xlsx');
%data=xlsread('Both_Eyes_Corrected.xlsx');

data=uigetfile('*.xlsx','Import RAW DATA file in xlsx format');
data=xlsread(data);

%input stimuli names file
stimuli_names=uigetfile('*.txt','Import STIMULI NAMES file in txt format');
stimuli_names=importdata(stimuli_names);
%stimuli_names=importdata('stimuli_names.txt');

%input AOIs (format: [X(pixels)][Y(pixels)][Stimuli ID][AOI ID])
AOI=uigetfile('*.xlsx','Import AOIs file in xlsx format');
AOI=xlsread(AOI);
%AOI=xlsread('AOIs_data.xlsx');

%input parameters
%eye tracking system & fixation detection parameters 
eye_tracker_system_parameters=uigetfile('*.xlsx','Import EYE TRACKING SYSTEM & FIXATION PARAMETERS file in xlsx format');
eye_tracker_system_parameters=xlsread(eye_tracker_system_parameters);
%eye tracker coordinate system parameters
maxx=eye_tracker_system_parameters(1);%horizontal
maxy=eye_tracker_system_parameters(2); %vertical
coord_start=eye_tracker_system_parameters(3); %starting point (0,0) of the coordinate system: down-left:1, up-left:2 
maxx_pixel=eye_tracker_system_parameters(4);
maxy_pixel=eye_tracker_system_parameters(5);
%fixation detection parameters (based on EyeMMV fixation detection algorithm)
t1=eye_tracker_system_parameters(6); %in pixels
t2=eye_tracker_system_parameters(7); %in pixels
minDur=eye_tracker_system_parameters(8); %minimum fixation duration in ms

max_r_fixation_visualization=eye_tracker_system_parameters(9); %maximum radius for fixations visualization in pixels

%parameters for heatmap generation
spacing_coef=eye_tracker_system_parameters(10); 
kernel_size_gaussian=eye_tracker_system_parameters(11);
s_gaussian=eye_tracker_system_parameters(12);

%parameters weights (file)
parameters_weights=uigetfile('*.xlsx','Import PARAMETERS WEIGHTS file in xlsx format');
parameters_weights=xlsread(parameters_weights);
%report file name
report_filename=uiputfile('*.txt','Give the file name to store the FINAL REPORT (txt format)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize report file
if (exist(report_filename))
    delete(report_filename);
end
diary (report_filename)
diary on

%data for all subjects and stimuli
%convert passing times in msec (for both eyes)
data(:,2)=1000*data(:,2);
data(:,5)=1000*data(:,5);

%convert gaze coordinates in pixel system (start point: up-left corner)
cc=maxx_pixel/maxx; %convert coefficient
if coord_start==1
    data(:,3)= cc*data(:,3);
    data(:,4)= maxy_pixel-cc*data(:,4);
    data(:,6)= cc*data(:,6);
    data(:,7)= maxy_pixel-cc*data(:,7);
elseif coord_start==2
    data(:,3)= cc*data(:,3);
    data(:,4)= cc*data(:,4);
    data(:,6)= cc*data(:,6);
    data(:,7)= cc*data(:,7);
end

%total number of records
n_total=size(data);
n_total=n_total(1,1);

id=data(:,1); %record id

%passing time, x, y according to dominant eye values (Eye A:1, Eye B:2,Average: 3)
dom_eye=data(:,10); %dominant eye values
for i=1:n_total
    if dom_eye(i)==1
        time(i)=data(i,2); %passing time (Eye A)
        x(i)=data(i,3); %x coordinate (Eye A)
        y(i)=data(i,4); %y coordinate (Eye A)
    elseif dom_eye(i)==2
        time(i)=data(i,5); %passing time (Eye B)
        x(i)=data(i,6); %x coordinate (Eye B)
        y(i)=data(i,7); %y coordinate (Eye B)
    elseif dom_eye(i)==3
        time(i)=(data(i,2)+data(i,5))/2; %passing time (Eye B)
        x(i)=(data(i,3)+data(i,6))/2; %x coordinate (Eye B)
        y(i)=(data(i,4)+data(i,7))/2; %y coordinate (Eye B)
    else
        fprintf('Please select Eye A (eye:1) or Eye B (eye:2)\n')
    end
end

%convert time,x,y to nx1 from 1xn vectors
time=time';
x=x';
y=y';

subjects=data(:,8); % subjects numbers
stimuli=data(:,9); %stimuli numbers

data=[id time x y subjects stimuli dom_eye]; %rebuild data matrix (after dominant eye selection)
subjects_number=max(subjects); %number of subjects
stimuli_number=max(stimuli); %number of stimuli

%convert stimuli names from the corresponded file to characters
stimuli_names=char(stimuli_names);

%initialize metrics matrix [subjects]x[stimuli]x[metrics number]
metrics_subj_stim=zeros(subjects_number,stimuli_number);

fprintf('LandRate toolbox: Eye movement analysis and landscape rating tool\n\n')
fprintf('Analysis Report\n\n')
fprintf('[1] Subject, Stimuli, and AOIs report:\n')
for i=1:subjects_number
    for j=1:stimuli_number
        
        % data that correspond to each combination of subjects and stimuli
        data_subj_stim=data(subjects==i & stimuli==j,:);
        
        %calculate fixations using EyeMMV' s fixation detection algorithm
        fixation_data_subj_stim=fixation_detection_EyeMMV_modified(data_subj_stim,t1,t2,minDur,maxx,maxy,'t2');
       
        %total number of fixations on visual scene
        n_fixation_data_subj_stim=size(fixation_data_subj_stim);
        n_fixation_data_subj_stim=n_fixation_data_subj_stim(1,1);
        
        fprintf('\nFixation analysis report [Subject: %.0f, Stimulus: %.0f]\n',i,j);
        fprintf('-Total number of fixations on the visual scene: %.0f\n',n_fixation_data_subj_stim)
        fprintf('-Minimum fixation duration on the visual scene: %.2f ms\n',min(fixation_data_subj_stim(:,7)))
        fprintf('-Maximum fixation duration on the visual scene: %.2f ms\n',max(fixation_data_subj_stim(:,7)))
        fprintf('-Average fixation duration on the visual scene: %.2f ms\n',mean(fixation_data_subj_stim(:,7)))
        fprintf('-Total duration of fixations on the visual scene: %.2f ms\n',sum(fixation_data_subj_stim(:,7)))
        fprintf('-Fixation list:\n')
        fprintf('ID-Xcenter(pixels)-Ycenter(pixels)-Nt1-Nt2/3s-StartTime(ms)-EndTime(ms)-Duration(ms)\n')
        for k=1:n_fixation_data_subj_stim
            fprintf('%.0f %.2f %.2f %.0f %.0f %.2f %.2f %.2f\n',k, fixation_data_subj_stim(k,1),fixation_data_subj_stim(k,2),fixation_data_subj_stim(k,3),fixation_data_subj_stim(k,4),fixation_data_subj_stim(k,5),fixation_data_subj_stim(k,6),fixation_data_subj_stim(k,7))
        end
              
        %calculate saccades
        saccades_data_subj_stim=saccade_analysis_EyeMMV_modified(fixation_data_subj_stim);
        n_saccades_data_subj_stim=size(saccades_data_subj_stim);
        n_saccades_data_subj_stim=n_saccades_data_subj_stim(1,1);
        fprintf('\nSaccade analysis report [Subject: %.0f, Stimulus: %.0f]\n',i,j);
        fprintf('-Total number of saccades on the visual scene: %.0f\n',n_saccades_data_subj_stim)
        fprintf('-Minimum saccade duration on the visual scene: %.2f ms\n',min(saccades_data_subj_stim(:,5)))
        fprintf('-Maximum saccade duration on the visual scene: %.2f ms\n',max(saccades_data_subj_stim(:,5)))
        fprintf('-Average saccade duration on the visual scene: %.2f ms\n',mean(saccades_data_subj_stim(:,5)))
        fprintf('-Total duration of saccades on the visual scene: %.2f ms\n',sum(saccades_data_subj_stim(:,5)))     
        fprintf('-Minimum saccade amplitude on the visual scene: %.2f pixels\n',min(saccades_data_subj_stim(:,6)))
        fprintf('-Maximum saccade amplitude on the visual scene: %.2f pixels\n',max(saccades_data_subj_stim(:,6)))
        fprintf('-Average saccade amplitude on the visual scene: %.2f pixels\n',mean(saccades_data_subj_stim(:,6)))        
        fprintf('-Minimum saccade direction angle on the visual scene: %.2f degrees\n',min(saccades_data_subj_stim(:,7)))
        fprintf('-Maximum saccade direction angle on the visual scene: %.2f degrees\n',max(saccades_data_subj_stim(:,7)))
        fprintf('-Average saccade direction angle on the visual scene: %.2f degrees\n',mean(saccades_data_subj_stim(:,7)))        
        fprintf('-Saccade list:\n')
        fprintf('ID-X_Start_Point(pixels)-Y_Start_Point(pixels)-X_End_Point(pixels)-Y_End_Point(pixels)\n-Duration(ms)-Amplitude(pixels)-Direction_angle(degrees)-Start_Fixation-End_Fixation)\n')        
        for s=1:n_saccades_data_subj_stim
            fprintf('%.f %.4f %.4f %.4f %.4f %.1f %.3f %.3f %.f %.f\n',s,saccades_data_subj_stim(s,1),saccades_data_subj_stim(s,2),saccades_data_subj_stim(s,3),saccades_data_subj_stim(s,4),saccades_data_subj_stim(s,5),saccades_data_subj_stim(s,6),saccades_data_subj_stim(s,7),saccades_data_subj_stim(s,8),saccades_data_subj_stim(s,9))
        end
        
        %scapath statistics
        fprintf('\nScanpath analysis report [Subject: %.0f, Stimulus: %.0f]\n',i,j)
        fprintf('-Total scanpath length: %.2f pixels\n',sum(saccades_data_subj_stim(:,6)))
        fprintf('-Total scanpath duration: %.2f ms\n',sum(fixation_data_subj_stim(:,7))+sum(saccades_data_subj_stim(:,5)))
        fprintf('-Ratio saccades/fixations (durations): %.3f\n',(sum(saccades_data_subj_stim(:,5)))/(sum(fixation_data_subj_stim(:,7))))
           
        %initialize scanpath plot
        figure ('Name','Scanpath visualization')
        imshow(stimuli_names(j,:),'InitialMagnification','fit')
  
        %generate metrics matrix (all subjects and stimuli combinations)
        metrics_subj_stim(i,j,1)=n_fixation_data_subj_stim; %Metric 1: Total number of fixations on visual scene
        metrics_subj_stim(i,j,2)=min(fixation_data_subj_stim(:,7)); %Metric 2: Minimum fixation duration on visual scene
        metrics_subj_stim(i,j,3)=max(fixation_data_subj_stim(:,7)); %Metric 3: Maximum fixation duration on visual scene
        metrics_subj_stim(i,j,4)=mean(fixation_data_subj_stim(:,7)); %Metric 4: Average fixation duration on visual scene
        metrics_subj_stim(i,j,5)=sum(fixation_data_subj_stim(:,7)); %Metric 5: Total duration of fixations on the visual scene
        metrics_subj_stim(i,j,6)=n_saccades_data_subj_stim;%Metric 6: Total number of saccades on the visual scene
        metrics_subj_stim(i,j,7)=min(saccades_data_subj_stim(:,5));%Metric 7: Minimum saccade duration on the visual scene 
        metrics_subj_stim(i,j,8)=max(saccades_data_subj_stim(:,5));%Metric 8: Maximum saccade duration on the visual scene 
        metrics_subj_stim(i,j,9)=mean(saccades_data_subj_stim(:,5));%Metric 9: Average saccade duration on the visual scene
        metrics_subj_stim(i,j,10)=sum(saccades_data_subj_stim(:,5));%Metric 10: Total duration of saccades on the visual scene
        metrics_subj_stim(i,j,11)=min(saccades_data_subj_stim(:,6));%Metric 11: Minimum saccade amplitude on the visual scene
        metrics_subj_stim(i,j,12)=max(saccades_data_subj_stim(:,6));%Metric 12: Maximum saccade amplitude on the visual scene
        metrics_subj_stim(i,j,13)=mean(saccades_data_subj_stim(:,6));%Metric 13: Average saccade amplitude on the visual scene
        metrics_subj_stim(i,j,14)=min(saccades_data_subj_stim(:,7));%Metric 14: Minimum saccade direction anlge on the visual scene
        metrics_subj_stim(i,j,15)=max(saccades_data_subj_stim(:,7));%Metric 15: Maximum saccade direction anlge on the visual scene
        metrics_subj_stim(i,j,16)=mean(saccades_data_subj_stim(:,7));%Metric 16: Average saccade direction anlge on the visual scene
        metrics_subj_stim(i,j,17)=sum(saccades_data_subj_stim(:,6));%Metric 17: Total scanpath length
        metrics_subj_stim(i,j,18)=sum(fixation_data_subj_stim(:,7))+sum(saccades_data_subj_stim(:,5));%Metric 18: Total scanpath duration
        metrics_subj_stim(i,j,19)=(sum(saccades_data_subj_stim(:,5)))/(sum(fixation_data_subj_stim(:,7)));%Metric 19: Ratio saccades/fixations (durations)
                      
        %AOI analysis 
        AOI_visual_scene=AOI(AOI(:,3)==j,:); %all AOIs on the visual scene
        number_AOI_visual_scene=max(AOI_visual_scene(:,4));
        for m=1:number_AOI_visual_scene
            AOI_visual_scene_id=AOI_visual_scene(AOI_visual_scene(:,4)==m,:); %Seperate AOIs on the visual scene
            hold on
            fill(AOI_visual_scene_id(:,1),AOI_visual_scene_id(:,2),'w','FaceColor','none','EdgeColor','w','LineStyle','--')
            hold on
            text(mean(AOI_visual_scene_id(:,1)),mean(AOI_visual_scene_id(:,2)),strcat('AOI:',num2str(m)),'Color','w','HorizontalAlignment','center')
            
            %Fixations inside each AOI
            fixations_in_AOI=inpolygon(fixation_data_subj_stim(:,1),fixation_data_subj_stim(:,2),AOI_visual_scene_id(:,1),AOI_visual_scene_id(:,2)); %output:logical values (0: out of AOI,1:inside or on the edge of AOI)
            %initialize fixation list in AOI
            fixation_list_in_AOI=[];
            for w=1:n_fixation_data_subj_stim
                if fixations_in_AOI(w)==1
                    fixation_list_in_AOI=[fixation_list_in_AOI;fixation_data_subj_stim(w,:)];
                end
            end
            
            fprintf('\nAOI analysis [Subject: %.0f, Stimulus: %.0f, AOI: %.0f]\n',i,j,m);                    
            %number of fixations in AOI
            n_fixation_list_in_AOI=size(fixation_list_in_AOI);
            n_fixation_list_in_AOI=n_fixation_list_in_AOI(1,1);
            if n_fixation_list_in_AOI==0
                fprintf('-There are no fixations inside or on the edge of the AOI\n')
            else
                fprintf('-Total number of fixations inside or on the edge of AOI: %.0f\n',n_fixation_list_in_AOI)          
                fprintf('-Minimum fixation duration inside or on the edge of AOI: %.2f ms\n',min(fixation_list_in_AOI(:,7)))
                fprintf('-Maximum fixation duration inside or on the edge of AOI: %.2f ms\n',max(fixation_list_in_AOI(:,7)))
                fprintf('-Average fixation duration inside or on the edge of AOIe: %.2f ms\n',mean(fixation_list_in_AOI(:,7)))
                fprintf('-Total duration of fixations inside or on the edge of AOI: %.2f ms\n',sum(fixation_list_in_AOI(:,7)))
                fprintf('-Fixation list inside or on the edge of AOI:\n')
                fprintf('ID-Xcenter(pixels)-Ycenter(pixels)-Nt1-Nt2/3s-StartTime(ms)-EndTime(ms)-Duration(ms)\n')
                for p=1:n_fixation_list_in_AOI
                    fprintf('%.0f %.2f %.2f %.0f %.0f %.2f %.2f %.2f\n',p, fixation_list_in_AOI(p,1),fixation_list_in_AOI(p,2),fixation_list_in_AOI(p,3),fixation_list_in_AOI(p,4),fixation_list_in_AOI(p,5),fixation_list_in_AOI(p,6),fixation_list_in_AOI(p,7))
                end
                                          
            end                

        end
        
        %finalize scanpath visualization
        scan_path_visualization_EyeMMV_modified(fixation_data_subj_stim,max_r_fixation_visualization)
        title(['Scanpath visualization [Subject:',num2str(i),', Stimulus:',num2str(j),']'])
        xlabel('* Fixations durations values are presented in ms', 'FontSize',10)
            
    end
end

%compute visual scene parameters as the average value of each metric based on all subjects 
p1=mean(metrics_subj_stim(:,:,1));
p2=mean(metrics_subj_stim(:,:,2));
p3=mean(metrics_subj_stim(:,:,3));
p4=mean(metrics_subj_stim(:,:,4));
p5=mean(metrics_subj_stim(:,:,5));
p6=mean(metrics_subj_stim(:,:,6));
p7=mean(metrics_subj_stim(:,:,7));
p8=mean(metrics_subj_stim(:,:,8));
p9=mean(metrics_subj_stim(:,:,9));
p10=mean(metrics_subj_stim(:,:,10));
p11=mean(metrics_subj_stim(:,:,11));
p12=mean(metrics_subj_stim(:,:,12));
p13=mean(metrics_subj_stim(:,:,13));
p14=mean(metrics_subj_stim(:,:,14));
p15=mean(metrics_subj_stim(:,:,15));
p16=mean(metrics_subj_stim(:,:,16));
p17=mean(metrics_subj_stim(:,:,17));
p18=mean(metrics_subj_stim(:,:,18));
p19=mean(metrics_subj_stim(:,:,19));
parameters_all_stimuli=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;p16;p17;p18;p19]';
parameters_all_stimuli_normalized=values_normalization(parameters_all_stimuli);

%compute Landscape Rate Index (LRI) for each stimuli
wp_parameters=(parameters_weights/sum(abs(parameters_weights)).*parameters_all_stimuli_normalized);
LRI=(sum(wp_parameters'))'; %all stimuli

fprintf('\n\n[2] Landscape Rate Index (LRI) report:\n\n') 
fprintf('LRI model\n')

fprintf('-Fixations based weights/parameters:')
fprintf('\nw1/p1: Total number of fixations on the visual scene')
fprintf('\nw2/p2: Minimum fixation duration on the visual scene')
fprintf('\nw3/p3: Maximum fixation duration on the visual scene')
fprintf('\nw4/p4: Average fixation duration on the visual scene')
fprintf('\nw5/p5: Total duration of fixations on the visual scene')

fprintf('\n-Saccades based weights/parameters:')
fprintf('\nw6/p6: Total number of saccades on the visual scene')
fprintf('\nw7/p7: Minimum saccade duration on the visual scene')
fprintf('\nw8/p8: Maximum saccade duration on the visual scene')
fprintf('\nw9/p9: Average saccade duration on the visual scene')
fprintf('\nw10/p10: Total duration of saccades on the visual scene')
fprintf('\nw11/p11: Minimum saccade amplitude on the visual scene')
fprintf('\nw12/p12: Maximum saccade amplitude on the visual scene')
fprintf('\nw13/p13: Average saccade amplitude on the visual scene')
fprintf('\nw14/p14: Minimum saccade direction angle on the visual scene')
fprintf('\nw15/p15: Maximum saccade direction angle on the visual scene')
fprintf('\nw16/p16: Average saccade direction angle on the visual scene')

fprintf('\n-Scanpath based wights/parameters:')
fprintf('\nw17/p17: Total scanpath length')
fprintf('\nw18/p18: Total scanpath duration')
fprintf('\nw19/p19: Ratio saccades/fixations (durations)\n')

%print all average parameters of each stimulus
fprintf('\nAverage values of parameters (before normalization) for each stimuli\n')
fprintf('p1-p2-p3-p4-p5-p6-p7-p8-p9-p10-p11-p12-p13-p14-p15-p16-p17-p18-p19\n')
for i_all_par=1:stimuli_number
fprintf('Stim.%0.f:%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n',i_all_par, parameters_all_stimuli(i_all_par,1),parameters_all_stimuli(i_all_par,2),parameters_all_stimuli(i_all_par,3),parameters_all_stimuli(i_all_par,4),parameters_all_stimuli(i_all_par,5),parameters_all_stimuli(i_all_par,6),parameters_all_stimuli(i_all_par,7),parameters_all_stimuli(i_all_par,8),parameters_all_stimuli(i_all_par,9),parameters_all_stimuli(i_all_par,10),parameters_all_stimuli(i_all_par,11),parameters_all_stimuli(i_all_par,12),parameters_all_stimuli(i_all_par,13),parameters_all_stimuli(i_all_par,14),parameters_all_stimuli(i_all_par,15),parameters_all_stimuli(i_all_par,16),parameters_all_stimuli(i_all_par,17),parameters_all_stimuli(i_all_par,18),parameters_all_stimuli(i_all_par,19))
end

%normalized weights
weights_normalized=parameters_weights/sum(abs(parameters_weights));
n_weights=length(weights_normalized);
fprintf('\nNormalized values of selected weights (between 0 and 1) for LRI computation:\n')
for i_par_weights=1:n_weights
fprintf('w%.0f: %.3f\n',i_par_weights,weights_normalized(i_par_weights))
end

fprintf('\nComputation of LRI for each stimuli\n')
for i_lri=1:stimuli_number
 fprintf('LRI (Stimulus %.0f): %.3f \n',i_lri,LRI(i_lri))
end

%Generate heatmap visualizations for each stimulus
for g=1:stimuli_number
    data_heatmap=data(stimuli==g,:);
    %Heatmap visualization  
    [heatmapRGB,heatmap_visualization]=heatmap_generator_EyeMMV_modified(data_heatmap(:,3:4),stimuli_names(g,:),spacing_coef,maxx_pixel,maxy_pixel,kernel_size_gaussian,s_gaussian,g);    
end

%Plot LRI disstribution in bar chart
figure('Name','LRI values distribution')
bar(categorical(cellstr(stimuli_names)), LRI)
title('Landscape Rate Index (LRI) values distribution')
xlabel('Stimulus No', 'FontSize',10)
ylabel('LRI value', 'FontSize',10)
alpha(0.8)
grid on

fprintf('\nReport created successfully\n')
date_time_of_report=datestr(now);
fprintf(date_time_of_report)
fprintf('\n')
toc
fprintf('\n\n----------------------------------------------------------------------\n')
fprintf('LandRate toolbox: Eye movement analysis and landscape rating tool\n')
fprintf('Copyright (C) 2017 Vassilios Krassanakis\n\n') 
fprintf('This program is free software: you can redistribute it and/or modify\n')
fprintf('it under the terms of the GNU General Public License as published by\n')
fprintf('the Free Software Foundation, either version 3 of the License, or\n')
fprintf('(at your option) any later version.\n\n')
fprintf('This program is distributed in the hope that it will be useful,\n')
fprintf('but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
fprintf('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
fprintf('GNU General Public License for more details.\n\n') 
fprintf('You should have received a copy of the GNU General Public License\n')
fprintf('along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\n') 
fprintf('For further information, please email me at [krasanakisv@gmail.com], [krasvas@mail.ntua.gr], or [krasvas@teiath.gr] \n\n')
%close report file
diary off
%close all
