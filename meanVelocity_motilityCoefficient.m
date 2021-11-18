clear all;
close all;
warning off;

%% parameters and input fields
% filename of lymph node tracks
data_filename = '20210828-M1L-ZT2-z220-259-green';

% filename of frame position
reference_filename = '20210828-M1L-ZT2-z220-259-ref';

% number of frames in video
frames = 232;

% time between frames
t_scale=15.641; % s/frame

% number of tracks for analysis
analy_tracks=30;

%% Read data from files
x = xlsread(data_filename,'A:A');
y = xlsread(data_filename,'B:B');
z = xlsread(data_filename,'C:C');
frame = xlsread(data_filename,'G:G');
ID = xlsread(data_filename,'H:H');
tracks = unique(ID);

x_ref = xlsread(reference_filename,'A:A');
y_ref = xlsread(reference_filename,'B:B');
z_ref = xlsread(reference_filename,'C:C');
frame_ref = xlsread(reference_filename,'G:G');

%% sort data by track and time order, subtract frame position offset
marker = -10000;
x_adj = marker*ones(length(tracks),frames);
y_adj = marker*ones(length(tracks),frames);
z_adj = marker*ones(length(tracks),frames);

x_ref_adj=zeros(length(x_ref),1);
y_ref_adj=zeros(length(x_ref),1);
z_ref_adj=zeros(length(x_ref),1);
    
% sort reference positions by time order
for index = 1:length(x_ref)
    x_ref_adj(frame_ref(index))=x_ref(index);
    y_ref_adj(frame_ref(index))=y_ref(index);
    z_ref_adj(frame_ref(index))=z_ref(index);
end

% Go through all data points
for index = 1:length(ID)
        
    % Check to see if the current data point is on a track we want to parse
    for track_index = 1:length(tracks)
        if (ID(index) == tracks(track_index))
            % sort and offset x,y,z coordinates by frame time and track
            x_adj(track_index,frame(index)) = x(index) - x_ref_adj(frame(index));
            y_adj(track_index,frame(index)) = y(index) - y_ref_adj(frame(index));
            z_adj(track_index,frame(index)) = z(index) - z_ref_adj(frame(index));
            
            break;
        end
    end
end

%% Calculate different displacements within each track and sort by travel time
displacements = marker*ones(analy_tracks,(frames-1));
errors = zeros(analy_tracks,(frames-1));
displacement_count = 0;
displacement = zeros(1,(frames-1));  

% Cycle through tracks
num_tracks=analy_tracks;
if analy_tracks> length(tracks), num_tracks=length(tracks);end
track_time = zeros(1,num_tracks);
v_mean = zeros(1,length(tracks));

for index = 1:num_tracks
    start_index = 1;
    while (x_adj(index,start_index) == marker)
        start_index = start_index + 1;
    end
    
    end_index = frames;
    while (x_adj(index,end_index) == marker)
        end_index = end_index - 1;
    end
    
    x_temp = x_adj(index,:);
    y_temp = y_adj(index,:);
    z_temp = z_adj(index,:);
    
    max_steps = end_index - start_index;
    for num_steps = 1:max_steps
        displacement_count = 0;
        for index_2 = start_index:(end_index-num_steps)
            if ((x_temp(index_2) ~= marker) && (x_temp(index_2+num_steps) ~= marker))
                displacement_count = displacement_count + 1;
                displacement(displacement_count) = sqrt((x_temp(index_2) - x_temp(index_2+num_steps))^2+(y_temp(index_2) - y_temp(index_2+num_steps))^2+(z_temp(index_2) - z_temp(index_2+num_steps))^2);
            end
        end
        
        displacements(index,num_steps) = mean(displacement(1:displacement_count));
        errors(index,num_steps) = std(displacement(1:displacement_count));
    end
    
    % Calculate tracking time for all tracks
    valid_points1=x_temp~=marker;
    track_time(index)=length(x_temp(valid_points1))-1;
    
    % Calculate the separation between points and mean velocity
    x_temp2 = x_adj(index,valid_points1);
    y_temp2 = y_adj(index,valid_points1);
    z_temp2 = z_adj(index,valid_points1);
    d = sqrt((x_temp2(2:end)-x_temp2(1:(end-1))).^2 + (y_temp2(2:end)-y_temp2(1:(end-1))).^2 + (z_temp2(2:end)-z_temp2(1:(end-1))).^2);
    v_mean(index) = 60/t_scale*sum(d)/(end_index-start_index);
end

%% Plot displacments for each tracks
% xt=1:19;
% xsqrt=sqrt(xt/2);
% for index=1:30
% temp2=displacements(index,:);
% valid_points=temp2~=marker;
% plot(xsqrt(1:length(temp2(valid_points))),temp2(valid_points));
% c=polyfit(xsqrt(1:length(temp2(valid_points))),temp2(valid_points),1);
% slope(index)=c(1);
% hold on;
% end
% hold off;

figure();
xt=1:(frames-1);
xsqrt=sqrt(xt*t_scale/60); % unit: min^(1/2)
mean_track_time=round(mean(track_time(1:num_tracks)));
slope=zeros(1,analy_tracks);
start=2; % start time for linear fitting

for index=1:num_tracks
temp2=displacements(index,:);
valid_points=temp2~=marker;
plot(xsqrt(1:length(temp2(valid_points))),temp2(valid_points));

% calculate slope by linear fitting from start to mean_track_time
valid_points(1:start-1)=0;
valid_points(mean_track_time+1:frames-1)=0;
fit_end=track_time(index);
if fit_end > mean_track_time, fit_end=mean_track_time;end
c=polyfit(xsqrt(start:fit_end),temp2(valid_points),1);
slope(index)=c(1);
hold on;
end
hold off;

%% save files
filename_result= data_filename(1:end-4)+"-result.xls";
xlswrite(filename_result,displacements,1);
xlswrite(filename_result,transpose(track_time),2);
xlswrite(filename_result,transpose(v_mean),3);

savefig(data_filename(1:end-4));


motility_coefficient=(slope.^2)/6; %motility coefficient = x^2/6t [um^2/min]
xlswrite(filename_result,transpose(motility_coefficient),4);
% xlswrite('Motility_Coefficient',transpose(slope),2);

% %% Calculate mean and std of displacement for each travel time
% 
% mean_track_time=round(mean(track_time(1:num_tracks))); 
% mean_displacement = zeros((mean_track_time),1);
% errors_2 = zeros((mean_track_time),1);
% for index = 1:mean_track_time
%     temp = displacements(:,index);
%     valid_points = temp~=marker;    
% 
%     mean_displacement(index) = mean(temp(valid_points));
%     errors_2(index) = std(temp(valid_points));
% 
% end
% 
% time = t_scale*(1:length(mean_displacement))/60;
% % time = t_scale*(1:(frames-1))/60;
% %% Linear fit
% X = [ones(length(mean_displacement),1), sqrt(time)'];
% Y = mean_displacement;
% C = inv(X'*X)*X'*Y;
% 
% %% Plot Results
% figure();
% hold on;
% scatter(sqrt(time),mean_displacement);
% handle = errorbar(sqrt(time),mean_displacement,errors_2);
% handle.LineStyle = 'none';
% plot([sqrt(time(1)),sqrt(time(end))],[(C(1)+C(2)*sqrt(time(1))),(C(1)+C(2)*sqrt(time(end)))],'k--','LineWidth',1.5);
% xlabel('squared root of time (min^{1/2})');
% ylabel(strcat('mean displacement (','\mu','m)'));
% % xlim([0,ceil(sqrt(time(end)))]);
% xlim([0,ceil(sqrt(t_scale*(frames-1)/60))]);
% 
% 
% % sets tick marks on top and right
% ax1=gca;
% ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
% set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
% set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
% 
% % removes tick labels for top and right
% OppTickLabels = {' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '};
% set(ax2, 'XTickLabel', OppTickLabels,'YTickLabel',OppTickLabels);
% hold off;
% 
% 
% disp('Slope of Fit (um/sqrt(s)): ');
% disp(C(2));
% 
% %% Save the results 
% savefig(data_filename(1:end-4));
% filename_result= data_filename(1:end-4)+"-result.xls";
% % filename_result= "Mean_displacement_rawdata";
% xlswrite(filename_result,mean_displacement,1,'A');
% xlswrite(filename_result,errors_2,1,'B');
% xlswrite(filename_result,transpose(track_time),1,'C');
% xlswrite(filename_result,C(2),1,'D');
% xlswrite(filename_result,transpose(displacements),1,'E');
