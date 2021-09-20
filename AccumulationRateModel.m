
%% Set-up

filename_in = sprintf('Thar_data.txt'); % data file
repeats = 1e4; % 1e6 suggested for final versions
resolution = 300; % 500 suggested for final versions
scaling = 1e8; % This is arbitrary. 1e8 works well for plotting and reporting numbers more easily

% Plot constraints
X_max = 60;
Y_max = 60;


%% -------import & re-format data--------------------------------------------

data = load(filename_in);
[rows, cols] = size(data);
data(:,1) = data(:,1)+0.001;
reform_data = zeros(rows,rows*3);
current_strat = 1;
current_row=1;

for i = 1:rows
    
    depth = data(i,1);
    age = data(i,2);
    error = data(i,3);
    
    if age == 0
        age = age + 0.0001;
    end
    
    if i == 1 
        reform_data(i,1) = depth; 
        reform_data(i,2) = age; 
        reform_data(i,3) = error;
        previous_depth = depth;
    end 
    if i > 1
        if depth > previous_depth
            current_row = numel(find(reform_data(:,current_strat))) + 1;
            reform_data(current_row, [current_strat, current_strat + 1, current_strat+2]) = [depth age error];
        else
            current_strat=current_strat + 3;
            current_row = 1;
            reform_data(current_row, [current_strat,current_strat + 1, current_strat + 2]) = [depth age error];
        end
        previous_depth = depth;
    end
end

%clean table
nc = sum(reform_data(1,:) > 0);
reform_data(:, nc+1 : rows * 3) = [];
occupied = reform_data > 0;
a = sum(occupied,1);
b = max(a);
reform_data(b + 1:rows,:) = [];

% optional save
%save 'test.txt' reform_data -ASCII -tabs -double

%% -------process data-------------------------------------------------------
data=reform_data; %columns - depth, age, error
[r,c] = size(data);

cores = c / 3;

upper_age = max(max(data)) * 1.5;
bin_centres = resolution / 2 : resolution:upper_age;
hist_bins = numel(bin_centres);
summed = zeros(hist_bins,1);

EachCore=zeros(hist_bins,cores);

for j=1:cores
    
    fprintf('Progress: %.0f %%  ', 100*j / cores);
    
    sampled = data(:,2 + ((j - 1) * 3));
    errors = data(:,3 + ((j - 1) * 3));
    rel_errors = errors ./ sampled;
    depths = data(:,1 + ((j - 1) * 3));
    ns = numel(sampled);
    rates_table = zeros(r,repeats);
    dates_table = zeros(r,repeats);
    ddepth_table = zeros(r,repeats);
    
    rand_array = random('Normal',0,1,ns * repeats,1);
    
    for i = 1:repeats
        rand_value = rand_array(((i-1) * ns) + 1 : ((i - 1) * ns) + ns);
        sampledr = sampled + (rand_value .* errors);
        diff_sampled = diff(sampledr);
        diff_depths = abs(diff(depths));
        
        rates = diff_depths ./ diff_sampled;       
        n_diff = numel(diff_sampled);
        rates_table(1:n_diff,i) = rates;
        dates_table(1:n_diff,i) = sampledr(2 : n_diff + 1);      
        ddepth_table(1:n_diff,i) = diff_depths;
    end
    
    n=numel(rates_table);
    rates_dat = reshape(rates_table,n,1);
    dates_dat = reshape(dates_table,n,1);
    ddepth_dat = reshape(ddepth_table,n,1);
    
    del = rates_dat == 0;
    rates_dat(del) = [];
    dates_dat(del) = [];
    ddepth_dat(del) = [];
    
    % Optional:
    %del=rates_dat>0.02;
    %rates_dat(del)=0.02; %truncate high values of accumulation rate at 0.02
    
    test_data = dates_dat; 
    rates_data = rates_dat;
    
    for k = 1:hist_bins
        a = test_data > ((k-1) * resolution) & test_data <= (k * resolution);
        std_dat = rates_data(a);
        n_a = sum(a);        
        total_acc = sum( ddepth_dat(a)) * n_a / repeats;
                
        if numel(std_dat) > 0
            EachCore(k,j) = ((iqr(std_dat) * total_acc) / repeats) / (resolution^2);
        end
        
    end
    
    fprintf(' (core %d complete)\n', j); 
    
    
end

m = (max(max(data)) / hist_bins) : (upper_age / hist_bins) : upper_age;
m = m';

stdevs_acc = sum(EachCore,2) * scaling;
smax = max(stdevs_acc);
ls = log(stdevs_acc); 

del = ls == -inf();
ls(del) = 1e10;
ls(del) = min(ls);
record = zeros(numel(m),4);
record(:,1:4) = [m stdevs_acc m ls];


%% PLOTS

%------ (auto X,Y) ------------
figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',20);
box(axes1,'on');
hold(axes1,'on');
ylabel('Accumulation intensity (au)');
xlabel('Age (ka)');
area(record(:,1)/1000,record(:,2),0);

%------ (constrained X, log Y)------------
figure2 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure2,'FontSize',20);

xlim(axes1,[0 X_max]);

box(axes1,'on');
hold(axes1,'on');
ylabel('Accumulation intensity (au)');
xlabel('Age (ka)');
plot(record(:,3)/1000,record(:,4));

%------ (constrained X,Y)------------
figure3 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure3,'FontSize',20);

xlim(axes1,[0 X_max]);
ylim(axes1,[0 Y_max]);

box(axes1,'on');
hold(axes1,'on');
ylabel('Accumulation intensity (au)');
xlabel('Age (ka)');
area(record(:,1)/1000,record(:,2),0);

%------ (constrained X)------------
figure4 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure4,'FontSize',20);

xlim(axes1,[0 X_max]);

box(axes1,'on');
hold(axes1,'on');
ylabel('Accumulation intensity (au)');
xlabel('Age (ka)');
area(record(:,1)/1000,record(:,2),0);


% disp(['Max = ', num2str(max(record(:,2))*1e3)]);

