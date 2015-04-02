%Extract postseismic signals from El Mayor-Cucapah GPS timeseries

clear all
close all
clc

%What to do with the outputs
options = struct(...
    'Plot_timeseries',true,...
    'Save_plots',true,...
    'Export_cumulative_postseismic_offsets',true,...
    'Export_coseismic_offsets',true,...
    'Export_decay_timeseries',true,...
    'Map_modeled_offsets',true,...
    'Save_maps',false...
    );

%Time of El Mayor-Cucapah earthquake in years
elmayor_time = 2010.2589;

%Names of components for plotting
components = {'north','east','up'};

%Get station names and select stations within specified radius
flocations=fopen('../gps/gps_within200km.dat');
location_data=textscan(flocations,'%s %f %f %f','commentstyle','#');
fclose(flocations);
stations = struct('staNam',[],'east',[],'north',[]);
stations.names = location_data{1};
stations.east = location_data{2};
stations.north = location_data{3};

N=length(stations.names);

%Set up output files of postseismic decay coefficients
coseismic_output = NaN(N,6);
%Set up output files of coseismic decay coefficients
postseismic_output = NaN(N,6);

%.mat file of evident offsets in timeseries that do not correspond to listed offsets in SOPAC timeseries or maintenance dates in station logs
load otheroffsets.mat

%We run multiple lsqcurvefit iterations with different initial values for the decay parameters. This is the # of iterations.
searchnum = 10;

%To assess the uncertainty in the fit, we add a vector of Gaussian noise * the timeseries uncertainty to the data and then run
%multiple lsqcurvefit iterations with different initial values on that "noise data." This is how many different "noise data" vectors we create.
noisenum = 4;

%Stations being run (to run all, choose 1 and N for start and finish)
start = 1;
finish = N;

for staNum=start:finish
    
    %Name of timeseries file
    fdata = ['../gps/WNAM_Filter_TrendNeuTimeSeries_sopac_20150212/' stations.names{staNum} 'FilterTrend.neu'];
    
    %If this station does not have a timeseries file
    if 0==exist(fdata,'file')
        display([stations.names{staNum} ': No data in SOPAC timeseries'])
        continue
    end
    
    fid = fopen(fdata);
    input = textscan(fid,'%f %*f %*f %f %f %f %f %f %f','commentstyle','#');
    fclose(fid);
    
    %If this station's timeseries stops before the El Mayor earthquake
    if input{1}(end) < elmayor_time
        display([stations.names{staNum} ': No postseismic data in sopac timeseries'])
        continue
    end
    
    %If this station's timeseries starts <2 years before the El Mayor earthquake
    if input{1}(1) > elmayor_time - 2
        display([stations.names{staNum} ': <2 years of interseismic data in sopac timeseries'])
        continue
    end
    
    data = struct('data',[],'sigma',[],'noise',[]);
    
    %Import time vector
    time = input{1} - elmayor_time;
    
    %Import position timeseries and uncertainties
    for component = 1:3
        %Import position values in m
        data.data(:,component) = 1e-3*input{component+1}(1:length(time));
        %Approximate implicit data uncertainty as 2x SOPAC formal error
        data.sigma(:,component) = 2e-3*input{component+4}(1:length(time));
    end
    
    %Initialize output postseismic timeseries
    decay_timeseries = NaN(length(time(time>=0)),7);
    decay_timeseries(:,1) = time(time>=0);
    
    clear input
    
    %% SOPAC automatically fits for amplitudes and dates of coseismic offsets in its timeseries and subtracts those out before publishing the timeseries.
    %This fitting only uses one postseismic decay term, either an exp or a log (we use both) and as a consequence visibly "cuts the corner"
    %at the beginning of the postseismic period in many timeseries, as visible in the SOPAC GPS Explorer app.
    %This means that the postseismic terms are also wrong. These pre-subtracted offsets are reported in the headers of the timeseries files;
    %this part of the code adds them back into the timeseries so that they can be fit with everything else.
    
    %Read in timeseries file again
    fid = fopen(fdata);
    offsetfile = textscan(fid,'%11s');
    fclose(fid);
    offsetfile = offsetfile{:};
    
    offsets = struct('indices',[],'amplitudes',[],'dates',[],'alldates',[],'sopac_coseismic',[]);
    
    %Find locations of offsets in file
    offsets.indices = find(1==strcmp(offsetfile,'offset'));
    
    %When offset amplitudes have <3 digits before the decimal, they are separated from the character string before them with a space in the SOPAC timeseries file;
    %when they have 3 or more digits, they are not. Make this separation.
    for m = 1:length(offsets.indices)-1
        offsets.indices = find(1==strcmp(offsetfile,'offset'));
        if length(char(offsetfile(offsets.indices(m)+1))) > 3
            splitentry = strsplit(char(offsetfile(offsets.indices(m)+1)),':');
            offsetfile = [offsetfile(1:offsets.indices(m));splitentry';offsetfile(offsets.indices(m)+2:end)];
        end
        if length(char(offsetfile(offsets.indices(m)+3))) > 3
            splitentry = strsplit(char(offsetfile(offsets.indices(m)+3)),'-');
            offsetfile = [offsetfile(1:offsets.indices(m)+2);splitentry';offsetfile(offsets.indices(m)+4:end)];
        end
    end
    
    %Find locations of offsets in modified file
    offsets.indices = find(1==strcmp(offsetfile,'offset'));
    offsets.indices(end) = [];
    %The offsets are reported in separate sections for north, east and up.
    %Find the section breaks (to separate by components)
    component_entries = find(1==strcmp(offsetfile,'component'));
    
    %Initialize amplitudes and dates of offsets separated by component
    offsets.amplitudes = cell(0,0);
    offsets.dates = cell(0,0);
    %offsets.otherdates = cell(3,1);
    offsets.savedates = NaN(0,1);
    
    %Offset amplitudes and dates in north component
    offsets.amplitudes{1} = offsetfile(offsets.indices(offsets.indices<component_entries(2))+2);
    offsets.dates{1} = strrep(strrep(offsetfile(offsets.indices(offsets.indices<component_entries(2))+7),'[',''),']','');
    
    %Offset amplitudes and dates in east component
    offsets.amplitudes{2} = offsetfile(offsets.indices(and(offsets.indices>component_entries(2),offsets.indices<component_entries(3)))+2);
    offsets.dates{2} = strrep(strrep(offsetfile(offsets.indices(and(offsets.indices>component_entries(2),offsets.indices<component_entries(3)))+7),'[',''),']','');
    
    %Offset amplitudes and dates in up component
    offsets.amplitudes{3} = offsetfile(offsets.indices(offsets.indices>component_entries(3))+2);
    offsets.dates{3} = strrep(strrep(offsetfile(offsets.indices(offsets.indices>component_entries(3))+7),'[',''),']','');
    
    %Add offsets into timeseries and save dates to fit later except if date is El Mayor mainshock date (which we always fit as a step function anyway)
    for component = 1:3
        for k = 1:length(offsets.amplitudes{component})
            data.data(:,component) = data.data(:,component) + 1e-3*str2num(offsets.amplitudes{component}{k})*heavi(time - (str2num(offsets.dates{component}{k}) - elmayor_time));
            if 0==(str2num(offsets.dates{component}{k})==elmayor_time)
                offsets.savedates(end+1) = str2num(offsets.dates{component}{k}) - elmayor_time;
            end
            
        end
    end
    
    offsets.savedates = unique(offsets.savedates)';
    
    clear offsetfile
    clear offsets.indices
    clear offsets.amplitudes
    clear offsets.dates
    
    %% Other offsets are caused by station maintenance, e.g. in April 2011. The maintenance dates are reported in the station logs.
    %This part of the code imports those station maintenance dates to fit them as offsets later.
    logdates = zeros(0,1);
    flog = (['../gps/stationlogs_2015/' stations.names{staNum} '.log']);
    if 2==exist(flog,'file')
        flog = fopen(flog);
        logfile = textscan(flog,'%10s');
        fclose(flog);
        logfile = logfile{:};
        
        %Find descriptions of maintenance in log file
        removed = find(1==strcmp(logfile,'Removed'));
        removed = removed(1==strcmp(logfile(removed-1),'Date'));
        removed = removed(0==strcmp(logfile(removed+2),'(CCYY-MM-D'));
        removed = removed(0==strcmp(logfile(removed+2),'CCYY-MM-DD'));
        
        %Grab dates of maintenance
        logdates = datenum(strrep(logfile(removed+2), '-', ','))/365.243 - elmayor_time;
        clear logfile
    end
    
    % Add station log dates and dates of unexplained offsets for each component
    
    offsets.savedates = [offsets.savedates; (otheroffsets{staNum} - elmayor_time);logdates(:)];
    
    %% Run nonlinear least-squares to find best-fitting parameters
    
    fit = struct('coeff_0',[],...
        'coeff_initial',[],...
        'individual_coeff',[],...
        'individual_allsignals',[],...
        'individual_allsignals_avg',zeros(length(time),3),...
        'individual_decay',[],...
        'individual_decay_int',[],...
        'individual_decay_avg',zeros(length(time),3),...
        'prefactor',[],...
        'regularized_coeff',[],...
        'regularized_allsignals',[],...
        'regularized_allsignals_avg',zeros(length(time),3),...
        'regularized_decay',[],...
        'regularized_decay_avg',zeros(length(time),3),...
        'regularized_coseismic',[],...
        'regularized_coseismic_avg',zeros(1,3),...
        'noise_coeff',[],...
        'coseismicstd',[],...
        'decaystd',[],...
        'totalstd',[]);
    
    % The regularized curve fitting used here minimizes the function 
    % (model-data).^2 + prefactor*(second derivative of cumulative decay function)
    
    % where the prefactor is a function of how "evident" the postseismic transient decay is in a given component's timeseries.
    % This "evidentness" is estimated by running a set of unregularized fits and assessing how much they agree/disagree on the postseismic decay.
    % The cumulative decay function from each unregularized fit is integrated (summed) over the postseismic period
    % and the "evidentness" is estimated by std(sums)/mean(sums).
    % The prefactor is also proportional to the average total residual of the unregularized fits (which basically penalizes
    % the second derivative in the vertical component of motion as that component is noisy and more prone to overfitting).
    % With this prefactor, the scheme penalizes the second derivative strongly for timeseries with a lot of scatter
    % and no obvious postseismic transient, preventing the curve fitting from erroneously fitting noise;
    % however it only lightly affects the solution for timeseries with an obvious postseismic signal.
    
    %MATLAB options for lsqnonlin
    lsqoptions = optimset('TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',1e4,'MaxIter',1e3,'Display','iter');
    
    %Initialize string that will become function: coeff = parameters being fit; t = time vector
    funstring  = 'coeff(1) + coeff(2)*time + coeff(3)*cos(2*pi*time) + coeff(4)*sin(2*pi*time) + coeff(5)*cos(4*pi*time) + coeff(6)*sin(4*pi*time) + coeff(7)*log(1 + heavi(time).*time/coeff(8)) + coeff(9)*(1 - exp(-(heavi(time).*time)/coeff(10))) + coeff(11)*(1 - (1 + heavi(time).*time/coeff(12)).^(-1/2.5)) + coeff(13)*heavi(time) + 0';
    
    %coeff(1) is a constant offset
    %coeff(2) is a linear trend
    %coeff(3) and coeff(4) are an annual oscillation
    %coeff(5) and coeff(6) are a semiannual oscillation
    %c(7) and c(8) are a logarithmic decay starting at the time of the earthquake (formula from Barbot et al 2009)
    %c(9) and c(10) are an exponential postseismic decay starting at the time of the earthquake (formula from Barbot et al 2009)
    %c(11) and c(12) are a power-law postseismic decay starting at the time of the earthquake (formula from Barbot et al 2009)
    %c(13) is step function for the coseismic offset
    
    %Initialize initial coefficients
    fit.coeff_0 = zeros(1,13);
    
    %Add step functions for other offsets (other earthquakes, station maintenance, etc.) to fit function
    for offsetdate = 1:length(offsets.savedates)
        fit.coeff_0(end+1) = 0;
        funstring = strrep(funstring,'+ 0', ['+ coeff(' num2str(length(fit.coeff_0)) ')*heavi(time - ' num2str(offsets.savedates(offsetdate)) ') + 0']);
    end
    
    %Finalize fit function
    funstring = strrep(funstring,'+ 0', '');
    fitfun = inline(funstring,'coeff','time');
    
    %Initialize initial coefficients
    L = 13 + length(offsets.savedates);
    
    %Bounds for coefficients
    LB = [-Inf*ones(1,6) -Inf 0.02 -Inf 0.02 -Inf 0.01 -Inf*ones(1,L-12)];
    UB = [Inf*ones(1,6) Inf 10 Inf 10 Inf 10 Inf*ones(1,L-12)];
    
    % Second-order finite difference approximation of the second derivative of the cumulative decay function
    % on a nonuniform grid (nonuniform in case data has missing days)
    secondderiv = @(coeff,time) ...
        (coeff(7)*log(1 + heavi(time(1:end-2)).*time(1:end-2)/coeff(8))...
        + coeff(9)*(1 - exp(-(heavi(time(1:end-2)).*time(1:end-2))/coeff(10)))...
        + coeff(11)*(1 - (1 + heavi(time(1:end-2)).*time(1:end-2)/coeff(12)).^(-1/2.5)))...
        *2./((time(2:end-1) - time(1:end-2)).*((time(2:end-1) - time(1:end-2)) + (time(3:end) - time(2:end-1))))...
        + (coeff(7)*log(1 + heavi(time(2:end-1)).*time(2:end-1)/coeff(8))...
        + coeff(9)*(1 - exp(-(heavi(time(2:end-1)).*time(2:end-1))/coeff(10)))...
        + coeff(11)*(1 - (1 + heavi(time(2:end-1)).*time(2:end-1)/coeff(12)).^(-1/2.5)))...
        *-2./((time(2:end-1) - time(1:end-2)).*((time(3:end) - time(2:end-1))))...
        + (coeff(7)*log(1 + heavi(time(3:end)).*time(3:end)/coeff(8))...
        + coeff(9)*(1 - exp(-(heavi(time(3:end)).*time(3:end))/coeff(10)))...
        + coeff(11)*(1 - (1 + heavi(time(3:end)).*time(3:end)/coeff(12)).^(-1/2.5)))...
        *2./((time(3:end) - time(2:end-1)).*((time(2:end-1) - time(1:end-2)) + (time(3:end) - time(2:end-1))));
    
    % Prepare sets of randomized initial values for characteristic times of decay terms
    for j = 1:searchnum
        fit.coeff_initial{j} = fit.coeff_0;
        fit.coeff_initial{j}(8) = 10.^(log10(0.02) + rand(1,1)*(log10(10) - log10(0.02)));
        fit.coeff_initial{j}(10) = 10.^(log10(0.02) + rand(1,1)*(log10(10) - log10(0.02)));
        fit.coeff_initial{j}(12) = 10.^(log10(0.01) + rand(1,1)*(log10(10) - log10(0.01)));
    end
    
    for component = 1:3
        
        % FIRST: run unregularized models to infer two things that will be used in a moment:
        % 1) the typical residual for this station and component (usually 10x higher in the vertical component than in the horizontal components).
        % This is used as the numerator in the prefactor of the damping function.
        
        % 2) a set of unregularized best-fitting decay functions (one for each set of randomized initial characteristic times).
        % The average of these decay functions will be used as a guess of the postseismic decay function
        % and summed in the denominator of the damping function. Although this is not perfectly accurate,
        % calculating the decay function on the fly while fitting the data not only takes more computational time
        % but also causes the inferred decay to become artificially large, as the function tries to maximize the denominator of
        % the smoothing function by maximizing the decay. Therefore we assess this beforehand instead.
        for j = 1:searchnum
            %Run unregularized lsqnonlin with each set of initial characteristic times
            [fit.individual_coeff{component,j}] = lsqnonlin(@(allcoeff)sqrt((fitfun(allcoeff,time) - data.data(:,component)).^2),fit.coeff_initial{j},LB,UB,lsqoptions);
            
            %The complete unregularized signal
            fit.individual_allsignals{component}(:,j) = fitfun(fit.individual_coeff{component,j},time);
            %Just the unregularized postseismic decay
            fit.individual_decay{component}(:,j) = fitfun([zeros(1,6),fit.individual_coeff{component,j}(7:12),zeros(1,L-12)],time);
            fit.individual_decay_int{component}(j) = sum(fit.individual_decay{component}(:,j));
            
            %Add to the average unregularized signal
            fit.individual_allsignals_avg(:,component) = fit.individual_allsignals_avg(:,component) + fitfun(fit.individual_coeff{component,j},time)/searchnum;
            
            %Add to the average unregularized postseismic decay
            fit.individual_decay_avg(:,component) = fit.individual_decay_avg(:,component) + fitfun([zeros(1,6),fit.individual_coeff{component,j}(7:12),zeros(1,L-12)],time)/searchnum;
        end
        
        %Calculate the prefactor (0.02 is an appropriate coefficient for this use case; a different value may be appropriate elsewhere)
        fit.prefactor(component) = 0.02*sqrt(sum((fit.individual_allsignals_avg(:,component) - data.data(:,component)).^2)/length(time))*std(fit.individual_decay_int{component})/abs(mean(fit.individual_decay_int{component}));
        
        %Run the regularized least-squares with the penalized second derivative
        for j = 1:searchnum
            [fit.regularized_coeff{component,j}] = lsqnonlin(@(allcoeff)sqrt((fitfun(allcoeff,time) - data.data(:,component)).^2 ...
                + fit.prefactor(component)*[zeros(length(time)-(length(time(time>=0))-2),1);abs(secondderiv(allcoeff,time(time>=0)))]), ...
                fit.coeff_initial{j},LB,UB,lsqoptions);
            
            %The complete regularized signal
            fit.regularized_allsignals{component}(:,j) = fitfun(fit.regularized_coeff{component,j},time);
            %Just the coseismic offset
            fit.regularized_coseismic{component}(j) = fit.regularized_coeff{component,j}(13);
            %Just the postseismic decay
            fit.regularized_decay{component}(:,j) = fitfun([zeros(1,6),fit.regularized_coeff{component,j}(7:12),zeros(1,L-12)],time);
            
            %Add to the average regularized signal
            fit.regularized_allsignals_avg(:,component) = fit.regularized_allsignals_avg(:,component) + fit.regularized_allsignals{component}(:,j)/searchnum;
            %Add to the average coseismic offset
            fit.regularized_coseismic_avg(component) = fit.regularized_coseismic_avg(component) + fit.regularized_coseismic{component}(j)/searchnum;
            %Add to the average decay function (what will get output as the north/east/up decay function)
            fit.regularized_decay_avg(:,component) = fit.regularized_decay_avg(:,component) + fit.regularized_decay{component}(:,j)/searchnum;
        end
        
        fit.noise_coseismic{component} = NaN(0,0);
        fit.noise_decay{component} = NaN(length(time),0);
        
        %To estimate the uncertainty in the fit, we add combinations of Gaussian random noise multiplied by the estimated implicit
        %uncertainty in the data (2x the SOPAC formal error) to the data and then rerun the procedure on those timeseries
        for p = 1:noisenum
            % Add noise to this component of timeseries data
            data.noise(:,component) = data.data(:,component) + randn(length(data.sigma(:,component)),1).*data.sigma(:,component);
            
            %Run the regularized least-squares with the penalized second derivative on the noisy data
            for j = 1:searchnum
                fit.noise_coeff = lsqnonlin(@(allcoeff)sqrt((fitfun(allcoeff,time) - data.noise(:,component)).^2 ...
                    + fit.prefactor(component)*[zeros(length(time)-(length(time(time>=0))-2),1);abs(secondderiv(allcoeff,time(time>=0)))]), ...
                    fit.coeff_initial{j},LB,UB,lsqoptions);
                
                %Collect just the coseismic offset
                fit.noise_coseismic{component}(end+1) = fit.noise_coeff(13);
                %Collect just the postseismic decay
                fit.noise_decay{component}(:,end+1) = fitfun([zeros(1,6),fit.noise_coeff(7:12),zeros(1,L-12)],time);
            end
        end
        
        %Estimate the total uncertainty in the coseismic offset as the standard deviation of all of the fit coseismic offsets,
        %with noise and without
        fit.coseismicstd(component) = std([fit.regularized_coseismic{component},fit.noise_coseismic{component}]);
        
        %Estimate the total uncertainty in the decay function as the standard deviation of all of the fit decay functions,
        %with noise and without
        for m = 1:length(time)
            fit.decaystd(m,component) = std([fit.regularized_decay{component}(m,:),fit.noise_decay{component}(m,:)]);
        end
        %Estimate the total uncertainty as a geometric sum of the coseismic and postseismic uncertainties
        fit.totalstd(:,component) = sqrt(fit.coseismicstd(:,component).^2 ...
            + fit.decaystd(:,component).^2);
        decay_timeseries(:,component+1) = fit.regularized_decay_avg(time>=0,component);
        decay_timeseries(:,component+4) = fit.totalstd(time>=0,component);
        
        %Prepare outputs of coseismic offsets
        coseismic_output(staNum,component) = fit.regularized_coseismic_avg(component);
        coseismic_output(staNum,component+3) = fit.coseismicstd(component);
        
        %Prepare outputs of fit postseismic decays
        if 0==(0==sum(time>4.5))
            postseismic_output(staNum,component) = fit.regularized_decay_avg(length(time(time<4.5))+1,component);
            postseismic_output(staNum,component+3) = fit.totalstd(length(time(time<4.5))+1,component);
        end
    end
    
    %% Export fits
    %=================================================================================================================================
    if options.Export_decay_timeseries
        %Export best-fit postseismic decay timeseries
        fid = fopen(['../fits/Decay_timeseries/' stations.names{staNum} '.dat'],'wt');
        fprintf(fid,'# time (yr) north(m) east(m) up(m) northsigma(m) eastsigma(m) upsigma(m)\n');
        for index = 1:length(decay_timeseries)
            fprintf(fid,'%f %f %f %f %f %f %f\n',decay_timeseries(index,:));
        end
        fclose(fid);
    end
    
    if options.Export_coseismic_offsets
        %Export fit coseismic offsets
        fid = fopen(['../fits/Coseismic/' stations.names{staNum} '.dat'],'wt');
        fprintf(fid,'# north(m) east(m) up(m)\n');
        fprintf(fid,'%f %f %f %f %f %f\n',coseismic_output(staNum,:));
        fclose(fid);
    end
    
    %% Plot timeseries
    %=================================================================================================================================
    if options.Plot_timeseries
        
        close all
        
        figure(staNum)
        
        for component = 1:3
            
            %Plot average regularized fit to all signals vs. data
            subplot(3,3,3*component-2)
            hold on
            box on
            plot(time+elmayor_time,data.data(:,component),'Color',[0.625 0.625 0.625]);
            plot(time+elmayor_time,fit.regularized_allsignals_avg(:,component),'Color','k','linewidth',1);
            
            title(['Station ' upper(stations.names{staNum}) ', ' components{component}])
            xlim([-1+elmayor_time 5+elmayor_time]);
            
            %Plot decay timeseries from unregularized, regularized and noise fits
            subplot(3,3,3*component-1)
            hold on
            box on
            
            plot(time+elmayor_time,fit.noise_decay{component},'Color',[0.625 0.625 0.625])
            plot(time+elmayor_time,fit.individual_decay{component},'Color',[1 0.25 0])
            plot(time+elmayor_time,fit.regularized_decay{component},'Color',[0 0.5 1])
            plot(time+elmayor_time,fit.regularized_decay_avg(:,component),'Color','k','linewidth',1);
            
            title(['Station ' upper(stations.names{staNum}) ', ' components{component}])
            xlim([elmayor_time 5+elmayor_time]);
            
            %Plot best-fit model + errorbars
            subplot(1,3,3)
            hold on
            box on
            
            errorbar(time(1:60:end)+elmayor_time,fit.regularized_decay_avg(1:60:end,component),fit.totalstd(1:60:end,component),'Color',[0.625 0.625 0.625]);
            plot(time+elmayor_time,fit.individual_decay_avg(:,component),'Color',[1 0 0],'linewidth',1);
            plot(time+elmayor_time,fit.regularized_decay_avg(:,component),'Color','k','linewidth',1);
            
            title(['Station ' upper(stations.names{staNum}) ', ' components{component}])
            xlim([elmayor_time 5+elmayor_time]);
            ylim([min([min(min(fit.regularized_decay_avg-fit.totalstd)),-0.02]) max([max(max(fit.regularized_decay_avg+fit.totalstd)),0.02])]);
            
            if options.Save_plots
                orient landscape
                saveas(gcf,['../fits/' stations.names{staNum} '.pdf'])
            end
            
        end
        
    end
    
end

if options.Export_cumulative_postseismic_offsets
    %Export cumulative 5-year postseismic offsets
    fexport=fopen(['../fits/postseismic_fit_5yr.dat'],'wt');
    fprintf(fexport,'# name east(km) north(km) veast(m) vnorth(m) vup(m) eastrms(m) northrms(m) uprms(m)\n');
    for staNum = start:finish
        if 0==isnan(postseismic_output(staNum,1))
            fprintf(fexport,'%s %f %f %f %f %f %f %f %f\n',stations.names{staNum},stations.east(staNum),...
                stations.north(staNum),...
                postseismic_output(staNum,2),...
                postseismic_output(staNum,1),...
                postseismic_output(staNum,3),...
                postseismic_output(staNum,5),...
                postseismic_output(staNum,4),...
                postseismic_output(staNum,6));
        end
    end
    fclose(fexport);
end

if options.Map_modeled_offsets
    %% Plot cumulative 5-year offsets in mapview
    
    figure(1000);
    hold on
    box on

    load ../plotting/uplift_colormap.mat %Uplift colormap
    
    %Plot fault, coastline and border data
    load ../plotting/ca_faults_dim.dat %Fault traces
    load ../plotting/ca_faults_km.dat %Fault traces    
    load ../plotting/ca_coast_dim.dat %Coastline
    load ../plotting/ca_coast_km.dat %Coastline
    load ../plotting/ca_border_dim.dat %Border
    load ../plotting/ca_border_km.dat %Border
    plot_faults(ca_faults_km,ca_faults_dim,xlim,ylim);
    plot_faults(ca_coast_km,ca_coast_dim,[-500 500],[-500 500]);
    plot_faults(ca_border_km,ca_border_dim,[-500 500],[-500 500]);
   
    %Plot fault trace
    ftrace = fopen(['../gmt/elmayor_trace_km.xyz']);
    trace = textscan(ftrace,'%f %f %f');
    fclose(ftrace);
    plot(trace{1},trace{2},'k','linewidth',2)
    
    %Plot GPS data
    sc=5e+2; %plotting scale
    pos=find((abs(postseismic_output(:,1))>0)); %stations with decays
    names_plot = stations.names(pos); %names of stations
    scatter(stations.east(pos),stations.north(pos),100,postseismic_output(pos,3),'filled') %Cumulative 5-year verticals
    quiver(stations.east(pos),stations.north(pos),sc*postseismic_output(pos,2),sc*postseismic_output(pos,1),0,'Color',[0.5 0 0.5],'linewidth',1); %5-year horizontals
    
    ellipse(sc*postseismic_output(pos,5),sc*postseismic_output(pos,4),zeros((length(pos)),1),...
        stations.east(pos) + sc*postseismic_output(pos,2),stations.north(pos) + sc*postseismic_output(pos,1),'k'); %Uncertainties in horizontals
    text(stations.east(pos),stations.north(pos),names_plot,'FontSize',8); %names of stations
    
    %Scale arrow
    text(-130,-45,'5 cm')
    quiver(-130,-50,sc*(5e-2),0,0,'k','linewidth',1)
    
    xlim=[-187.4 93.7];
    ylim=[-55.5 166.7];
    axis equal
    
    set(gca,'xlim',xlim,'ylim',ylim,'clim',[-0.02 0.02])
    colormap(scec_colormap);
    colorbar('YTickLabel',{'-2.0 cm','-1.5 cm','-1.0 cm','-0.5 cm','0.0 cm','+0.5 cm','+1.0 cm','+1.5 cm','+2.0 cm'})
    
    if options.Save_maps
        orient landscape
        saveas(gcf,['../fits/postseis_horiz_Mar252015.pdf'])
    end 
end