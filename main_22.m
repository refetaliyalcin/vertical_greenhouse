clc
clear all
close all

economizer = 0; % 1 or 0

% Case 1 - IR ✔ Flo ✔ CD ✔ LED ✔ % All In
% Case 2 - IR X Flo ✔ CD ✔ LED ✔ % Investigate effect of IR
% Case 3 - IR ✔ Flo X CD ✔ LED ✔ % Investigate effect of Flo
% Case 4 - IR X Flo X CD X LED ✔ % LED Only
% Case 5 - IR ✔ Flo ✔ CD ✔ LED X % Sun Only

%design parameters starts%

solar_collector_area = 120; % at the roof in m^2
greenhouse_area = 77; % from top view in m^2
hour_light = 20; %20h / 24h lighting will be used according to https://doi.org/10.1626/JCS.58.689
number_of_shelves = round(8*120/77); %subject to change to find an optimum
desired_PPFD_per_shelf = 120; % in mmol*m^-2*s^-1 we can consider to adjust 
ita_cd_light = 0.5; %collection distribution efficiency
C_thermal = 0.7; % i think it should be smaller than ita_cd_light as absorption in fiber cables also contributes to heating beside transmitted radiation
t_eff=1.13; %fluorescent effect, effective transmittance of the fluorescent coating
PPF_W_coeff = 3.3; % PPFD to wall-plug Watt conversion of growth LEDs https://www.assets.signify.com/is/content/Signify/Assets/philips-lighting/global/20211217-production-module.pdf
doy_to_on_IR_filter=-1; %day of the year to enable IR filter. make -1 to enable always, 90 is approximately 1 April subject to change wrt climate
doy_to_off_IR_filter=270; %day of the year to disable IR filter  270 is approximately 1 October subject to change wrt climate

%design parameters ends%

load('solar_data_22.mat')
ir_filter
lamda=(280:4000)';
IR_filter=interp1(IR_filter_raw(:,1),IR_filter_raw(:,2),lamda);
% solar_data(wavelength,tip,doy,hour)
% tip (1=direct_horiz,2=diff_horiz,3=direct_PPFD)
par_start_ind=400-lamda(1)+1;
par_end_ind=700-lamda(1)+1;
direct_solar=zeros(365,24);
direct_solar_IR_filter=zeros(365,24);
PPFD=zeros(365,24);
PPFD_IR_filter=zeros(365,24);
growth_day_hour_solar_only=zeros(365,24);
for doy=1:365
    for hour=1:24
        % when to use IR filter starts
        if doy_to_on_IR_filter > 0
            if doy >= doy_to_on_IR_filter && doy <= doy_to_off_IR_filter
               IR_filter_used = IR_filter;
            else
               IR_filter_used = ones(length(IR_filter),1);
            end
        else
            IR_filter_used = IR_filter; % always enabled
        end
        % when to use IR filter ends
        direct_solar_lamda=solar_data(:,1,doy,hour);
        filtered_direct_solar_lamda=IR_filter_used.*direct_solar_lamda;
        direct_par_lamda=solar_data(:,2,doy,hour);
        filtered_direct_par_lamda=direct_par_lamda.*IR_filter_used;
        direct_solar(doy,hour)=trapz(lamda,direct_solar_lamda);
        direct_solar_IR_filter(doy,hour)=trapz(lamda,filtered_direct_solar_lamda);
        PPFD(doy,hour)=trapz(lamda(par_start_ind:par_end_ind),direct_par_lamda(par_start_ind:par_end_ind));
        PPFD_IR_filter(doy,hour)=trapz(lamda(par_start_ind:par_end_ind),filtered_direct_par_lamda(par_start_ind:par_end_ind));
    end
end

% PPF Calculation Starts
PPF_case_1 = PPFD_IR_filter * ita_cd_light * solar_collector_area * t_eff; % cases 1 and 5
PPF_case_2 = PPFD * ita_cd_light * solar_collector_area * t_eff; % case 2
PPF_case_3 = PPFD_IR_filter * ita_cd_light * solar_collector_area; % case 3
PPF_case_5 = PPF_case_1;

total_PPF_desired = desired_PPFD_per_shelf * number_of_shelves * greenhouse_area; % cases 1, 2, 3, 4

LED_compansate_case_1 = total_PPF_desired - PPF_case_1; % find the amount to be compansated by LEDs
LED_compansate_case_1(LED_compansate_case_1<0)=0; % can't compansate negatives, so make them zero
LED_compansate_case_1(:,hour_light+1:end)=0; %close the lights after 20h

LED_compansate_case_2 = total_PPF_desired - PPF_case_2; % find the amount to be compansated by LEDs
LED_compansate_case_2(LED_compansate_case_2<0)=0; % can't compansate negatives, so make them zero
LED_compansate_case_2(:,hour_light+1:end)=0; %close the lights after 20h

LED_compansate_case_3 = total_PPF_desired - PPF_case_3; % find the amount to be compansated by LEDs
LED_compansate_case_3(LED_compansate_case_3<0)=0; % can't compansate negatives, so make them zero
LED_compansate_case_3(:,hour_light+1:end)=0; %close the lights after 20h

LED_compansate_case_4 = total_PPF_desired*ones(365,24); % find the amount to be compansated by LEDs
LED_compansate_case_4(:,hour_light+1:end)=0; %close the lights after 20h

LED_compansate_case_5 = zeros(365,24); % find the amount to be compansated by LEDs

% PAR Calculation Ends

%Growth Calculation starts%

for doy=1:365
    for hour=1:24
        growth_day_hour_solar_only(doy,hour)=growth_fn_2(t_eff*PPF_case_5(doy,hour)/(number_of_shelves*greenhouse_area))*greenhouse_area*number_of_shelves;
    end
end

yearly_lettuce_kg_cases_1_2_3_4 = growth_fn_2(desired_PPFD_per_shelf)*greenhouse_area*number_of_shelves*365*hour_light/1000 % Case 1,2,3,4

yearly_lettuce_kg_case_5 = sum(growth_day_hour_solar_only,'All')/1000 % Case 5


%Growth calculation ends%

% Heat Load Calculations Starts
Q_solar = direct_solar * C_thermal * solar_collector_area; %in Watt Case 2
Q_solar_IR_filter = direct_solar_IR_filter * C_thermal * solar_collector_area; %in Watt Cases 1, 3 and 5
Q_solar_case_1 = Q_solar_IR_filter;
Q_solar_case_2 = Q_solar;
Q_solar_case_3 = Q_solar_IR_filter;
Q_solar_case_4 = zeros(365,24);
Q_solar_case_5 = Q_solar_IR_filter;

Q_led_case_1 = LED_compansate_case_1 / PPF_W_coeff; %in Watt
Q_led_case_2 = LED_compansate_case_2 / PPF_W_coeff; %in Watt
Q_led_case_3 = LED_compansate_case_3 / PPF_W_coeff; %in Watt
Q_led_case_4 = LED_compansate_case_4 / PPF_W_coeff; %in Watt
Q_led_case_5 = zeros(365,24);

Q_total_load_case_1 = Q_solar_case_1 + Q_led_case_1; %in Watt
Q_total_load_case_2 = Q_solar_case_2 + Q_led_case_2; %in Watt
Q_total_load_case_3 = Q_solar_case_3 + Q_led_case_3; %in Watt
Q_total_load_case_4 = Q_solar_case_4 + Q_led_case_4; %in Watt
Q_total_load_case_5 = Q_solar_case_5 + Q_led_case_5; %in Watt

maximum_q_1=max(max(Q_total_load_case_1))
maximum_q_2=max(max(Q_total_load_case_2))
maximum_q_3=max(max(Q_total_load_case_3))
maximum_q_4=max(max(Q_total_load_case_4))
maximum_q_5=max(max(Q_total_load_case_5))

fileID = fopen('idf_output_case_1.txt','w');
for i=1:365
    fprintf(fileID,'\n\nSchedule:Day:Interval,\n');
    fprintf(fileID,'    day%d,                   !- Name\n',i);
    fprintf(fileID,'    Fraction,                !- Schedule Type Limits Name\n');
    fprintf(fileID,'    No,                      !- Interpolate to Timestep\n');
    for j=1:23
        ratio=Q_total_load_case_1(i,j)/maximum_q_1;
        fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',j,j);
        fprintf(fileID,'    %.5f,                 !- Value Until Time %d\n',ratio,j);
    end
    ratio=Q_total_load_case_1(i,24)/maximum_q_1;
    fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',24,24);
    fprintf(fileID,'    %.5f;                 !- Value Until Time %d\n',ratio,24);
end
fclose(fileID);

fileID = fopen('idf_output_case_2.txt','w');
for i=1:365
    fprintf(fileID,'\n\nSchedule:Day:Interval,\n');
    fprintf(fileID,'    day%d,                   !- Name\n',i);
    fprintf(fileID,'    Fraction,                !- Schedule Type Limits Name\n');
    fprintf(fileID,'    No,                      !- Interpolate to Timestep\n');
    for j=1:23
        ratio=Q_total_load_case_2(i,j)/maximum_q_2;
        fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',j,j);
        fprintf(fileID,'    %.5f,                 !- Value Until Time %d\n',ratio,j);
    end
    ratio=Q_total_load_case_2(i,24)/maximum_q_2;
    fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',24,24);
    fprintf(fileID,'    %.5f;                 !- Value Until Time %d\n',ratio,24);
end
fclose(fileID);

fileID = fopen('idf_output_case_3.txt','w');
for i=1:365
    fprintf(fileID,'\n\nSchedule:Day:Interval,\n');
    fprintf(fileID,'    day%d,                   !- Name\n',i);
    fprintf(fileID,'    Fraction,                !- Schedule Type Limits Name\n');
    fprintf(fileID,'    No,                      !- Interpolate to Timestep\n');
    for j=1:23
        ratio=Q_total_load_case_3(i,j)/maximum_q_3;
        fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',j,j);
        fprintf(fileID,'    %.5f,                 !- Value Until Time %d\n',ratio,j);
    end
    ratio=Q_total_load_case_3(i,24)/maximum_q_3;
    fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',24,24);
    fprintf(fileID,'    %.5f;                 !- Value Until Time %d\n',ratio,24);
end
fclose(fileID);

fileID = fopen('idf_output_case_4.txt','w');
for i=1:365
    fprintf(fileID,'\n\nSchedule:Day:Interval,\n');
    fprintf(fileID,'    day%d,                   !- Name\n',i);
    fprintf(fileID,'    Fraction,                !- Schedule Type Limits Name\n');
    fprintf(fileID,'    No,                      !- Interpolate to Timestep\n');
    for j=1:23
        ratio=Q_total_load_case_4(i,j)/maximum_q_4;
        fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',j,j);
        fprintf(fileID,'    %.5f,                 !- Value Until Time %d\n',ratio,j);
    end
    ratio=Q_total_load_case_4(i,24)/maximum_q_4;
    fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',24,24);
    fprintf(fileID,'    %.5f;                 !- Value Until Time %d\n',ratio,24);
end
fclose(fileID);

fileID = fopen('idf_output_case_5.txt','w');
for i=1:365
    fprintf(fileID,'\n\nSchedule:Day:Interval,\n');
    fprintf(fileID,'    day%d,                   !- Name\n',i);
    fprintf(fileID,'    Fraction,                !- Schedule Type Limits Name\n');
    fprintf(fileID,'    No,                      !- Interpolate to Timestep\n');
    for j=1:23
        ratio=Q_total_load_case_5(i,j)/maximum_q_5;
        fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',j,j);
        fprintf(fileID,'    %.5f,                 !- Value Until Time %d\n',ratio,j);
    end
    ratio=Q_total_load_case_5(i,24)/maximum_q_5;
    fprintf(fileID,'    %02d:00,                   !- Time %d {hh:mm}\n',24,24);
    fprintf(fileID,'    %.5f;                 !- Value Until Time %d\n',ratio,24);
end
fclose(fileID);


fileID = fopen('idf_output_common.txt','w');
for i=1:365
    [yy,mm,dd,HH,MM] = datevec(datenum(2012,1,doy));
    if mod(i,7)==1
        fprintf(fileID,'\nSchedule:Week:Daily,\n');
        fprintf(fileID,'Jan_Week%d,               !- Name\n',(i+6)/7);
    end
    if mod(i,7)==1 fprintf(fileID,'    day%d,                    !- Sunday Schedule:Day Name\n',i);, end
    if mod(i,7)==2 fprintf(fileID,'    day%d,                    !- Monday Schedule:Day Name\n',i);, end
    if mod(i,7)==3 fprintf(fileID,'    day%d,                    !- Tuesday Schedule:Day Name\n',i);, end
    if mod(i,7)==4 fprintf(fileID,'    day%d,                    !- Wednesday Schedule:Day Name\n',i);, end
    if mod(i,7)==5 fprintf(fileID,'    day%d,                    !- Thursday Schedule:Day Name\n',i);, end
    if mod(i,7)==6 fprintf(fileID,'    day%d,                    !- Friday Schedule:Day Name\n',i);, end
    if mod(i,7)==0 fprintf(fileID,'    day%d,                    !- Saturday Schedule:Day Name\n',i);, end
    if mod(i,7)==0
    fprintf(fileID,'    day%d,                    !- Holiday Schedule:Day Name\n',i-6);
    fprintf(fileID,'    day%d,                    !- SummerDesignDay Schedule:Day Name\n',i-6);
    fprintf(fileID,'    day%d,                    !- WinterDesignDay Schedule:Day Name\n',i-6);
    fprintf(fileID,'    day%d,                    !- CustomDay1 Schedule:Day Name\n',i-6);
    fprintf(fileID,'    day%d;                    !- CustomDay2 Schedule:Day Name\n',i-6);
    end
end
fclose(fileID);

Q_led_monthly_case_1=zeros(1,12);
Q_led_monthly_case_2=zeros(1,12);
Q_led_monthly_case_3=zeros(1,12);
Q_led_monthly_case_4=zeros(1,12);
Q_led_monthly_case_5=zeros(1,12);
Q_solar_monthly_case_1=zeros(1,12);
Q_solar_monthly_case_2=zeros(1,12);
Q_solar_monthly_case_3=zeros(1,12);
Q_solar_monthly_case_4=zeros(1,12);
Q_solar_monthly_case_5=zeros(1,12);
for i=1:365
   t=datetime(2012,1,i);
    [ y , m , d ] = ymd( t );
    Q_led_monthly_case_1(m)=Q_led_monthly_case_1(m)+sum(Q_led_case_1(i,:));
    Q_led_monthly_case_2(m)=Q_led_monthly_case_2(m)+sum(Q_led_case_2(i,:));
    Q_led_monthly_case_3(m)=Q_led_monthly_case_3(m)+sum(Q_led_case_3(i,:));
    Q_led_monthly_case_4(m)=Q_led_monthly_case_4(m)+sum(Q_led_case_4(i,:));
    Q_led_monthly_case_5(m)=Q_led_monthly_case_5(m)+sum(Q_led_case_5(i,:));
    Q_solar_monthly_case_1(m)=Q_solar_monthly_case_1(m)+sum(Q_solar_case_1(i,:));
    Q_solar_monthly_case_2(m)=Q_solar_monthly_case_2(m)+sum(Q_solar_case_2(i,:));
    Q_solar_monthly_case_3(m)=Q_solar_monthly_case_3(m)+sum(Q_solar_case_3(i,:));
    Q_solar_monthly_case_4(m)=Q_solar_monthly_case_4(m)+sum(Q_solar_case_4(i,:));
    Q_solar_monthly_case_5(m)=Q_solar_monthly_case_5(m)+sum(Q_solar_case_5(i,:));
end



% Yearly_solar_1 = sum(Q_solar_case_1,'all') %this is yearly total we normally do not need it
% Yearly_solar_2 = sum(Q_solar_case_2,'all')
% Yearly_solar_3 = sum(Q_solar_case_3,'all')
% Yearly_solar_4 = sum(Q_solar_case_4,'all')
% Yearly_solar_5 = sum(Q_solar_case_5,'all')
% 
Yearly_led_1 = sum(Q_led_case_1,'all')
Yearly_led_2 = sum(Q_led_case_2,'all')
Yearly_led_3 = sum(Q_led_case_3,'all')
Yearly_led_4 = sum(Q_led_case_4,'all')
Yearly_led_5 = sum(Q_led_case_5,'all')
% 
% Yearly_total_1 = Yearly_solar_1 + Yearly_led_1;
% Yearly_total_2 = Yearly_solar_2 + Yearly_led_2;
% Yearly_total_3 = Yearly_solar_3 + Yearly_led_3;
% Yearly_total_4 = Yearly_solar_4 + Yearly_led_4;
% Yearly_total_5 = Yearly_solar_5 + Yearly_led_5;

% eff_1 = yearly_lettuce_kg_case_1_2_3_4 / Yearly_total_1 % this is not the eff. we seek for
% eff_2 = yearly_lettuce_kg_case_1_2_3_4 / Yearly_total_2
% eff_3 = yearly_lettuce_kg_case_1_2_3_4 / Yearly_total_3
% eff_4 = yearly_lettuce_kg_case_1_2_3_4 / Yearly_total_4
% eff_5 = yearly_lettuce_kg_case_5 / Yearly_total_5


% Heat Load Calculations Ends


save('solar_load_data.mat','Q_solar_case_1','Q_solar_case_2','Q_solar_case_3','Q_solar_case_4','Q_solar_case_5')
save('led_load_data.mat','Q_led_case_1','Q_led_case_2','Q_led_case_3','Q_led_case_4','Q_led_case_5')
save('Q_max_for_idf.mat','maximum_q_1','maximum_q_2','maximum_q_3','maximum_q_4','maximum_q_5')

formatSpec = '%.0f';
if economizer == 1
temptext = fileread('../Case_X_temp_v24.1_economizer.idf');
else
temptext = fileread('../Case_X_temp_v24.1_first.idf');
end
text1 = strrep(temptext,'change_here',fileread('idf_output_case_1.txt'));
text1 = strrep(text1,'q_max_number',num2str(maximum_q_1,formatSpec));
text2 = strrep(temptext,'change_here',fileread('idf_output_case_2.txt'));
text2 = strrep(text2,'q_max_number',num2str(maximum_q_2,formatSpec));
text3 = strrep(temptext,'change_here',fileread('idf_output_case_3.txt'));
text3 = strrep(text3,'q_max_number',num2str(maximum_q_3,formatSpec));
text4 = strrep(temptext,'change_here',fileread('idf_output_case_4.txt'));
text4 = strrep(text4,'q_max_number',num2str(maximum_q_4,formatSpec));
text5 = strrep(temptext,'change_here',fileread('idf_output_case_5.txt'));
text5 = strrep(text5,'q_max_number',num2str(maximum_q_5,formatSpec));



fileID = fopen('Case_1_22.idf','w');
fprintf(fileID,'%s',text1);
fclose(fileID);
fileID = fopen('Case_2_22.idf','w');
fprintf(fileID,'%s',text2);
fclose(fileID);
fileID = fopen('Case_3_22.idf','w');
fprintf(fileID,'%s',text3);
fclose(fileID);
fileID = fopen('Case_4_22.idf','w');
fprintf(fileID,'%s',text4);
fclose(fileID);
fileID = fopen('Case_5_22.idf','w');
fprintf(fileID,'%s',text5);
fclose(fileID);

error('stop')


figure
contourf([1:20,20.1,22,23,24],1:365,LED_compansate_case_1/1000,(0:0.1:120),'edgecolor','none')
cb1 = colorbar;
hAx=gca;
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
hLg.EdgeColor = [0 0 0];
xlabel('Hour, h [h]')
ylh=ylabel('Day of the year');
ylh.VerticalAlignment	= 'bottom'; %if it is not alligned well, try 'top' and 'bottom' too
cb1.Label.String = 'LED Power [kW]';
cb1.Label.Position(1) = 3;
caxis([0 120])
hold on;
[M,c] = contour(LED_compansate_case_1/1000,[10,50,100],'-k',"ShowText",true)
c.LineWidth = 3;
xticks([1 6 12 18 24])
yticks([1 60 120 180 240 300 365])
set(gca,'FontSize',13)
saveas(gcf,'22_fig_6b.png')


figure
contourf(PPFD,(0:1:2000),'edgecolor','none')
cb1 = colorbar;
hAx=gca;
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
hLg.EdgeColor = [0 0 0];
xlabel('Hour, h [h]')
ylh=ylabel('Day of the year');
ylh.VerticalAlignment	= 'bottom'; %if it is not alligned well, try 'top' and 'bottom' too
cb1.Label.String = 'PPFD [\mumol m^-^2 s^-^1]';
cb1.Label.Position(1) = 3;
caxis([0 2000])
hold on;
[M,c] = contour(PPFD,[100,1000,1500,1800,1900],'-k',"ShowText",true);
c.LineWidth = 3;
xticks([1 6 12 18 24])
yticks([1 60 120 180 240 300 365])
set(gca,'FontSize',13)
saveas(gcf,'22_fig_6a.png')

% error('stop')


ep_result_1 = call_ep('Case_1_22','CHN_Hong.Kong.SAR.450070_CityUHK');
ep_result_2 = call_ep('Case_2_22','CHN_Hong.Kong.SAR.450070_CityUHK');
ep_result_3 = call_ep('Case_3_22','CHN_Hong.Kong.SAR.450070_CityUHK');
ep_result_4 = call_ep('Case_4_22','CHN_Hong.Kong.SAR.450070_CityUHK');
ep_result_5 = call_ep('Case_5_22','CHN_Hong.Kong.SAR.450070_CityUHK');


cooling_1_ton_hrs=0;
month_cooling_1_ton_hrs = zeros(12,1);




if economizer == 1
    start=68;
else
    start=14;
end

for i=start:start+8
    cooling_1_ton_hrs = cooling_1_ton_hrs + ep_result_1{i}{2}{15}{2};
    for j=1:12
        month_cooling_1_ton_hrs(j)=month_cooling_1_ton_hrs(j)+ ep_result_1{i}{2}{j+1}{2};
    end
end
heating_1_kBtu=0;
month_heating_1_kBtu = zeros(12,1);
for i=start+9:start+17
    heating_1_kBtu = heating_1_kBtu + ep_result_1{i}{2}{15}{2};
    for j=1:12
        month_heating_1_kBtu(j)=month_heating_1_kBtu(j)+ ep_result_1{i}{2}{j+1}{2};
    end
end
electricity_1_kWh=0;
for i=start+18:start+26
    electricity_1_kWh = electricity_1_kWh + ep_result_1{i}{2}{15}{2};
end

cooling_2_ton_hrs=0;
month_cooling_2_ton_hrs = zeros(12,1);
for i=start:start+8
    cooling_2_ton_hrs = cooling_2_ton_hrs + ep_result_2{i}{2}{15}{2};
    for j=1:12
        month_cooling_2_ton_hrs(j)=month_cooling_2_ton_hrs(j)+ ep_result_2{i}{2}{j+1}{2};
    end
end
heating_2_kBtu=0;
month_heating_2_kBtu = zeros(12,1);
for i=start+9:start+17
    heating_2_kBtu = heating_2_kBtu + ep_result_2{i}{2}{15}{2};
    for j=1:12
        month_heating_2_kBtu(j)=month_heating_2_kBtu(j)+ ep_result_2{i}{2}{j+1}{2};
    end
end
electricity_2_kWh=0;
for i=start+18:start+26
    electricity_2_kWh = electricity_2_kWh + ep_result_2{i}{2}{15}{2};
end

cooling_3_ton_hrs=0;
month_cooling_3_ton_hrs = zeros(12,1);
for i=start:start+8
    cooling_3_ton_hrs = cooling_3_ton_hrs + ep_result_3{i}{2}{15}{2};
    for j=1:12
        month_cooling_3_ton_hrs(j)=month_cooling_3_ton_hrs(j)+ ep_result_3{i}{2}{j+1}{2};
    end
end
heating_3_kBtu=0;
month_heating_3_kBtu = zeros(12,1);
for i=start+9:start+17
    heating_3_kBtu = heating_3_kBtu + ep_result_3{i}{2}{15}{2};
    for j=1:12
        month_heating_3_kBtu(j)=month_heating_3_kBtu(j)+ ep_result_3{i}{2}{j+1}{2};
    end
end
electricity_3_kWh=0;
for i=start+18:start+26
    electricity_3_kWh = electricity_3_kWh + ep_result_3{i}{2}{15}{2};
end

cooling_4_ton_hrs=0;
month_cooling_4_ton_hrs = zeros(12,1);
for i=start:start+8
    cooling_4_ton_hrs = cooling_4_ton_hrs + ep_result_4{i}{2}{15}{2};
    for j=1:12
        month_cooling_4_ton_hrs(j)=month_cooling_4_ton_hrs(j)+ ep_result_4{i}{2}{j+1}{2};
    end
end
heating_4_kBtu=0;
month_heating_4_kBtu = zeros(12,1);
for i=start+9:start+17
    heating_4_kBtu = heating_4_kBtu + ep_result_4{i}{2}{15}{2};
    for j=1:12
        month_heating_4_kBtu(j)=month_heating_4_kBtu(j)+ ep_result_4{i}{2}{j+1}{2};
    end
end
electricity_4_kWh=0;
for i=start+18:start+26
    electricity_4_kWh = electricity_4_kWh + ep_result_4{i}{2}{15}{2};
end

cooling_5_ton_hrs=0;
month_cooling_5_ton_hrs = zeros(12,1);
for i=start:start+8
    cooling_5_ton_hrs = cooling_5_ton_hrs + ep_result_5{i}{2}{15}{2};
    for j=1:12
        month_cooling_5_ton_hrs(j)=month_cooling_5_ton_hrs(j)+ ep_result_5{i}{2}{j+1}{2};
    end
end
heating_5_kBtu=0;
month_heating_5_kBtu = zeros(12,1);
for i=start+9:start+17
    heating_5_kBtu = heating_5_kBtu + ep_result_5{i}{2}{15}{2};
    for j=1:12
        month_heating_5_kBtu(j)=month_heating_5_kBtu(j)+ ep_result_5{i}{2}{j+1}{2};
    end
end
electricity_5_kWh=0;
for i=start+18:start+26
    electricity_5_kWh = electricity_5_kWh + ep_result_5{i}{2}{15}{2};
end

cooling_1_MWh=cooling_1_ton_hrs*3.5168528421/1000;
cooling_2_MWh=cooling_2_ton_hrs*3.5168528421/1000;
cooling_3_MWh=cooling_3_ton_hrs*3.5168528421/1000;
cooling_4_MWh=cooling_4_ton_hrs*3.5168528421/1000;
cooling_5_MWh=cooling_5_ton_hrs*3.5168528421/1000;
heating_1_MWh=heating_1_kBtu*0.29301453352086/1000;
heating_2_MWh=heating_2_kBtu*0.29301453352086/1000;
heating_3_MWh=heating_3_kBtu*0.29301453352086/1000;
heating_4_MWh=heating_4_kBtu*0.29301453352086/1000;
heating_5_MWh=heating_5_kBtu*0.29301453352086/1000;
cooling_1_MWh=ep_result_1{4}{2}{17}{12}*0.29301453352086/1000;
heating_1_MWh=ep_result_1{4}{2}{17}{13}*0.29301453352086/1000;
cooling_2_MWh=ep_result_2{4}{2}{17}{12}*0.29301453352086/1000;
heating_2_MWh=ep_result_2{4}{2}{17}{13}*0.29301453352086/1000;
cooling_3_MWh=ep_result_3{4}{2}{17}{12}*0.29301453352086/1000;
heating_3_MWh=ep_result_3{4}{2}{17}{13}*0.29301453352086/1000;
cooling_4_MWh=ep_result_4{4}{2}{17}{12}*0.29301453352086/1000;
heating_4_MWh=ep_result_4{4}{2}{17}{13}*0.29301453352086/1000;
cooling_5_MWh=ep_result_5{4}{2}{17}{12}*0.29301453352086/1000;
heating_5_MWh=ep_result_5{4}{2}{17}{13}*0.29301453352086/1000;
led_1_MWh=Yearly_led_1/10^6;
led_2_MWh=Yearly_led_2/10^6;
led_3_MWh=Yearly_led_3/10^6;
led_4_MWh=Yearly_led_4/10^6;
led_5_MWh=Yearly_led_5/10^6;
total_1=cooling_1_MWh+heating_1_MWh+led_1_MWh;
total_2=cooling_2_MWh+heating_2_MWh+led_2_MWh;
total_3=cooling_3_MWh+heating_3_MWh+led_3_MWh;
total_4=cooling_4_MWh+heating_4_MWh+led_4_MWh;
total_5=cooling_5_MWh+heating_5_MWh+led_5_MWh;
En_per_crop_1=1000*total_1/yearly_lettuce_kg_cases_1_2_3_4;
En_per_crop_2=1000*total_2/yearly_lettuce_kg_cases_1_2_3_4;
En_per_crop_3=1000*total_3/yearly_lettuce_kg_cases_1_2_3_4;
En_per_crop_4=1000*total_4/yearly_lettuce_kg_cases_1_2_3_4;
En_per_crop_5=1000*total_5/yearly_lettuce_kg_case_5;


Type = ["Heating load (MWh)";"Cooling load (MWh)";"LED electrical load (MWh)";"Total load (MWh)";"Growth (ton)";"Energy consumption per kg lettuce (kWh/kg)";"Change in energy cons."];
Case_1=[heating_1_MWh;cooling_1_MWh;led_1_MWh;total_1;yearly_lettuce_kg_cases_1_2_3_4/1000;En_per_crop_1;"-"];
Case_2=[heating_2_MWh;cooling_2_MWh;led_2_MWh;total_2;yearly_lettuce_kg_cases_1_2_3_4/1000;En_per_crop_2;100*(En_per_crop_2-En_per_crop_1)/En_per_crop_1];
Case_3=[heating_3_MWh;cooling_3_MWh;led_3_MWh;total_3;yearly_lettuce_kg_cases_1_2_3_4/1000;En_per_crop_3;100*(En_per_crop_3-En_per_crop_1)/En_per_crop_1];
Case_4=[heating_4_MWh;cooling_4_MWh;led_4_MWh;total_4;yearly_lettuce_kg_cases_1_2_3_4/1000;En_per_crop_4;100*(En_per_crop_4-En_per_crop_1)/En_per_crop_1];
Case_5=[heating_5_MWh;cooling_5_MWh;led_5_MWh;total_5;yearly_lettuce_kg_case_5/1000;En_per_crop_5;100*(En_per_crop_5-En_per_crop_1)/En_per_crop_1];
table_5= table(Type,Case_1,Case_2,Case_3,Case_4,Case_5)



Type = ["Jan";"Feb";"Mar";"Apr";"May";"Jun";"Jul";"Aug";"Sep";"Oct";"Nov";"Dec"];
Case_1=[month_cooling_1_ton_hrs*3.5168528421/1000];
Case_2=[month_cooling_2_ton_hrs*3.5168528421/1000];
Case_3=[month_cooling_3_ton_hrs*3.5168528421/1000];
Case_4=[month_cooling_4_ton_hrs*3.5168528421/1000];
Case_5=[month_cooling_5_ton_hrs*3.5168528421/1000];
fig_7_cooling= table(Type,Case_1,Case_2,Case_3,Case_4,Case_5)


Type = ["Jan";"Feb";"Mar";"Apr";"May";"Jun";"Jul";"Aug";"Sep";"Oct";"Nov";"Dec"];
Case_1=[month_heating_1_kBtu*0.29301453352086/1000];
Case_2=[month_heating_2_kBtu*0.29301453352086/1000];
Case_3=[month_heating_3_kBtu*0.29301453352086/1000];
Case_4=[month_heating_4_kBtu*0.29301453352086/1000];
Case_5=[month_heating_5_kBtu*0.29301453352086/1000];
fig_7_heating= table(Type,Case_1,Case_2,Case_3,Case_4,Case_5)


Type = ["Jan";"Feb";"Mar";"Apr";"May";"Jun";"Jul";"Aug";"Sep";"Oct";"Nov";"Dec"];
Case_1=[Q_led_monthly_case_1'/10^6];
Case_2=[Q_led_monthly_case_2'/10^6];
Case_3=[Q_led_monthly_case_3'/10^6];
Case_4=[Q_led_monthly_case_4'/10^6];
Case_5=[Q_led_monthly_case_5'/10^6];
fig_7_led= table(Type,Case_1,Case_2,Case_3,Case_4,Case_5)

