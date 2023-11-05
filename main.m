clc
clear all
close all

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

load('solar_data.mat')
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
        direct_par_lamda=solar_data(:,3,doy,hour);
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

maximum_q_1=max(max(Q_total_load_case_1));
fileID = fopen('deneme.txt','w');
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

for i=1:365
    [yy,mm,dd,HH,MM] = datevec(datenum(2012,1,doy));
    if mod(i,7)==1
        fprintf(fileID,'\n\nSchedule:Week:Daily,\n');
        fprintf(fileID,'Jan_Week%d,               !- Name\n',(i+6)/7);
    end
    fprintf(fileID,'    day%d,                    !- Sunday Schedule:Day Name\n',i);
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
% hold on;
% [M,c] = contour(LED_compansate_case_1/1000,[10,50,100],'-k',"ShowText",true)
% c.LineWidth = 3;
xticks([1 6 12 18 24])
yticks([1 60 120 180 240 300 365])
set(gca,'FontSize',13)
saveas(gcf,'fig_6b_.png')


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
[M,c] = contour(PPFD,[100,1000,1500,1800],'-k',"ShowText",true);
c.LineWidth = 3;
xticks([1 6 12 18 24])
yticks([1 60 120 180 240 300 365])
set(gca,'FontSize',13)
saveas(gcf,'fig_6a.png')




% figure
% contourf(direct_solar_IR_filter)
% cb2 = colorbar;
% figure
% contourf(PAR)
% cb3 = colorbar;
% figure
% contourf(PAR_IR_filter)
% cb4 = colorbar;