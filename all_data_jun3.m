% all_data_jun3.m

% take all files and make a single datafile containing all data adding the
% arms' lengths

% new version: takes measurements with small difference between pressures
% by different methods
mainfolder = cd;
%junA3 = load('junA3.mat');    % load datafile
%junA3 = junA3.junA3;

tolerance = 0.02;   %% 10% difference between pressure in different methods is allowed

junctionfiles = dir('*junction*');
dirFlags = [junctionfiles.isdir];
junctionfolders = junctionfiles(dirFlags);

%{
files = dir;
dirFlags = [files.isdir];
folders = files(dirFlags);
tubefiles = dir('*tube*');
dirFlags = [tubefiles.isdir];
tubefolders = tubefiles(dirFlags);
%}
a = [];

%% adds arm length?
for ind = 1 : size(junctionfolders,1)
    %folder = junctionfolders(ind).name.
    cd(junctionfolders(ind).name)
    files_wanted = dir('junction with arm length =*.txt');
    for ind2 = 1 : size(files_wanted,1)
        datafile = files_wanted(ind2).name;
        atmp = load(datafile);
        % find the arms' length in datafile name
        C = strsplit(datafile);
        k = strfind(C,"=");
        atmp(:,8) = str2num(C{find(not(cellfun('isempty',k)),1)+1});  %% adds arms' lengths

        a = [a;atmp];
    end
cd(mainfolder);   
%[datafile, folder] = uigetfile('*.txt');%,'MultiSelect','on');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the pressure using the tube's radius
%dEdr = 2*P*r.^3-2*gamma*r.^2+kappa ;
gamma = 1; kappa = 1; 
r = zeros(1,length(0.0001 : 0.0001 : ((2/3)*(gamma/kappa^(1/3)))^1.5)) ;
ind = 1;
for P = 0.0001 : 0.0001 : ((2/3)*(gamma/kappa^(1/3)))^1.5 
    nr  = roots([2*P -2*gamma 0 kappa]);
    r(ind) = nr(2);
    ind = ind+1;
end
P = 0.0001 : 0.0001 : ((2/3)*(gamma/kappa^(1/3)))^1.5 ;
a(:,7) = interp1(r,P,a(:,5)./sqrt(kappa./a(:,2))).*a(:,2).^(1.5)/(kappa^0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% copy data
%%
A =[];
a(isnan(a)) = 0;  %% switch NaN (from interpolation) to zeros
a(abs(2*(a(:,6)-a(:,7))./(a(:,6)+a(:,7))) > tolerance , :) = [];  % keep only those with similar pressure in both methods
%{
for ind = 1 : length(a(:,1))
    if not(ismember(junA3,a(ind,:),'rows')) 
        if abs(2*(a(ind,6)-a(ind,7))/(a(ind,6)+a(ind,7)))<tolerance
            A = [A; a(ind,:)]; % keep only rows that is not found in the datafile and pressure is similar in both methods
        end
    end
end
newjunA =   [junA3; A];
newrows =   size(A,1);%length(A(:,1)); 
numofrows = size(newjunA,1);%length(newjunA(:,1));
%%
%}

%%
%a = newjunA;
%a = junA;
%a = A;
arm = a(:,8);
r = a(:,5);
Ve = a(:,3)-3*pi*arm.*r.^2;
P = mean([a(:,7),a(:,6)],2);
gamma = a(:,2);
kappa = 1;
lengthscale = sqrt(kappa./gamma);
Pscale = gamma.^(3/2)/kappa^(1/2);
fig1 = figure;
%NVe = Ve/lengthscale^3;
scatter(Ve./lengthscale.^3,P./Pscale)
title('pressure vs volume excess'); xlabel('norV_{e} [(\kappa/\gamma)^{3/2}]'); ylabel('norP [\gamma^{3/2}/\kappa^{1/2}]');
movegui(fig1,'northeast');
%hold on;
%ind = numofrows-newrows+1 : numofrows;
%scatter(Ve(ind)./lengthscale(ind).^3,P(ind)./Pscale(ind),'sr')
x = Ve./lengthscale.^3; y = P./Pscale;

fig2 = figure;
Fj = a(:,1)-3*pi*(kappa./r+2*gamma.*r).*arm;
%scatter(Ve,A(:,1));
scatter(Ve./lengthscale.^3,Fj);
title('junction - Energy vs volume excess'); xlabel('norV_{e} [(\kappa/\gamma)^{3/2}]'); ylabel('F_{exc-P} [\kappa]');
movegui(fig2,'northwest');
%hold on;
%scatter(Ve(ind)./lengthscale(ind).^3,Fj(ind),'sr');
%%
answer = questdlg('keep new data?');
if answer == "Yes"
    junA3 = a;
    save('junA3.mat','junA3');
end

close(fig1)
close(fig2)

%%
%answer = questdlg('continue uploading?');
%end


%{ 
% make an interpolating curve

[xi,yi] = getpts;
[xi,ix,~] = unique(xi);
yi = yi(ix);
VexJ = 5 : 1750;
PJ   = spline(xi,yi,VexJ); 
[VeN, ind] = sort(Ve./lengthscale.^3);
[VeN, ind2] = unique(VeN);
FeJ   = Fj(ind);
FeJ   = interp1(VeN,FeJ(ind2),VexJ);
%%%figure; scatter(VexJ,FeJ); title('junction - Energy vs volume excess'); xlabel('norV_{e} [(\kappa/\gamma)^{3/2}]'); ylabel('F_{exc-P} [\kappa]');
%hold on; plot(VexJ,PJ);
[~,imax] = max(PJ);
VexJL = VexJ(1:imax);
PJL = PJ(1:imax);
VexJR = VexJ(imax:end);
PJR = PJ(imax:end);
FeJL = FeJ(1:imax);
FeJR = FeJ(imax:end);
%}


%%
%%{
%junA3 = load('junA3.mat');    % load datafile
%junA3 = junA3.junA3;
%a = junA3;
arm = a(:,8);
r = a(:,5);
Ve = a(:,3)-3*pi*arm.*r.^2;
P = mean([a(:,7),a(:,6)],2);
gamma = a(:,2);
kappa = 1;
lengthscale = sqrt(kappa./gamma);
Pscale = gamma.^(3/2)/kappa^(1/2);
r = a(:,5);
fig3 = figure;
%NVe = Ve/lengthscale^3;
scatter(Ve./lengthscale.^3,P./Pscale)
title('pressure vs volume excess'); xlabel('norV_{e} [(\kappa/\gamma)^{3/2}]'); ylabel('norP [\gamma^{3/2}/\kappa^{1/2}]');
movegui(fig3,'northeast');
hold on;
ind = numofrows-newrows+1 : numofrows;
scatter(Ve(ind)./lengthscale(ind).^3,P(ind)./Pscale(ind),'sr')
x = Ve./lengthscale.^3; y = P./Pscale;

fig4 = figure;
Fj = a(:,1)-3*pi*(kappa./r+2*gamma.*r).*arm;
%scatter(Ve,A(:,1));
scatter(Ve./lengthscale.^3,Fj);
title('junction - Energy vs volume excess'); xlabel('norV_{e} [(\kappa/\gamma)^{3/2}]'); ylabel('F_{exc-P} [\kappa]');
movegui(fig4,'northwest');
hold on;
scatter(Ve(ind)./lengthscale(ind).^3,Fj(ind),'sr');
%%
%%}