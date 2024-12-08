% interpolate_jun_swelling_2.m

%% uses all results of pressure-volume relations (and other) and makes an interpolation.
%% takes the raw data from "junA3" interpolates and saving in --> "jun_inter".
%% version 2: takes the point P=0, calculated seperately and averaged over a few simulations, and costraining the interpolation through it (it more accurate than other simulations)
%% version 2: Finally, the point V=0;P=0, is very close to volume of natural (free volume) junctions -> this will be used

junA3 = load('junA3.mat');    % load datafile
junA3 = junA3.junA3;
sizeoffont = 14;

a = junA3;
arm = a(:,8);
r = a(:,5);
Ve = a(:,3)-3*pi*arm.*r.^2;
gamma = a(:,2);
kappa = 1;
lengthscale = sqrt(kappa./gamma);
Pscale = gamma.^(3/2)/kappa^(1/2);
P = mean([a(:,7),a(:,6)],2);
Fj = a(:,1)-3*pi*(kappa./r+2*gamma.*r).*arm;
Vex_P0 = 3.0;  % this is the exteranly calculated excess volume for P=0 but not very accurate
Vex_P0 = 0;  %temporarily
Fex_P0 = -15;  % this is the exteranly calculated excess energy for P=0
%%{
% make an interpolating curve
VexJ = 0.001 : 1450;
VexJ = Vex_P0 : 1000;
xi = Ve./lengthscale.^3;
yi = a(:,7)./Pscale;
%figure; scatter(xi,yi);

%%interpolate the energies for each volume %%
fE = fit([Vex_P0 ; xi],[Fex_P0 ; Fj],'smoothingspline','SmoothingParam',0.001,'Weights',[1000 ; ones(length(xi),1)]);   %% 
figure; graph1=plot(fE,'k--',[Vex_P0 ; xi],[Fex_P0 ; Fj],'.r'); xlabel('$$V_{exc} [\lambda^3]$$','Interpreter','latex'); ylabel('$$F [\kappa]$$','Interpreter','latex');
delfE = fit([Vex_P0 ; xi],[0 ; Fj-Fex_P0],'smoothingspline','SmoothingParam',0.001,'Weights',[1000 ; ones(length(xi),1)]);   %% 
figure; graph1=plot(delfE,'k--',[Vex_P0 ; xi],[0 ; Fj-Fex_P0],'.r'); xlabel('$$V_{exc} [\lambda^3]$$','Interpreter','latex'); ylabel('$$ \Delta F [\kappa]$$','Interpreter','latex');
FeJ   = feval(fE,[VexJ]);
graph1(2).LineWidth = 2;
%figure;plot(fP);
PJ   = feval(fP,VexJ);
%leg = legend([{'Simulation'} {'Fit'} {'Fit'}]);
legend('off');
xlim([0 1000]); xlim([0 635]); %xlim([0 max(VexJ)]);
set(gca,'FontSize',sizeoffont);


%% interpolate the pressure for each volume%%
% need to add weight to the first part (small volumes) for interpolation
weights = ones(length(xi),1); %weights(find(xi<200)) = 1;
decay_length = 100;
weights(xi>150) = exp((-xi(xi>150)+150)/decay_length);  % I made the weights exponentially decay as it gets to larger volumes than 150 because of numerical noise
fP = fit([Vex_P0 ; xi],[0 ; yi],'smoothingspline','SmoothingParam',0.01,'Weights',[1000 ; weights]);  % here I am using smoothing spline. the smoothing parameter is ? (how differ for the data vs how mow curviness (second derivative).
graph2 = plot(fP,'k--',[Vex_P0 ; xi],[0 ; yi],'.r'); xlabel('V_{exc} [(\kappa/\gamma)^{3/2}]'); ylabel('P [\gamma^{3/2}/\kappa^{1/2}]');
graph2(2).LineWidth = 2;
PJ   = feval(fP,VexJ);
leg
legend('off');
xlim([0 1000]); xlim([0 635]); xlim([0 max(VexJ)]);
ylim([0 0.56]);
set(gca,'FontSize',sizeoffont);
xlabel('$$V_{exc} [\lambda^3]$$','Interpreter','latex'); ylabel('$$p$$','Interpreter','latex');

% save the interpolation curves: VexJ,FeJ(VexJ),PJ(VexJ)
answer = questdlg('save new interpolation?');
if answer == "Yes"
    jun_inter = [VexJ ; FeJ' ; PJ'];                                   % excess volume ; excess energy ; pressure
    save('jun_inter.mat','jun_inter');
end

%figure; plot(jun_inter(1,:),jun_inter(2,:)); xlabel('normalized volume'); ylabel('excess energy');
%figure; plot(jun_inter(1,:),jun_inter(3,:)); xlabel('normalized volume'); ylabel('pressure');

%%%%%%%%%%%%%%%%%%%%%
%More master curves
%% radius.....
fR = fit(Ve./lengthscale.^3,r./lengthscale,'smoothingspline','SmoothingParam',0.0001); 
%figure; scatter(Ve./lengthscale.^3,r./lengthscale,'.r')
figure; graph2 = plot(fR,'k--',Ve./lengthscale.^3,r./lengthscale,'.r'); graph2(2).LineWidth = 2;
xlim([0 1000]); xlim([0 635]); %xlim([0 max(VexJ)]); ylim([0 0.56]);
set(gca,'FontSize',sizeoffont);
xlabel('$$V_{exc} [\lambda^3]$$','Interpreter','latex'); ylabel('$$R [\lambda]$$','Interpreter','latex');

