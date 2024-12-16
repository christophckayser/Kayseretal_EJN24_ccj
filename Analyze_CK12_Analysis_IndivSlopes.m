clear;

EXP =2
% AV-A experiment with two auditory reliabilities and CCJ judgments. Study CK12 Exp1/2
% -> Create figures showing bias and CCJ vs DVA 
%    Slopes for individual participants.
%    Analyze slopes for effects of Rel and CCJ
%    Analyze CCJ itself

% Data
% audrel 1 is high, 2 is low
%  CCJ is 1 for no, 2 yes
% dataformat in AllData{subject}
% [VEbias, VAEbias, Vpos, Apos, DVA, Rel, APosA, respAV, respA, CCJ, Trial  sub]

if EXP==1
  load('F:/CKDATA/Projects/projects/Multisensory_decoding/CK12_ccj/CK12_Exp1_Alldata.mat','AllData')
else
  load('F:/CKDATA/Projects/projects/Multisensory_decoding/CK12_ccj/CK12_EXP2_Alldata.mat','AllData')
end

% ----------------------------------------------------------------------------------
% collect data and compute individual slopes
% ----------------------------------------------------------------------------------
dva = [-33:11:33];
apos = [-22:11:22];
ndva = length(dva);
nsub =length(AllData);
for s=1:nsub
  data = AllData{s};
  %------------------------------------------------------------
  % unisensory SD  % for each position compute SD and average across positions
  %------------------------------------------------------------
  tmpx=[]; tmpy=[];
  for a=1:length(apos)
    % a-high
    j = find( (data(:,6)==1).*(data(:,4)==apos(a)));
    tmpx(a,1) = std(data(j,8));
    tmpy(a,1) = mean(data(j,8));
    % a-low
    j = find( (data(:,6)==2).*(data(:,4)==apos(a)));
    tmpx(a,2) = std(data(j,8));
    tmpy(a,2) = mean(data(j,8));
  end
  UniSD(s,:) = mean(tmpx);

  %------------------------------------------------------------
  % compute mean CCJ  for each dva and reliability
  %------------------------------------------------------------
  for m=1:ndva
    for rel=1:2
      j = find( (data(:,5)==dva(m)).*( data(:,6)==rel));
      CCJ(s,rel,m) = mean(data(j,10),'omitnan')-1;
    end
  end

  % -----------------------------------------------------
  % biases per CCJ
  % -----------------------------------------------------
  % VE per DVA and reliability and CCJ
  for m=1:ndva
    for rel=1:2
      for c=1:2
        j = find( (data(:,5)==dva(m)).*( data(:,6)==rel).*( data(:,10)==c));
        VE_bias(s,m,rel,c) = mean(data(j,1),'omitnan');
      end
    end
  end
  % VE slope vs dva
  for rel=1:2
    for c=1:2
      j = find(( data(:,6)==rel).*( data(:,10)==c));
      m = data(j,5);
      m(:,2) = 1;
      slope = regress(data(j,1),m);
      VE_slope(s,rel,c) = slope(1); % slope vs DVA only
    end
  end

  % VAE per DVA and reliability and CCJ
  for m=1:ndva
    for rel=1:2
      for c=1:2
        j = find( (data(:,5)==dva(m)).*( data(:,6)==rel).*( data(:,10)==c));
        VAE_bias(s,m,rel,c) = mean(data(j,2),'omitnan');
      end
    end
  end
  % VAE slope
  for rel=1:2
    for c=1:2
      j = find(( data(:,6)==rel).*( data(:,10)==c));
      m = data(j,5);
      m(:,2) = 1;
      slope = regress(data(j,2),m);
      VAE_slope(s,rel,c) = slope(1);
    end
  end
end

% ----------------------------------------------------------------------------------
%  analyze unisensory SDs
% ----------------------------------------------------------------------------------
for t=1:2
  x = UniSD(:,t);
  % percentile bootstrap CI
  [Lo,Up,p]=ck_stat_bootCI_percentile(x,'median',0.01,8000,0);
  CIs_SD(t,:) =[median(x) Lo Up];
end
fprintf('-------------------------------------------------------- \n')
fprintf('Unisensory SD (Median and CI): \n');
for t=1:2
  fprintf('%1.3f [%1.3f, %1.3f]  \n',CIs_SD(t,1),CIs_SD(t,2),CIs_SD(t,3));
end

% paired t-test for effect of A reliability on unisensory response SD
% THIS REQUIRES THE BF TOOLBOX FOR MATLAB 
% % See https://klabhub.github.io/bayesFactor/ 
[bf10,pValue,CI,stats] =  bf.ttest(UniSD(:,1),UniSD(:,2));
fprintf('Effect of A rel: t=%1.2f p=%1.3f   BF=%1.2f \n',stats.tstat, pValue,bf10)

% ----------------------------------------------------------------------------------
% test indiv. conditions vs zero
% ----------------------------------------------------------------------------------
for rel=1:2
  for c=1:2
    [h,p,~,stat] = ttest(VE_slope(:,rel,c),0);
    VE_test(rel,c,:) = [p,stat.tstat];
    [h,p,~,stat] = ttest(VAE_slope(:,rel,c),0);
    VAE_test(rel,c,:) = [p,stat.tstat];
  end
end



% ----------------------------------------------------------------------------------
% figure
% ----------------------------------------------------------------------------------

Colorvector(1,:) = [250,150,40]/250; %
Colorvector(2,:) = [180,40,40]/250; %
Colorvector(3,:) = [40,150,250]/250; %
Colorvector(4,:) = [40,40,180]/250; %
Rel_L = {'A+','A-'};
CC_L = {'CCJ-','CCJ+'};
colind = [1,2;3,4];


figure(EXP);clf;

%-------------------------------------------
% VE for each reliability and CCJ
subplot(2,3,2); hold on
line([-40 40],[0 0],'LineWidth',1,'color',ones(1,3)*0.7)
for d=1:2 % reliability
  for c=1:2
    co = colind(c,d);
    x = VE_bias(:,:,d,c);
    h = errorbar(dva,nanmean(x),sem(x));
    set(h, 'Color', Colorvector(co,:),'LineWidth',2);
    h = text(15,-25+(co-1)*4,sprintf('%s%s',Rel_L{d},CC_L{c}),'FontSize',10,'Fontweight','bold'); set(h, 'Color', Colorvector(co,:));
  end
end
xlabel('discrepancy (\DeltaVA)[\circ]');
axis([-40 40 -30 30]);
set(gca,'YTick',[-30:10:30])
ylabel('judgement error [deg]')
title('ventriloquism effect','FontSize',13,'FontWeight','Bold');
text(-58,35,'B','FontSize',20,'FontWeight','bold');


% VE slopes 
subplot(2,3,5); hold on
a = reshape(VE_slope(:,:,:),[size(VE_slope,1),4,1]);
ckmeanplotcompact(a,4,1,0,Colorvector);
axis([0.5 4.75 -0.5 1.5]);
set(gca,'XTick',[]);
set(gca,'YTick',[-.5:0.5:1.5])
line([-0.5 4.75],[0 0],'LineWidth',1,'color',ones(1,3)*0.7)
ylabel('slope')
for d=1:2 % reliability
  for c=1:2
    co = colind(c,d);
    h = text(0.7+(co-1)*1,1.4,sprintf('%s%s',Rel_L{d},CC_L{c}),'FontSize',10,'Fontweight','bold'); set(h, 'Color', Colorvector(co,:));
    if VE_test(d,c,1)>0.0001
      h = text(0.7+(co-1)*1,-0.3,sprintf('p=%1.4f',VE_test(d,c,1)),'FontSize',8);;
    else
      h = text(0.7+(co-1)*1,-0.3,sprintf('p<0.0001'),'FontSize',8);;
    end
  end
end


%-------------------------------------------
% VAE against dva for each reliability 
subplot(2,3,3); hold on
line([-40 40],[0 0],'LineWidth',1,'color',ones(1,3)*0.7)

for d=1:2 % reliability
  for c=1:2
    co = colind(c,d);
    x = VAE_bias(:,:,d,c);
    h = errorbar(dva,nanmean(x),sem(x));
    set(h, 'Color', Colorvector(co,:),'LineWidth',2);
    h = text(15,-12.5+(co-1)*2,sprintf('%s%s',Rel_L{d},CC_L{c}),'FontSize',10,'Fontweight','bold'); set(h, 'Color', Colorvector(co,:));
  end
end
xlabel('discrepancy (\DeltaVA)[\circ]');
axis([-40 40 -15 15]);
set(gca,'YTick',[-15:5:15]);
ylabel('judgement error [deg]')
title('aftereffect','FontSize',13,'FontWeight','Bold');
text(-58,18,'C','FontSize',20,'FontWeight','bold');
set(gca,'Position',[0.72    0.5838    0.2134    0.3412])


% VAE slopes 
subplot(2,3,6); hold on
a = reshape(VAE_slope,[size(VAE_slope,1),4,1]);
ckmeanplotcompact(a,4,1,0,Colorvector);
axis([0.5 4.75 -0.1 0.4]);
line([-0.5 4.75],[0 0],'LineWidth',1,'color',ones(1,3)*0.7)
set(gca,'XTick',[]);
set(gca,'YTick',[-0.1:0.1:0.4])
ylabel('slope')
set(gca,'Position',[0.72    0.1100    0.2134    0.3412])

for d=1:2 % reliability
  for c=1:2
    co = colind(c,d);
    h = text(0.7+(co-1)*1,0.37,sprintf('%s%s',Rel_L{d},CC_L{c}),'FontSize',10,'Fontweight','bold'); set(h, 'Color', Colorvector(co,:));
    if VAE_test(d,c,1)>0.0001
      h = text(0.7+(co-1)*1,-0.06,sprintf('p=%1.4f',VAE_test(d,c,1)),'FontSize',8);;
    else
      h = text(0.7+(co-1)*1,-0.06,sprintf('p<0.0001'),'FontSize',8);;
    end
  end
end

subplot(2,3,1); hold on
for c=1:2
  x = sq(CCJ(:,c,:));
  h = errorbar(dva,nanmean(x),sem(x));
  set(h, 'Color', Colorvector(c,:),'LineWidth',2);
  h = text(-30+(c-1)*10,0.9,Rel_L{c},'FontSize',12,'Fontweight','bold'); set(h, 'Color', Colorvector(c,:));
end
xlabel('discrepancy (\DeltaVA)[\circ]');
ylabel('ccj')
axis([-40 40 0 1]);
set(gca,'YTick',[0:0.2:1])
title('common cause judgement','FontSize',13,'FontWeight','Bold');
text(-56,1.2,'A','FontSize',20,'FontWeight','bold');
set(gca,'Position',[0.1000    0.3838    0.2134    0.3412])


% cosmetics for the figures
ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',12);
set(gcf,'Position',[      125         537        1306         735]);
FigureDir = 'F:\CKDATA\Projects\projects\Multisensory_decoding\CK12_ccj\figures';

snamef = sprintf('%s/CK12_Result_exp%d.jpg',  FigureDir,EXP);
print('-djpeg','-r300',snamef);



% ----------------------------------------------------------------------------------
%% lmes
% ----------------------------------------------------------------------------------

% collect data into table
X=[];
for s=1:size(VE_slope,1)  
  for rel=1:2
    for c=1:2
      
      X = cat(1,X,[VE_slope(s,rel,c),VAE_slope(s,rel,c),rel,c,s]);
    end
  end
end

Table = array2table(X,'VariableNames',{'VE','VAE','Rel','ccj','S'});

%-------------------------------------------------------------------------------
% Ve slopes effect of Reliability

Model=[]; bic=[]; Evar=[];
Model{1} = 'VE ~ Rel*ccj + (1|S) + (ccj|S) + (Rel|S)';
Model{2} = 'VE ~ ccj + (1|S) + (ccj|S)';
for m=1:length(Model)
  Result{m} = fitglme(Table,Model{m}, 'Distribution','normal','Link','identity','FitMethod','MPL','DummyVarCoding','reference','InitPLIterations',20,'PLIterations',500);
  bic(m) = Result{m}.ModelCriterion.BIC;
end
anova(Result{1})
BF = exp( (bic(2)-bic(1))/2); if BF<1, BF = -1./BF; end;
bic = bic-bic(1); % positive bic2 is evidence FOR Reliabiliyt
fprintf('VE deltaBIC for Rel: %1.f  BF %5.1f \n',bic(2),BF)


% code for table in paper
if 1==1
  fprintf('Table VE \n');
  for f=1:4
    fprintf('%1.3f + %1.3f',Result{1}.Coefficients.Estimate(f),Result{1}.Coefficients.SE(f));
    fprintf('\t');
    fprintf('%1.3f',Result{1}.Coefficients.tStat(f));
    fprintf('\t');
    if Result{1}.Coefficients.pValue(f)>0.0001
      fprintf('%1.4f',Result{1}.Coefficients.pValue(f));
    else
      fprintf('<0.0001');
    end
    fprintf('\n')
  end
end


%-------------------------------------------------------------------------------
% Ve slopes effect of CCJ
Model=[]; bic=[]; Evar=[];
Model{1} = 'VE ~ Rel*ccj + (1|S) + (ccj|S) + (Rel|S)';
Model{2} = 'VE ~ Rel + (1|S) + (Rel|S)';
for m=1:length(Model)
  Result{m} = fitglme(Table,Model{m}, 'Distribution','normal','Link','identity','FitMethod','MPL','DummyVarCoding','reference','InitPLIterations',20,'PLIterations',500);
  bic(m) = Result{m}.ModelCriterion.BIC;
end
BF = exp( (bic(2)-bic(1))/2); if BF<1, BF = -1./BF; end;
bic = bic-bic(1); % positive bic2 is evidence FOR Reliabiliyt
fprintf('VE deltaBIC for CCJ: %1.f  BF %5.1f \n',bic(2),BF)


%-------------------------------------------------------------------------------
% Vae slopes effect of Reliability

Model=[]; bic=[]; Evar=[];
Model{1} = 'VAE ~ Rel*ccj + (1|S) + (ccj|S) + (Rel|S)';
Model{2} = 'VAE ~ ccj + (1|S) + (ccj|S)';
for m=1:length(Model)
  Result{m} = fitglme(Table,Model{m}, 'Distribution','normal','Link','identity','FitMethod','MPL','DummyVarCoding','reference','InitPLIterations',20,'PLIterations',500);
  bic(m) = Result{m}.ModelCriterion.BIC;
end
anova(Result{1})
BF = exp( (bic(2)-bic(1))/2); if BF<1, BF = -1./BF; end;
bic = bic-bic(1); % positive bic2 is evidence FOR Reliabiliyt
fprintf('VAE deltaBIC for Rel: %1.f  BF %5.1f \n',bic(2),BF)

% code for table in paper
if 1==1
  fprintf('Table VAE \n');
  for f=1:4
    fprintf('%1.3f + %1.3f',Result{1}.Coefficients.Estimate(f),Result{1}.Coefficients.SE(f));
    fprintf('\t');
    fprintf('%1.3f',Result{1}.Coefficients.tStat(f));
    fprintf('\t');
    if Result{1}.Coefficients.pValue(f)>0.0001
      fprintf('%1.4f',Result{1}.Coefficients.pValue(f));
    else
      fprintf('<0.0001');
    end
    fprintf('\n')
  end
end

%-------------------------------------------------------------------------------
% VAe slopes effect of CCJ
Model=[]; bic=[]; Evar=[];
Model{1} = 'VAE ~ Rel*ccj + (1|S) + (ccj|S) + (Rel|S)';
Model{2} = 'VAE ~ Rel + (1|S) + (Rel|S)';
for m=1:length(Model)
  Result{m} = fitglme(Table,Model{m}, 'Distribution','normal','Link','identity','FitMethod','MPL','DummyVarCoding','reference','InitPLIterations',20,'PLIterations',500);
  bic(m) = Result{m}.ModelCriterion.BIC;
  Evar(m) = Result{m}.Rsquared.Adjusted;
end
BF = exp( (bic(2)-bic(1))/2); if BF<1, BF = -1./BF; end;
bic = bic-bic(1); % positive bic2 is evidence FOR Reliabiliyt
fprintf('VAE deltaBIC for CCJ: %1.f  BF %5.1f \n',bic(2),BF)






%% -------------------------------------------------------------------------
% statistics on CCJ:
% we model the ccj against the magnitude of the discrepancy
X=[];
for s=1:length(AllData)
  X = cat(1,X,AllData{s});
end
j = ~isnan(X(:,4));
X(:,13) = sign(X(:,5)).* sqrt(abs(X(:,5)));

Table = array2table(X(j,[5,6,10,12]),'VariableNames',{'DVA','Rel','ccj','S'});
Table.DVA = abs(Table.DVA);
Table.ccj = Table.ccj-1;

% effect of reliabiliy
Model=[]; bic=[];
Model{1} = 'ccj ~ DVA + Rel:DVA + Rel + (1|S)';
Model{2} = 'ccj ~ DVA +  (1|S) ';
for m=1:length(Model)
  Result{m} = fitglme(Table,Model{m}, 'Distribution','binomial','Link','logit','FitMethod','ApproximateLaplace','DummyVarCoding','reference','InitPLIterations',20,'PLIterations',500);
  bic(m) = Result{m}.ModelCriterion.BIC;
end
BF = exp( (bic(2)-bic(1))/2);
if BF<1, BF = -1./BF; end;

bic = bic-bic(1); % positive bic2 is evidence FOR 
fprintf('-------------------------------------------------------- \n')
Result{1}.Coefficients
fprintf('CCJ deltaBIC for Rel: %1.f  BF %5.1f \n',bic(2),BF)


% effect of DVA
Model=[]; bic=[];
Model{1} = 'ccj ~ DVA + Rel:DVA + Rel + (1|S) ';
Model{2} = 'ccj ~ Rel + (1|S) + (Rel | S)';
for m=1:length(Model)
  Result{m} = fitglme(Table,Model{m}, 'Distribution','binomial','Link','logit','FitMethod','ApproximateLaplace','DummyVarCoding','reference','InitPLIterations',20,'PLIterations',500);
  bic(m) = Result{m}.ModelCriterion.BIC;
end
BF = exp( (bic(2)-bic(1))/2);
if BF<1, BF = -1./BF; end;

bic = bic-bic(1); % positive bic2 is evidence FOR 
fprintf('CCJ deltaBIC for DVA: %1.f  BF %5.1f \n',bic(2),BF)

% code for table in paper
if 1==1
  fprintf('Table CCJ \n');
  for f=1:4
    fprintf('%1.3f + %1.3f',Result{1}.Coefficients.Estimate(f),Result{1}.Coefficients.SE(f));
    fprintf('\t');
    fprintf('%1.3f',Result{1}.Coefficients.tStat(f));
    fprintf('\t');
    if Result{1}.Coefficients.pValue(f)>0.0001
      fprintf('%1.4f',Result{1}.Coefficients.pValue(f));
    else
      fprintf('<0.0001');
    end
    fprintf('\n')
  end
end

