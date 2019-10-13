
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis of ICA time series according to Bilek 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Includes regression of artificial components, creating of the task
% design, MRA, component selection, cross-correlations between real and
% artificial pairs, randomisation and compatison between blocks
% Coop = cooperation; Comp = competition; CN = Concurrent; TB = Turn- based
% Inputs include results from gICA and timing file. Timing file used here
% is a struct encompassing {pair,session,role,condition}; e.g.
% timings.data{5,1,1,1} is a three-column array (onset-duration-weight - as
% used in FSL) listing all the successfull placements in Pair 5 in
% Concurrent block of the Builder in Cooperation




% ICA directory
icadir = ''; 
% ICA analysis prefix (was set in GIFT)
icaname = '';

% fMRI TR [sec]
TR = 2;

% dummy scans
dummy = 4;

% timings file
tf = '';

% artificial components
artcomps = [];
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of area that needs to be edited
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ICA info
icaf = spm_select('fplist',icadir,sprintf('%s.*_parameter_info.mat$',icaname));
load(icaf)

% load subjects' ICA timeseries
nS = sesInfo.numOfDataSets;
nT = sesInfo.numOfScans;
nC = sesInfo.numComp;
TS = zeros(nT,nC,nS);
compnames = '';
for ii = 1:nC
    compnames(ii,:) = sprintf('C%2.2d',ii);
end

for ii = 1:nS
    brfs{ii} = spm_select('fplist',icadir,sprintf('.*_br%d.mat$',ii));
    brtemp = load(deblank(brfs{ii}));
    TS(:,:,ii) = brtemp.compSet.tc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard artificial components
okcomps = setdiff(1:nC,artcomps);
nC = length(okcomps);
compnames = compnames(okcomps,:);
TS = TS(:,okcomps,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose components that show significant task effect
load(tf) 
for ss = 1:nS
    
    % get onsets and durations
    
    [pa na ex] = fileparts(sesInfo.userInput.files(ss).name(2,:));
    temp = strread(na,'%s',-1,'delimiter','_');
    
    pp(ss) = str2double(temp{1}(5:end));
    if ~isempty(strfind(temp{2},'Blue'))
        rr(ss) = 1;
    else
        rr(ss) = 2;
    end
    if ~isempty(strfind(temp{3},'CN'))
        ses(ss) = 1;
    else
        ses(ss) = 2;
    end
    
    ons1 = timings_HyPaGame.data{pp(ss),ses(ss),rr(ss),1}(:,1);
    ons2 = timings_HyPaGame.data{pp(ss),ses(ss),rr(ss),2}(:,1);
    ons3 = timings_HyPaGame.data{pp(ss),ses(ss),rr(ss),3}(:,1);
    durs1 = timings_HyPaGame.data{pp(ss),ses(ss),rr(ss),1}(:,2);
    durs2 = timings_HyPaGame.data{pp(ss),ses(ss),rr(ss),2}(:,2);
    durs3 = timings_HyPaGame.data{pp(ss),ses(ss),rr(ss),3}(:,2);
    
    % construct dependent variables (task design)
    td = [];
    td.xY.RT = TR;
    td.nscan = sesInfo.userInput.diffTimePoints(ss) + dummy;
    td.xBF.name = 'hrf';
    td.xBF.UNITS = 'secs';
    td.xBF.Volterra = 1;
    td.Sess.U(1).name = {'ev1'};
    td.Sess.U(1).ons = ons1;
    td.Sess.U(1).dur = durs1;
    td.Sess.U(1).P.name = 'none';
    td.Sess.C.C = [];
    td.Sess.C.name = cell(0);
    td.Sess.U(2).name = {'ev2'};
    td.Sess.U(2).ons = ons2;
    td.Sess.U(2).dur = durs2;
    td.Sess.U(2).P.name = 'none';
    td.Sess.U(3).name = {'ev3'};
    td.Sess.U(3).ons = ons3;
    td.Sess.U(3).dur = durs3;
    td.Sess.U(3).P.name = 'none';
    td = spm_fMRI_design(td,0);
    
    % compute beta values
    for cc = 1:3
        
        % model of response to condition cc
        Y = td.xX.X(dummy+1:end,cc);
        ys(:,ss) = Y(:,1);
        % design matrix composed by components' time series
        X = squeeze(TS(:,:,ss));
        % remove mean and add constant
        X = [X-repmat(mean(X),size(X,1),1) ones(size(X,1),1)];
        % get beta parameters
        b(ss,cc,:) = (X'*X)\(X'*Y);
        
    end   
end

% T-test on coop versus comp contrast 
con1 = squeeze(b(:,1,1:end-1)) - squeeze(b(:,2,1:end-1));
CN_con1 = con1(ses==1 ,:);
[CNh1,CNp1,CNci1,CNst1] = ttest(CN_con1);
TB_con1 = con1(ses==2 ,:);
[TBh1,TBp1,TBci1,TBst1] = ttest(TB_con1);

%T-test on coop vs ctrl
con2 = squeeze(b(:,1,1:end-1)) - squeeze(b(:,3,1:end-1));
CN_con2 = con2(ses==1 ,:);
[CNh2,CNp2,CNci2,CNst2] = ttest(CN_con2);
TB_con2 = con2(ses==2 ,:);
[TBh2,TBp2,TBci2,TBst2] = ttest(TB_con2);

%T-test on comp vs ctrl
con3 = squeeze(b(:,2,1:end-1)) - squeeze(b(:,3,1:end-1));
CN_con3 = con3(ses==1 ,:);
[CNh3,CNp3,CNci3,CNst3] = ttest(CN_con3);
TB_con3 = con3(ses==2 ,:);
[TBh3,TBp3,TBci3,TBst3] = ttest(TB_con3);

%T-test on comp vs coop
con4 = squeeze(b(:,2,1:end-1)) - squeeze(b(:,1,1:end-1));
CN_con4 = con4(ses==1 ,:);
[CNh4,CNp4,CNci4,CNst4] = ttest(CN_con4);
TB_con4 = con4(ses==2 ,:);
[TBh4,TBp4,TBci4,TBst4] = ttest(TB_con4);

mCN_bev1 = mean(squeeze(b(ses==1,1,1:end-1)));
mTB_bev1 = mean(squeeze(b(ses==2,1,1:end-1)));
mCN_bev2 = mean(squeeze(b(ses==1,2,1:end-1)));
mTB_bev2 = mean(squeeze(b(ses==2,2,1:end-1)));
mCN_bev3 = mean(squeeze(b(ses==1,3,1:end-1)));
mTB_bev3 = mean(squeeze(b(ses==2,3,1:end-1)));
mCN_bev = mean(squeeze(mean(b(ses==1,:,1:end-1))));
mTB_bev = mean(squeeze(mean(b(ses==2,:,1:end-1))));

% choose components of interest
%Coop > Comp
pthr = 0.05/nC; % Bonferoni p<0.05
wh1 = CNp1<pthr & CNst1.tstat>0 & TBp1<pthr & TBst1.tstat>0 & mCN_bev1>0 & mTB_bev1>0;
nC1 = sum(wh1);
compnames1 = compnames(wh1,:);

%COOP
pthr = 0.05/nC; % Bonferoni p<0.05
wh2 = CNp2<pthr & CNst2.tstat>0 & mCN_bev1>0 & TBp2<pthr & TBst2.tstat>0 &  mTB_bev1>0;
nC_Coop = sum(wh2);
compnames_Coop = compnames(wh2,:);

%COMP
pthr = 0.05/nC; % Bonferoni p<0.05
wh3 = CNp3<pthr & CNst3.tstat>0 &  mCN_bev2>0 & TBp3<pthr & TBst3.tstat>0 & mTB_bev2>0; 
nC_Comp = sum(wh3);
compnamesComp = compnames(wh3,:);

%Comp > Coop
pthr = 0.05/nC; % Bonferoni p<0.05
wh4 = CNp4<pthr & CNst4.tstat>0 & TBp4<pthr & TBst4.tstat>0 & mCN_bev2>0 & mTB_bev2>0;
nC4 = sum(wh4);
compnames4 = compnames(wh4,:);

wh = wh2 | wh3 | wh4 | wh1;
wh_x= wh2 | wh3
nC = sum(wh);
compnamesx=compnames(wh,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get real-pairs correlations
%
upp = unique(pp);
nP = length(upp);
Ctrue = zeros(nC,nC,nP,2);

for ii = 1:nP
    for ss = 1:2
        tsA = squeeze(TS(:,wh,pp==upp(ii) & ses==ss & rr==1));
        tsB = squeeze(TS(:,wh,pp==upp(ii) & ses==ss & rr==2));
        cc = corrcoef([tsA tsB]);
        Ctrue(:,:,ii,ss) = cc(1:nC,nC+1:end);
    end
end

zCtrue = .5*log((1+Ctrue)./(1-Ctrue))*sqrt(nT-3);
mzCtrue = squeeze(median(reshape(zCtrue,nC,nC,[]),3));
mCtrue = squeeze(median(reshape(Ctrue,nC,nC,[]),3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get non-pairs correlations distribution

nNP = 2*nP;
nPerm = 10000;
Cnontrue = zeros(nC,nC,nNP,nPerm);
allperm = [reshape(repmat(1:nP,nP,1),[],1) repmat((1:nP)',nP,1)];
allperm = allperm((allperm(:,1)-allperm(:,2))~=0,:);
allperm = [[allperm; allperm] reshape(repmat(1:2,size(allperm,1),1),[],1)];
nallperm = size(allperm,1);

fprintf('\n');
nch = fprintf('%2.2f%% permutations done...',0);
for ii = 1:nPerm
    
    curpermord = sortrows([rand(nallperm,1) (1:nallperm)'],1);
    
    for jj = 1:nNP
        perms(:,ii,jj) = allperm(curpermord(jj,2),:)';
        tsA = squeeze(TS(:,wh,pp==upp(perms(1,ii,jj)) & ses==perms(3,ii,jj) & rr==1));
        tsB = squeeze(TS(:,wh,pp==upp(perms(2,ii,jj)) & ses==perms(3,ii,jj) & rr==2));
        cc = corrcoef([tsA tsB]);
        Cnontrue(:,:,jj,ii) = cc(1:nC,nC+1:end);
    end
    fprintf(repmat(char(8),1,nch));
    nch = fprintf('%2.2f%% permutations done...',ii/nPerm*100);
end

zCnontrue = .5*log((1+Cnontrue)./(1-Cnontrue))*sqrt(nT-3);
mzCnontrue = squeeze(median(zCnontrue,3));
mCnontrue = squeeze(median(Cnontrue,3));

%%%%%%%%%%%%%%%%%%
% estimate p-value
for  ii = 1:nC
    for jj = 1:nC
        if median(squeeze(mzCnontrue(ii,jj,:))) < 0
            d = -squeeze(mzCtrue(ii,jj)-mzCnontrue(ii,jj,:));
        else
            d = squeeze(mzCtrue(ii,jj)-mzCnontrue(ii,jj,:));
        end
        p(ii,jj) = sum(d<0)/nPerm;
    end
end

%%%%%%%%%%
% Bonferoni 0.05
uc = 0.05/numel(p);
[aa bb] = ind2sub([nC nC],find(p<=uc));

% show res
figure(10), clf
for ii = 1:nC
    for jj = 1:nC
        subplot(nC,nC,sub2ind(size(p),jj,ii)), cla
        hist(squeeze(mCnontrue(ii,jj,:)),30)
        title(sprintf('Blue%s vs Yellow%s',compnamesx(ii,:),compnamesx(jj,:)))
        set(get(gca,'Title'),'FontSize',6)
        yL = get(gca,'Ylim');
        line(mCtrue(ii,jj)*[1 1],[0 .8*yL(2)],'color','r')
    end
end
        
for ii = 1:length(aa) 
    subplot(nC,nC,sub2ind(size(p),bb(ii),aa(ii)))
    %set(get(gca,'Title'),'String',sprintf('Prop%s vs Resp%s',compnames(aa(ii),:),compnames(bb(ii),:)))
    set(get(gca,'Title'),'FontSize',8)
    set(get(gca,'Title'),'Color',[0.8 0 0])
end
    
% % uncorrected 0.05
% [aa bb] = ind2sub([nC nC],find(p<=0.05));
% %[aa bb]
% 
% % show res
% figure(11), clf
% for ii = 1:size(aa,1)
%     subplot(size(aa,1),1,ii)
%     hist(squeeze(mCnontrue(aa(ii),bb(ii),:)),floor(nPerm/100))
%     title(sprintf('uncorr 0.05   C%2.2d  vs C%2.2d',aa(ii),bb(ii)))
%     yL = get(gca,'Ylim');
%     line(mCtrue(aa(ii),bb(ii))*[1 1],[0 .8*yL(2)],'color','r')
% end    
%}

% testing CN vs TB for true pairs
for ii = 1:nC
    for jj = 1:nC
        
        [p1(ii,jj) h(ii,jj) stats(ii,jj)]= signrank(squeeze(Ctrue(ii,jj,:,1)),squeeze(Ctrue(ii,jj,:,2)));
        
    end
end


uc1 = 0.05/numel(p1);
[aa1 bb1] = ind2sub([nC nC],find(p1<=0.05));

for ii = 1:length(aa1)
    figure(10+ii), clf, hold on
    predata = squeeze(Ctrue(aa1(ii),bb1(ii),:,1));
    postdata = squeeze(Ctrue(aa1(ii),bb1(ii),:,2));
    plot(ones(length(predata),1),predata,'or')
    plot(2*ones(length(predata),1),postdata,'or')
    for jj = 1:length(predata)
        line([1 2],[predata(jj) postdata(jj)])
    end
    set(gca,'xlim',[0.75 2.25])
    m1 = plot(0.9,mean(predata),'o');
    set(m1,'MarkerSize',15) 
    set(m1,'Color',[1 0 0]) 
    set(m1,'LineWidth',2) 
    m2 = plot(2.1,mean(postdata),'o');
    set(m2,'MarkerSize',15) 
    set(m2,'Color',[1 0 0]) 
    set(m2,'LineWidth',2) 
    L = line([0.9 2.1],[mean(predata) mean(postdata)]);
    set(L,'LineWidth',3)
    set(L,'LineStyle','--')
    title(sprintf('CN vs TB  Blue%s vs Yellow%s',compnamesx(aa1(ii),:),compnamesx(bb1(ii),:)))
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'CN' 'TB'})
    yl = get(gca,'ylabel');
    set(yl,'String','Correlation Coefficient')
    set(yl,'FontSize',15)
    
end
%}
