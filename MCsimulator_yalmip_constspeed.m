yalmip('clear');
clear all; clc; close all;
%Add gurobi solver path (require the second order cone constraint package)
global m g dt N_sim_end
rng(100);
m = 1; g = 9.81; dt = 0.05;
Ac = zeros(8,8); Ac(1,2) = 1; Ac(3,4) = 1;  Ac(5,6) = 1;  Ac(7,8) = 1; 
Bc = zeros(8,4); Bc(2,3) = -g; Bc(4,2) = g; Bc(6,1) = 1/m; Bc(8,4) = 1;
Cd = zeros(4,8); Cd(1,1) = 1; Cd(2,3) = 1; Cd(3,5) = 1; Cd(4,7) = 1;
N_sim_end = 200;
%% Obs define
ttraj = 0:dt:N_sim_end*dt;
%%
noiselevel = 0.25;
Boostap_yn = 1; Nboostraps = 40; Nconvergemax = 3;
rmaxboot = 8;
traininingadd = 5;
NSSA_training = 100;
L = floor(NSSA_training/4)-1; 
rmax = min(floor(L/2),10);
isBootstrap = 0;
Ueq = [m*g;0;0;0];
ucomd = zeros(4,1);
rbuff = 0.35;
nhorizon = 10; nstates = 8; ncontrols = 4;
% x_d(t), y_d(t), z_d(t), psi_d(t)
temptemp = expm([Ac,Bc;zeros(ncontrols,ncontrols+nstates)]*dt);
Ad = temptemp(1:nstates,1:nstates);
Bd = temptemp(1:nstates,nstates+1:nstates+ncontrols);
alpha = 50;
alpha_backup = 1000;
Qcost = diag([0.1,0.1,0.1,1]);
Rcost = diag([0.1,0.1,0.1,0.1]);  beta = 0.25;
Cobs = zeros(3,8); Cobs(1,1) = 1; Cobs(2,3) = 1; Cobs(3,5) = 1;
constantM = [Cobs,eye(3);-Cobs,eye(3)];
xyzINDEX = [1,3,5];
%% Controller with Straps
tic
    uvar_straps = sdpvar(repmat(ncontrols,1,nhorizon),repmat(1,1,nhorizon));
    xvar_straps = sdpvar(repmat(nstates,1,nhorizon+1),repmat(1,1,nhorizon+1));
    svar_straps = sdpvar(repmat(3,1,nhorizon+1),repmat(1,1,nhorizon+1));
    rvar_straps =  sdpvar(repmat(ncontrols,1,nhorizon+1),repmat(1,1,nhorizon+1));
    Aineq_straps = sdpvar(repmat(3+3,1,nhorizon+1),repmat(1,1,nhorizon+1));
    sigma_ak_straps = sdpvar(repmat(3,1,nhorizon+1),repmat(3,1,nhorizon+1));
    bineq_straps = sdpvar(repmat(3*2+1,1,nhorizon+1),repmat(1,1,nhorizon+1)); 
    pastu_straps = sdpvar(ncontrols,1);
    xbarvar_straps = sdpvar(repmat(nstates,1,nhorizon+1),repmat(1,1,nhorizon+1));
    trust_param_straps = sdpvar(1);
    objective_straps = 0; constraints_straps = [];
    for k = 1:nhorizon
        objective_straps = objective_straps + (Cd*xvar_straps{k}-rvar_straps{k})'*Qcost*(Cd*xvar_straps{k}-rvar_straps{k}) + uvar_straps{k}'*Rcost*uvar_straps{k};
        constraints_straps = [constraints_straps, xvar_straps{k+1} == Ad*xvar_straps{k}+Bd*uvar_straps{k}];
        constraints_straps = [constraints_straps, [-50;-pi;-pi;-pi] <= uvar_straps{k}<= [50;pi;pi;pi]];
    end
    for k = 2:nhorizon
        constraints_straps = [constraints_straps, [[[Aineq_straps{k}(1:3)'*Cobs,Aineq_straps{k}(4:end)'];sigma_ak_straps{k}*Cobs,-eye(3);-sigma_ak_straps{k}*Cobs,-eye(3)]*[xvar_straps{k};svar_straps{k}] <= bineq_straps{k}]];
    end  
    for k = 1:nhorizon-1
        constraints_straps = [constraints_straps, (xvar_straps{k} - xbarvar_straps{k})'*(xvar_straps{k} - xbarvar_straps{k}) <= trust_param_straps];
    end
    objective_straps = objective_straps + (Cd*xvar_straps{nhorizon+1}-rvar_straps{nhorizon+1})'*Qcost*(Cd*xvar_straps{nhorizon+1}-rvar_straps{nhorizon+1});
    parameters_in_straps = {xvar_straps{1},rvar_straps{:},Aineq_straps{:},sigma_ak_straps{:},bineq_straps{:},xbarvar_straps{:},trust_param_straps};
    solutions_out_straps = {[uvar_straps{:}], [xvar_straps{:}]};
   controller_straps = optimizer(constraints_straps, objective_straps,sdpsettings('solver','gurobi','gurobi.qcpdual',1, 'verbose',1,'debug',1),parameters_in_straps,solutions_out_straps);
%% Controller No Straps
    uvar_nostraps = sdpvar(repmat(ncontrols,1,nhorizon),repmat(1,1,nhorizon));
    xvar_nostraps = sdpvar(repmat(nstates,1,nhorizon+1),repmat(1,1,nhorizon+1));
    rvar_nostraps =  sdpvar(repmat(ncontrols,1,nhorizon+1),repmat(1,1,nhorizon+1));
    Aineq_nostraps = sdpvar(repmat(3,1,nhorizon+1),repmat(1,1,nhorizon+1));
    bineq_nostraps = sdpvar(repmat(1,1,nhorizon+1),repmat(1,1,nhorizon+1)); 
    pastu_nostraps = sdpvar(ncontrols,1);
    xbarvar_nostraps = sdpvar(repmat(nstates,1,nhorizon+1),repmat(1,1,nhorizon+1));
    trust_param_nostraps = sdpvar(1);
    objective_nostraps = 0; constraints_nostraps = [];
    for k = 1:nhorizon
        objective_nostraps = objective_nostraps + (Cd*xvar_nostraps{k}-rvar_nostraps{k})'*Qcost*(Cd*xvar_nostraps{k}-rvar_nostraps{k}) + uvar_nostraps{k}'*Rcost*uvar_nostraps{k};
        constraints_nostraps = [constraints_nostraps, xvar_nostraps{k+1} == Ad*xvar_nostraps{k}+Bd*uvar_nostraps{k}];
        constraints_nostraps = [constraints_nostraps, [-50;-pi;-pi;-pi] <= uvar_straps{k}<= [50;pi;pi;pi]];
    end
    for k = 2:nhorizon
        constraints_nostraps = [constraints_nostraps, [Aineq_nostraps{k}'*Cobs*xvar_nostraps{k} <= bineq_nostraps{k}]];
    end  
    for k = 1:nhorizon-1
        constraints_nostraps = [constraints_nostraps, (xvar_nostraps{k} - xbarvar_nostraps{k})'*(xvar_nostraps{k} - xbarvar_nostraps{k}) <= trust_param_nostraps];
    end
    objective_nostraps = objective_nostraps + (Cd*xvar_nostraps{nhorizon+1}-rvar_nostraps{nhorizon+1})'*Qcost*(Cd*xvar_nostraps{nhorizon+1}-rvar_nostraps{nhorizon+1});
    parameters_in_nostraps = {xvar_nostraps{1},rvar_nostraps{:},Aineq_nostraps{:},bineq_nostraps{:},xbarvar_nostraps{:},trust_param_nostraps};
    solutions_out_nostraps = {[uvar_nostraps{:}], [xvar_nostraps{:}]}; 
    %%
   controller_nostraps = optimizer(constraints_nostraps, objective_nostraps,sdpsettings('solver','gurobi','gurobi.qcpdual',1,'verbose',1,'debug',1),parameters_in_nostraps,solutions_out_nostraps);
   t_complie = toc
%% risk level
epsilon= 0.5;
v_epsilon = sqrt((1-epsilon)/epsilon);
%% MC simulator start here
nMonteCarlo = 1000;
iMonteCarlo = 1;
failedconverge = 0;
while iMonteCarlo <= nMonteCarlo
    crashed = 0;
    omittrun_badtrajectory =0;
    failedconverge = 0;
    icounter = 1; tode = 0;
    obsdist = [];
    r_ref = []; xtraj = []; q0 =[]; genflag = 0;
    while genflag == 0
        try
        [r_ref,xtraj,q0,indexcollosion] = fun_const_spd_MC();
        catch 
            genflag = 0;
        end
        genflag = 1;
        if isempty(r_ref)
            genflag = 0;
        end
        if isempty(xtraj)
            genflag = 0;
        end
        if isempty(q0)
            genflag = 0;
        end
    end
    noise_xtraj = addwhitenoise(xtraj, noiselevel,ttraj);
    noise_xtraj_List(iMonteCarlo,:,:) = noise_xtraj;
    xtraj_List(iMonteCarlo,:,:) = xtraj;
    indexcollosion_list(iMonteCarlo) = indexcollosion;
    isBootstrap = 0;
     Phistraps = []; 
     r_straps = [];
    looprelaxcounter = [];
    xu0 = [q0;ucomd];   
    uhist = []; qode_plot = q0'; 
    x_bar = [r_ref(1:nhorizon+1,1)';zeros(1,nhorizon+1);r_ref(1:nhorizon+1,2)';zeros(1,nhorizon+1);r_ref(1:nhorizon+1,3)';zeros(1,nhorizon+1);r_ref(1:nhorizon+1,4)';zeros(1,nhorizon+1)];
    Hankelx = [];
    Hankely = [];
    Hankelz = [];
    s_obs_traj = [];
    Urlistx = []; Urlisty = []; Urlistz = [];
    rlistx = []; rlisty = []; rlistz = [];
    Philistx = []; Philisty = []; Philistz = [];
    NSSA_training = 100;
        while icounter <= N_sim_end - nhorizon-1    
    Ex_obs_forecast = [];
    s_obs_measure = noise_xtraj(icounter,:)';
    s_obs_traj =  [s_obs_traj,s_obs_measure];
    inputsstrap = cell(1,6+5*nhorizon+1);
    inputs_nostrap = cell(1,5+4*nhorizon+1);
    tic
    Hankelx = hankel_add(Hankelx,s_obs_measure(1),L,icounter);
    Hankely = hankel_add(Hankely,s_obs_measure(2),L,icounter);
    Hankelz = hankel_add(Hankelz,s_obs_measure(3),L,icounter);

    looprelaxcounter(icounter) = 0;
    if icounter >= NSSA_training && isBootstrap == 0
        [Urlistxnew,rlistxnew,Philistxnew] = bootstrap_model(s_obs_traj(1,1:icounter),Hankelx,L,icounter,rmaxboot);
        [Urlistynew,rlistynew,Philistynew] = bootstrap_model(s_obs_traj(2,1:icounter),Hankely,L,icounter,rmaxboot);
        [Urlistznew,rlistznew,Philistznew] = bootstrap_model(s_obs_traj(3,1:icounter),Hankelz,L,icounter,rmaxboot);
        Urlistx = [Urlistx,Urlistxnew];
        rlistx = [rlistx,rlistxnew];
        Philistx = [Philistx;Philistxnew];
        Urlisty = [Urlisty,Urlistynew];
        rlisty = [rlisty,rlistynew];
        Philisty = [Philisty;Philistynew];
        Urlistz = [Urlistz,Urlistznew];
        rlistz = [rlistz,rlistznew];
        Philistz = [Philistz;Philistznew];
        NSSA_training = NSSA_training+traininingadd;
        if length(rlistx) >= Nboostraps
            isBootstrap = 1;
        end    
    end
    %%
    if isBootstrap == 0
        for istrap = 1:Nboostraps
            Ex_obs_forecast(istrap,:,:) = s_obs_measure*ones(1,nhorizon+1);
        end
    else
        for istrap = 1:Nboostraps
           Fstrapsfullx = [reconstructSSA_partial(Urlistx{istrap}*(Hankelx(:,end-L+1:end).'*Urlistx{istrap})',L), zeros(1,nhorizon+1)];  
           Fstrapsfully = [reconstructSSA_partial(Urlisty{istrap}*(Hankely(:,end-L+1:end).'*Urlisty{istrap})',L), zeros(1,nhorizon+1)];  
           Fstrapsfullz = [reconstructSSA_partial(Urlistz{istrap}*(Hankelz(:,end-L+1:end).'*Urlistz{istrap})',L), zeros(1,nhorizon+1)];  
           for ii = L : L+nhorizon
                for jj = 1:L-1
                    Fstrapsfullx(ii) = Fstrapsfullx(ii) + Philistx(istrap,jj)*Fstrapsfullx(ii-jj);
                    Fstrapsfully(ii) = Fstrapsfully(ii) + Philisty(istrap,jj)*Fstrapsfully(ii-jj);
                    Fstrapsfullz(ii) = Fstrapsfullz(ii) + Philistz(istrap,jj)*Fstrapsfullz(ii-jj);
                end
           end
           Ex_obs_forecast(istrap,:,:) = [Fstrapsfullx(:,L:end);Fstrapsfully(:,L:end);Fstrapsfullz(:,L:end)];         
        end           
    end
    nconverge = 1; trust = beta;
    sqrt_ak_stack = []; 
    while nconverge <= Nconvergemax      
        if isBootstrap == 1 % with strap (back-up) 
            for iii = 1:nhorizon+1
                 pbar_k = x_bar(xyzINDEX,iii);    
                 Alists = -(pbar_k- Ex_obs_forecast(:,:,iii)');
                 Blists = diag((pbar_k'- Ex_obs_forecast(:,:,iii))*Ex_obs_forecast(:,:,iii)')' + rbuff*vecnorm(pbar_k- Ex_obs_forecast(:,:,iii)',2);
                 ABlist = Alists.*Blists;
                 Ak_matrix = mean(Alists');
                 Bk_matrix = mean(Blists');
                 Ck_matrix = cov(Alists');
                 invCk = inv(Ck_matrix);
                 sqrtmCk = sqrtm(Ck_matrix);
                 Dk_matrix = mean(ABlist') - Ak_matrix*Bk_matrix; 
                 Ek_matrix = cov(Blists');
                 Kk_matrix = Ek_matrix - (Dk_matrix*invCk*Dk_matrix');
                 if Kk_matrix <=0
                     sqrt_Kk_matrix = 0;
                 else
                     sqrt_Kk_matrix = sqrt(Kk_matrix);
                 end
                 Hk_matrix = -invCk*Dk_matrix';
                 Aout_stack(:,iii) = [Ak_matrix';ones(3,1)*v_epsilon];
                 bineqvecstacks(:,iii) = [-Bk_matrix - v_epsilon*sqrt_Kk_matrix*sqrt(3);sqrtmCk*Hk_matrix;-sqrtmCk*Hk_matrix];
                 sqrt_ak_stack(iii,:,:) = sqrtmCk;
            end
            
            inputsstrap{1} = q0;
            for iii = 1: nhorizon+1
                inputsstrap{1+iii} = r_ref(iii+icounter,:)';
            end  
            for iii = 1:nhorizon+1
                inputsstrap{2+nhorizon+iii} = Aout_stack(:,iii);
            end    
            for iii = 1:nhorizon+1
                inputsstrap{3+2*nhorizon+iii} = squeeze(sqrt_ak_stack(iii,:,:));
            end   
            for iii = 1:nhorizon+1
                inputsstrap{4+3*nhorizon+iii} = bineqvecstacks(:,iii);
            end   

            for iii= 1: nhorizon+1
                inputsstrap{5+4*nhorizon+iii} = x_bar(:,iii);
            end   
            inputsstrap{6+5*nhorizon+1} = alpha*trust;
            [solutions,diagnostics] = controller_straps{inputsstrap};   
        else % no strap (back-up)
            rvec = (x_bar(xyzINDEX,:) - squeeze(Ex_obs_forecast(1,:,:)))./vecnorm(x_bar(xyzINDEX,:)-squeeze(Ex_obs_forecast(1,:,:)),2);
            Bk_matrix = rbuff + diag(rvec'*squeeze(Ex_obs_forecast(1,:,:)));
            inputs_nostrap{1} = q0;
            for ii = 1: nhorizon+1
                inputs_nostrap{1+ii} = r_ref(ii+icounter,:)';
            end   
                for iii = 1:nhorizon+1
                 inputs_nostrap{2+nhorizon+iii} = -rvec(:,iii);
                end    
                for iii = 1:nhorizon+1
                 Bk_matrix = rbuff + rvec'*Ex_obs_forecast(1,:,iii)';
                 inputs_nostrap{3+nhorizon*2+iii} = -Bk_matrix(iii);
                end   
                for iii= 1: nhorizon+1
                inputs_nostrap{4+nhorizon*3+iii} = x_bar(:,iii);
                end   
                inputs_nostrap{5+4*nhorizon+1} = alpha*trust;
                [solutions,diagnostics] = controller_nostraps{inputs_nostrap};   

        end
    if diagnostics~=0
        trust = trust/(beta); 
        looprelaxcounter(icounter) = looprelaxcounter(icounter) +1;
        if looprelaxcounter(icounter) > 10
            iMonteCarlo;
            if icounter < NSSA_training
                omittrun_badtrajectory = 1;
                break;
            else
                sucessflat(iMonteCarlo) = 0;    
                failedconverge = 1
                break;
            end
        end
    else
        u_list = solutions{1};
        x_bar = solutions{2};
        ucomd = u_list(:,1);
        trust = trust*beta;
        nconverge = nconverge+1;
    end 
    end
    tictoclist(icounter) = toc;
        if omittrun_badtrajectory == 1
            break;
        elseif  failedconverge == 1
           ucomd = u_bar(:,2);
           u_bar = [u_bar(:,2:end),zeros(4,1)];
           x_bar = [x_bar_save(:,2:end),x_bar(:,end)];
           break;
        else
           u_bar = u_list;
        end
   tspan = [ttraj(icounter), ttraj(icounter)+dt];
   [tode_out,qode_out] = ode45(@(t,qin) quad_dyn(t,qin,ucomd+Ueq,Ac,Bc), tspan, q0);
   qode_plot = [qode_plot;qode_out(end,:)];
   q0 = qode_out(end,:)';
   uhist = [uhist,ucomd+Ueq];
   tode = [tode;tode_out(end)];
   x_bar = [q0,x_bar(:,3:end),x_bar(:,end)];
   x_bar_save = x_bar;
   obsdist(icounter) = norm(q0(xyzINDEX) - xtraj(icounter+1,:)',2);
   if norm(q0(xyzINDEX) - xtraj(icounter+1,:)',2) <= rbuff
        sucessflat(iMonteCarlo) = 0;
        crashed = 1;
        iMonteCarlo = iMonteCarlo+1;
        break;
   end

   icounter = icounter + 1;
        end
    discounterlist{iMonteCarlo} =  obsdist;
    if omittrun_badtrajectory ~= 0
        omittrun_badtrajectory;
    else
        if crashed == 1
        sucessflat(iMonteCarlo) = 0;
        iMonteCarlo = iMonteCarlo + 1;
        elseif failedconverge == 1
            sucessflat(iMonteCarlo) = 0;
            failedconverge_list(iMonteCarlo) = 1;
        else
            sucessflat(iMonteCarlo) = 1;
            sucesstimeinfo(iMonteCarlo,:) = [mean(tictoclist),std(tictoclist)];    
            qode_List(iMonteCarlo,:,:) = qode_plot(:,xyzINDEX);        
            uode_List(iMonteCarlo,:,:) = [zeros(1,4);uhist'];  
            r_ref_List(iMonteCarlo,:,:) = r_ref; 
        end    
        iMonteCarlo = iMonteCarlo + 1;
    end
    iMonteCarlo
end
delete(gcp('nocreate'))
percent_sucess = sum(sucessflat)*100/nMonteCarlo
%%
save('MC_result_conspd_1000Run_ep05.mat','percent_sucess','sucesstimeinfo','qode_List','uode_List','r_ref_List',...
    'discounterlist','failedconverge_list','indexcollosion_list','-v7.3');
plotteron = 0;
%%
if plotteron == 1
%%
figure();
for iMonteCarlo = 1:nMonteCarlo
qode_plotplot = squeeze(qode_List(iMonteCarlo,:,:));
noisetraj_plotplot = squeeze(noise_xtraj_List(iMonteCarlo,:,:));
    for jlist = 1:3
        subplot(nMonteCarlo,3,jlist+(iMonteCarlo-1)*3);
        plot(noisetraj_plotplot(:,jlist),'--r');
        hold on;
        plot(qode_plotplot(:,jlist),'k');
    end
end

figure();
for iMonteCarlo = 1:nMonteCarlo
qode_plotplot = squeeze(qode_List(iMonteCarlo,:,:));
traj_plotplot = squeeze(xtraj_List(iMonteCarlo,:,:));
ref_plotplot = squeeze(r_ref_List(iMonteCarlo,:,:));
    for jlist = 1:3
        subplot(nMonteCarlo,3,jlist+(iMonteCarlo-1)*3);
        plot(traj_plotplot(:,jlist),'--r');
        hold on;
        plot(qode_plotplot(:,jlist),'k');
        hold on;
        plot(ref_plotplot(:,jlist),'--g');
    end
end
%% plot single traj with animator
nrunplot =1;
noise_xtraj = squeeze(noise_xtraj_List(nrunplot,:,:));
qode_plot = squeeze(qode_List(nrunplot,:,:));
maxplotrange = max(qode_plot) + 10;
minplotrange = min(qode_plot) - 10;
[Xsphere,Ysphere,Zsphere] = sphere;
robs = 0.25; ragent = 0.1;
figure();
for i = 100:1:N_sim_end - nhorizon-1
axis equal
surf(Xsphere*robs+noise_xtraj(i,1),Ysphere*robs+noise_xtraj(i,2),Zsphere*robs+noise_xtraj(i,3))
axis([minplotrange(1) maxplotrange(1) minplotrange(2) maxplotrange(2) minplotrange(3) maxplotrange(3)])
hold on;
surf(Xsphere*ragent+qode_plot(i,1),Ysphere*ragent+qode_plot(i,2),Zsphere*ragent+qode_plot(i,3))
pause(0.1)
drawnow
hold off;
end
end


%% post processing 
MC_postprocessing_flag = 0;
if MC_postprocessing_flag == 1
    %%
    countcount = 1;
    for i = 1:nMonteCarlo
        if isempty(discounterlist{i})
        else
                    mindist(countcount) = min(discounterlist{i});
                    countcount = countcount+1;
        end
    end
    meanmindist = mean(mindist)
    stdmindist = std(mindist)
end



