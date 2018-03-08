function varreset_optimize()


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fit two-stage model with imperfect reset
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % BEGIN FUNCTION FIT_TWOSTAGE_VAR_RESET


    for i=1:(length(P_noreset))
        startOptimization=tic;
        %sets initial guess
        x0 = P_noreset(i,:);
        fprintf("optimizing seed %d: m1=%f s1=%f m2=%f s2=%f r=%f\n", i, x0(1),x0(2),x0(3),x0(4),x0(5));
        f=@(x, m1,s1,m2,s2,r)convolv_2invG_var_reset(x,m1,s1,m2,s2,r,.01);

        fminsearch_options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',1, 'TolX', 0.01);
        myll=@(params)loglikelihood(datatrain, f, 5, params);
        objfun=@(params)penalize(myll, 5, params, [realmin  realmax;realmin  realmax;realmin  realmax;realmin  realmax;0.001  0.999])
        p=fminsearch(objfun,x0,fminsearch_options);

        %[p,conf]=mle(datatrain,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf 0.9],'lowerbound',[0 0 0 0 0],'options',options);

        fprintf("optimized: m1=%f s1=%f m2=%f s2=%f r=%f\n", p(1),p(2),p(3),p(4),p(5));
        %save parameters
        pd_varreset(i,:,kk)=p;
        %gets the likelihhod of the parameters
        [l,hp,flag_varreset(kk,i),E]=convolv_2invG_var_reset(datatrain,p(1),p(2),p(3),p(4),p(5),.01);
        l=sum(log(l));
        fprintf("log-liklihood=%f\n",l);
        if flag_varreset(i) == 1
            fprintf("used dirac delta approximation in final result\n");
        end
        ld_varreset(kk,i)=l;
        toc(startOptimization)
        fprintf("\n");
    end

            % we previously optimized with a larger step size, recalculate with
            % a smaller stepsize after the fact
    %         ld_true_noreset=zeros(length(ld),1);
    %         for i=1:length(ld)
    %             [l,hp_true(i),flag_true_noreset(i),E_true(i)]=convolv_2invG_noreset(datatrain,pd_noreset(i,1),pd_noreset(i,2),pd_noreset(i,3),pd_noreset(i,4),pd_noreset(i,5),.001);
    %             ld_true_noreset(i)=sum(log(l));
    %         end

    % skip recalculation for now

    ld_true_varreset=ld_varreset(kk,:);
    flag_true_varreset=flag_varreset(kk,:);

    %find the best flagged model and the best model that is not flagged

    %indices of flagged models
    indflag_varreset=find(flag_true_varreset==1);

    %ibndices of no flag models
    indnoflag_varreset=find(flag_true_varreset==0);


    ld_trueflag_varreset=ld_true_varreset(indflag_varreset);
    ld_true_varreset=ld_true_varreset(indnoflag_varreset);

    pflag_varreset=pd_varreset(indflag_varreset,:);
    p_varreset=pd_varreset(indnoflag_varreset,:);

    % common to each fit, consider factoring out
    [max_ld_varreset(kk,1),row_ld_varreset(kk,1)]=max(ld_true_varreset)

    %best nonflagged model
    pd_max_varreset(kk,:) = p_varreset(row_ld_varreset(kk,1),:,kk)

    %cross validate
    [ll]=convolv_2invG_varreset(datacross,pd_max_varreset(kk,1),pd_max_varreset(kk,2),pd_max_varreset(kk,3),pd_max_varreset(kk,4),pd_max_varreset(kk,5),.01);
    lcross_varreset(kk,1)=sum(log(ll));

    if isempty(indflag_varreset)==0
        [max_ldflag_varreset(kk,1),row_ldflag_varreset(kk,1)]=max(ld_trueflag_varreset)

        %best flagged model
        pd_maxflag_varreset(kk,:) = pflag_varreset(row_ldflag_varreset(kk,1),:,kk)

        %perform cross validation
        [ll]=convolv_2invG_varreset(datacross,pd_maxflag_varreset(1),pd_maxflag_varreset(2),pd_maxflag_varreset(3),pd_maxflag_varreset(4),pd_maxflag_varreset(5),.01);
        lcrossflag_varreset(kk,1)=sum(log(ll));
    end

    % END FUNCTION FIT_TWOSTAGE_VAR_RESET

end