function [aw,err] = fcn_match_strength(ar,w,s,temp,dec,energyfcn,tolerance,maxiter)
n = length(ar);
[u,v] = find(triu(ar,1));
indx = (v - 1)*n + u;
aw = zeros(n);
aw(indx) = w;
aw = aw + aw';
sw = sum(aw,2);
errvec = abs(sw - s);
err = sum(errvec);
m = length(u);
errbest = err;
wbest = w;
iter = 0;
while errbest > tolerance
    iter = iter + 1;
    x = randi(m);
    y = randi(m);
    
    ux = u(x);
    vx = v(x);
    uy = u(y);
    vy = v(y);
    
    swp = sw;
    swp([ux,vx]) = swp([ux,vx]) - w(x) + w(y);
    swp([uy,vy]) = swp([uy,vy]) + w(x) - w(y);
    
    errvecp = errvec;
    errvecp([ux,vx,uy,vy]) = abs(swp([ux,vx,uy,vy]) - s([ux,vx,uy,vy]));
    switch energyfcn
        case 'maxpercentchange'
            errp = max(errvecp./s);
        case 'meandiff'
            errp = mean(errvecp);
        case 'maxdiff'
            errp = max(errvecp);
    end
    if errp < err || rand < exp(-(errp - err)./temp)
        w([x,y]) = w([y,x]);
        err = errp;
        errvec = errvecp;
        sw = swp;
        if err < errbest
            errbest = err;
            wbest = w;
        end
    end
    temp = temp*dec;
    if mod(iter,10000) == 0
        fprintf('%i - %.3f\n',iter,errbest);
    end
    if iter > maxiter
        errbest = 0;
    end
end
aw = zeros(n);
aw(indx) = wbest;
aw = aw + aw';