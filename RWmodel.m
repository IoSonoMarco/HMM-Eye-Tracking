function o = RWmodel(par, var, array)

beta = par(1);
pi_par = 0;
aa = [0 par(2)];
bb = [0 par(3)];

init = [0.5 0.5];

LL = [];

for subj = 1:length(array)
    
    L = [];
    
    y = array(subj).AFO*10;
    x = array(subj).target_pos; 
    z = array(subj).target_ind;
    
    for session = 1:size(z,1)
        
        Y = y(session,:);
        nan_pos = find(isnan(Y));
        X = x(session,:); 
        Z = z(session,:);
        T = length(Z);
        
        if sum(nan_pos == 1) == 1
            diag_ = diag([1 1]);
        else
            diag_ = diag([truncatedVMpdf(Y(1),0,var(1),0,pi) truncatedVMpdf(Y(1),0,var(2),-pi,0)]);
        end
        
        alpha = init * diag_;
        phi = alpha ./ (sum(alpha)+1e-5);
        l = log(sum(alpha)+1e-5);
        
        omega(1) = 0.5;
        
        G = [0.5 0.5;
             0.5 0.5];
        
        if sum(nan_pos == 2) == 1
            diag_ = diag([1 1]);
        else
            diag_ = diag([truncatedVMpdf(Y(2),0,var(1),0,pi) truncatedVMpdf(Y(2),0,var(2),-pi,0)]);
        end
        
        phi = phi * G * diag_;
        l = l + log(sum(phi)+1e-5);
        phi = phi ./ (sum(phi)+1e-5);
         
        for t = 2:(T-1)    
            
            if sum(nan_pos == t) == 1
                Y_ = 0;
            else
                Y_ = Y(t);
            end
            
            ALT(t) = X(t) + X(t-1) - 2*X(t)*X(t-1);
            omega(t) = (1-ALT(t))*beta + (1 - beta)*omega(t-1);
            omega(t) = max([0.001 min([omega(t) 0.999])]);
            g_21(t) = tanh( aa(2)*(omega(t) + aa(1))/(1 - omega(t)) );
            g_12(t) = tanh( bb(2)*(1 - omega(t) + bb(1))/omega(t) );
            G = [1-g_12(t) g_12(t);
                 g_21(t) 1-g_21(t)]; 
            
            if sum(nan_pos == t+1) == 1
                diag_ = diag([1 1]);
            else
                diag_ = diag([truncatedVMpdf(Y(t+1),0,var(1),0,pi) truncatedVMpdf(Y(t+1),0,var(2),-pi,0)]);
            end
            
            phi = phi * G * diag_;
            l = l + log(sum(phi)+1e-5);
            phi = phi ./ (sum(phi)+1e-5);        
        end
         
        L(session) = l;        
    
    end
    
    LL(subj) = sum(L);

end

O = -sum(LL);
if isnan(O) | ~isreal(O)
    oo = 999999;
else
    oo = O;
end

o = oo;

end
