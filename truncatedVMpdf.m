function pdf = truncatedVMpdf(theta,mu,k,a,b)

num = exp(k .* cos(theta-mu));
den = integral(@(theta)( exp(k .* cos(theta-mu)) ),a,b);

if theta < a
    pdf = 0;
elseif theta > b
    pdf = 0;
else
    pdf = num/den;
end

end