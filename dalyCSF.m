function output = dalyCSF(u,teta,l,i2,d,e)

% Change in bandwidth due to accomodation level
bwa    = 0.856 * d^0.14;

% Change in bandwidth due to eccentricity
k      = 0.6 * 0.24;
bwe    = 1/(1+k*e);

% Change in bandwidth due to orientation
ob     = 0.78;
bwteta = ((1 - ob)/2)*cos(4*teta) + (1 + ob)/2;

% the absolute peak sensitivity of the CSF
P      = 250 * bwe;

    function S = sensitivity(u,l,i2)
        
        A1 = 0.801 * (1+0.7/l)^-0.2;
        B1 = 0.3 * (1 + 100/l)^0.15;
        epsilon = 0.9;
        
        S  = ((3.23 .* (u.^2 .* i2).^-0.3).^5 + 1).^-0.2 ...
             .* A1 .* epsilon .* u .* exp(-B1 .* epsilon .* u)...
             .* sqrt(1 + 0.06 .* exp(B1 .* epsilon .* u)); 
    end

output = P * sensitivity(u./(bwa*bwe*bwteta),l,i2);

end
