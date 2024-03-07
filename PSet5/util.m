function u = util(x, sigma)
    if sigma == 1
        u=log(x);
    else 
        u = -x^(1-sigma)/(1-sigma);
    end
end