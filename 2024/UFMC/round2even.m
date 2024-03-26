function result = round2even(x)
    resto = mod(x, 2);

    result = x;
    
    if resto > 0
        result = result + 1;
    end
