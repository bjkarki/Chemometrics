function [rank, stdrank] = draug (d)
% d = data matrix

[r, c] = size (d);

sm = c; lg = r;
if r < c 
    sm = r; lg = c; 
end

amt = 0.015 * sqrt (sum (diag (d' * d)) / (r * c));

for ii = 1 : 100
    dx = d + amt * rand (r, 1) * rand (1, c);
    [u, s, v] = svd (d, 0);
    [ux, sx, vx] = svd (dx, 0);
    
    for j = 1 : sm
        ev (j) = s (j, j) * s (j, j);
        evx (j) = sx (j, j) * sx (j, j);
    end
    
    for k = 1 : sm
        sev (k) = sum (ev (k : sm));
        sevx (k) = sum (evx (k : sm));
    end
    
    df (1) = lg * sm;
    var (1) = sev (1) / df (1);
    varx (1) = sevx (1) / df (1);
    
    for i = 1 : (sm - 1)
        df (i + 1) = (lg - i) * (sm - i);
        var (i + 1) = sev (i + 1) / df (i);
        varx (i + 1) = sevx (i + 1) / df (i);
    end
    
    for ij = 1 : (sm - 1)
        f (ij) = var (ij) / varx (ij + 1);
        sl (ij) = 100 * (1 - fcdf (f (ij), df (ij), df (ij + 1)));
    end
    
    fl = 1 : (sm - 1);
    ranki (ii) = 0;
    
    for j = 1 : (sm - 1)
        if sl (j) > 1, ranki (ii) = j - 1; break, end
    end
end

rank = mean (ranki)
stdrank = std (ranki)

end
        