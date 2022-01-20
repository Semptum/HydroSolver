function invert_2(matrix)
    a,c,b,d=matrix
    return [d  -b 
             -c  a]/(a*d-b*c)
end