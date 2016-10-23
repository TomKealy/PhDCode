eta = [0.6,0.3]
R =  [0.4, 0.6 ; 0.1 , 0.9]

t = eye(2) - R

tm = inv(t)

lambda = eta*tm 
