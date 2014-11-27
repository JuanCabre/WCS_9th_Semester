function [total]=simple_OPEX(vf,N,t_init,r,infl)

total=0;
for t=t_init:N
    total = total + vf*(1+infl/100)^t/((1+r/100)^t);
end
total
end