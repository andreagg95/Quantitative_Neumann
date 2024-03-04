function [z] = auxil(K2,c)
% characteristic function of the set 
% { (K2,c) \in \R^2 : -1 <= c <= 0, -c-1 <= K2 <= -c }

if (K2>-c || K2 < -c-1 || c>0 || c <-1)
    z = 0;
   return
else
    z = 1;
end

end
