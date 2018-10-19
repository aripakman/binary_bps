 classdef MRF < handle  % inherit from handle so that we can pass by reference

        properties(SetObservable = true)
        % These properties are public by default
           d;
           c1;
           c2;
           mlp=-Inf;    %this is a running 
           M; 
           r;           

        end
                
        %properties(GetAccess = 'public', SetAccess = 'private')
        % These properties can be accessed but not set from outside the class
        
      %  end
        
             
        methods
        % These methods are public by default. 
        
            function obj = MRF(d,c1,c2, frust)
            % class constructor            
                obj.c1 = c1;
                obj.c2 = c2;
                obj.d = d;
                m = c2*normrnd(0,1,d);
                obj.M = (m+m')/2;
                for i=1:d
                    obj.M(i,i)=0;
                end
                if ~frust
                    obj.r = c1*normrnd(0,1,d,1);                
                else
                    me = sum(obj.M,1);
                    obj.r = c1*normrnd(me',2);                
                end
            end
            
            
            % Conventions:
            % B is a binary vector (0,1)
            % S is a signs vector  (+1,-1)
            
             function lp = logp(obj,S)                 
                B = (S + ones(obj.d,1))/2;             
                lp = 0.5*B'*obj.M*B + obj.r'*B;
                if lp > obj.mlp
                    obj.mlp = lp;
                end                
                    
             end
             
             
             
             function lc = logp_change(obj,S,j)
                 S(j) = +1;
                 lc = obj.logp(S);
                 S(j) = -1;
                 lc = lc - obj.logp(S);
             end
             
         
             function sig2 = sigma2(obj,S)
                
                lp = obj.logp(S)-obj.mlp;
                fS = exp(lp);
                sig2 = 2/pi*(fS)^(2/obj.d);
                
             end
             
             
             function p = pprob(obj,B)
                 
                lp = -0.5*B'*obj.M*B - obj.r'*B;
                p=exp(lp);
             end
                 

             
        end
        
        
        
end