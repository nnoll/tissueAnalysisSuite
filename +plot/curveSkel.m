function [ ] = curveSkel( Struct, color, mem )
%INFERCURVE Summary of this function goes here
%   Detailed explanation goes here

    if (nargin == 1)
        color = 'b';
    end

   hold on
   for b = 1:length(Struct.Bdat)
       if (length(Struct.Bdat(b).verts) == 2 && sum(Struct.Bdat(b).verts > 0) == 2)
           if (Struct.Bdat(b).radius < inf)
               r = double([[Struct.Vdat(Struct.Bdat(b).verts).vertxcoord]; ...
                    [Struct.Vdat(Struct.Bdat(b).verts).vertycoord]]);
               r = r';
               rho = Struct.Bdat(b).rBar;
               S = sign( (r(1,1) - rho(1))*(r(2,2) - rho(2)) - (r(1,2) - rho(2))*(r(2,1) - rho(1)) );
               if ( sign(S) < 0 )
                   r = r([2,1],:);
               end
               R = Struct.Bdat(b).radius;

               theta = atan2(r(:,2)-rho(2),r(:,1)-rho(1));
               theta(theta<0) = theta(theta<0) + 2*pi;
               theta = sort(theta);
               if (theta(2) - theta(1) > pi)
                   theta(2) = theta(2)-2*pi;
                   theta = theta([2,1]);
               end
               theta = linspace(theta(1),theta(2),20);
               x = rho(1) + R*cos(theta);
               y = rho(2) + R*sin(theta);
               plot(x,y,'LineWidth',2,'Color',color)
 
           else
               r = [[Struct.Vdat(Struct.Bdat(b).verts).vertxcoord]; ...
                    [Struct.Vdat(Struct.Bdat(b).verts).vertycoord]];
               r = r';
               plot(r(:,1),r(:,2),'LineWidth',2,'Color',color)
           end
       end
   end
    
end

