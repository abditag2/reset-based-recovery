% In order to run this 

% 1. Run ContQuanser with desired parameters first
% 2. Alternativelyt, in order to generate the figures of the paper, load 
% the matrices stored in the "matrices" folder to the workspace. 



Hp = [  0   0   1   0;
        0   0   -1  0;
        0   0   0   1;
        0   0   0   -1]

    hp = [  0.191;
            -0.19;
            0.181;
            -0.18]

fixedPlane = polytope(Hp, hp)

ChiPolyHedron = polytope(Hx, hx)

plane2 = intersect(fixedPlane, ChiPolyHedron)

[Hi hi]=double(plane2);

plane = polytope(Hi, hi);

projPlane = plane.projection([1, 2]);
projPlane.plot();