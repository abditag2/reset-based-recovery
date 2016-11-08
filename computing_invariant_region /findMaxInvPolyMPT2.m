function P = findMaxInvPolyMPT2(n, m, A, B,  F, g, Hu, maxD, minD, max_it)
  %% Prepare iteration
  opt.abs_tol=1e-10;
  opt.rel_tol=1e-10;
  opt.verbose=2;
  opt.facecolor=[0.2 0.4 0.2];
  %opt.projection=3;

  plotting=1;


  % input constraints
  hu=ones(2*m,1);


  %one time to construct the initial constraints on the state space variables
  P = polytope(F, g)
  for i = 1:max_it
 
    disp(['iteration ',num2str(i)])

%     figure 
%     projChi = P.projection([1, 2, 4]);
%     plot(projChi)
    
%     [xx,yy] = double(projChi);
%     pppp = polytope(xx,yy)
%     pppp.plot();

    
    [H h]=double(P);
    clear PP;
%     %%fardins code
%     HAB = H* [A, B];
%     size_HAB = size(HAB); 
%     size_zeros_filling = [2*m, size_HAB(2)-m];
%     poly_before = polytope([H*[A B]; [zeros(size_zeros_filling), Hu]],[h;hu]); 
%     PP = projection(poly_before ,1:n,opt);
%     PP = intersect(PP,P)

    
    %%% mathiases code for \tau_r
    Ht=blkdiag(H,Hu);
    ht=[h;hu];
    PPT = polytope(Ht*[A B; zeros(m,n) eye(m)], ht)
    
%     %%% mathiases code for \tau_r and 2\tau_r
%     Ht=blkdiag(H, H, Hu);
%     ht=[h; h; hu];
%     PPT = polytope(Ht*[A*A A*B+B; A B; zeros(m,n) eye(m)], ht)

% for \tau_r and \tau_c
%     Ht=blkdiag(H, H, Hu);
%     ht=[h; h; hu];
%     PPT = polytope(Ht*[A_tc B_tc; A B; zeros(m,n) eye(m)], ht)

% for \2tau_r
%     Ht=blkdiag(H,Hu);
%     ht=[h;hu];
%     PPT = polytope(Ht*[A_tc B_tc; zeros(m,n) eye(m)], ht)


%     PPT2 = toPolyhedron(PPT);
%     PPT2.computeVRep();
    
%     PP=projection(PPT,1:n);
    PP = PPT.projection(1:n)
    

%     PP.computeHRep()
%     PP.A
%     PP = polytope(PP.A, PP.b);
    
    PP=intersect(P,PP) %CPre(P) * P  
    
    
    if isempty(PP) || ~isfulldim(PP)
        display('The resulting polytope is empty');
        break;
    end
    
    if le(P,PP,opt)
      P=PP;
      disp(['P_i \subset P_i+1'])
 %     if plotting
 %       plot(projection(PP,[1 2 3]),opt);
 %       drawnow
 %     end
       break;
    end 
    P=PP;
    toPolyhedron(P)
  end
end

