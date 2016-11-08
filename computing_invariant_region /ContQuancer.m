
%%
format long
clc;
clear;
clf;


%%%%%%%%%%%%%%%%%%%% Constants Definition

% maximum iterations for finding the invariant set
MAX_ITERATIONS_INVARIANT_REGION = 100

%state variables: 
n = 4;
%control inputs
m = 2;  
% number of disturbance variables
d = 1;

RESTART_TIME = 25;
FIRST_POINT_INVARIANCE = 5;
DT_CONTINUOUS = 0.010

N_STEPS_FOR_COMPUTING_READJUSTED_SAFETY_REGION = 5

STEPS_FOR_CONTROLLED_OUTCOME = 100  

%%%%%%%%%%%%%%%%%%%% Computing the timing parameters of the System


SECOND_POINT_INVARIANCE= RESTART_TIME + FIRST_POINT_INVARIANCE;

delta_t = RESTART_TIME * DT_CONTINUOUS;

t1 = FIRST_POINT_INVARIANCE * DT_CONTINUOUS;
t2 = SECOND_POINT_INVARIANCE * DT_CONTINUOUS;


%%%%%%%%%%%%%%%%%%%% Defining the Safety Constraints

%limits on the state
% ELEVATION_MAX = 0.2
% ELEVATION_MIN = - ELEVATION_MAX

THIRD_ELEV_PLUS_PITCH_MIN = 0.3

% PITCH_MAX = 1.2
% PITCH_MIN = -PITCH_MAX

% TRAVEL_ANGEL_MAX = 1/4
% TRAVEL_ANGEL_MIN = -TRAVEL_ANGEL_MAX

D_ELEVATION_MAX = 0.4
D_ELEVATION_MIN = -D_ELEVATION_MAX

D_PITCH_MAX = 1.25
D_PITCH_MIN = -D_PITCH_MAX
% 
% D_TRAVEL_ANGEL_MAX = 1
% D_TRAVEL_ANGEL_MIN = -D_TRAVEL_ANGEL_MAX

%limits on the control variables
U1_MAX = 1.1
U1_MIN = -U1_MAX

U2_MAX = U1_MAX
U2_MIN = -U2_MAX

%bounds on the disturbance
maxD = [1] 
minD = [-1]

%Definition of boundaries
%bounds on the state
F = [
     -1  -1/3    0   0   ;
%      1   0       0   0   ;
%      0   1       0   0   ;
     0   0       1   0   ;
     0   0       0   1   ;

     -1  +1/3    0    0   ;
%      -1  0       0    0  ;
%      0   -1      0    0  ;
     0   0       -1   0   ;
     0   0       0   -1   ;
     ] %should be 2n*m

%elements of acceleration in g cannot be large because then the freedom is
%limited in the permitted angle of operation for the system.
g = [
        THIRD_ELEV_PLUS_PITCH_MIN;
%         ELEVATION_MAX
%         PITCH_MAX;
        D_ELEVATION_MAX;
        D_PITCH_MAX;

        THIRD_ELEV_PLUS_PITCH_MIN;
%         -ELEVATION_MIN;
%         -PITCH_MIN;
        -D_ELEVATION_MIN;
        -D_PITCH_MIN;
        ] %should be 2n*1
    
    
%bounds on the control input

Hu = [  1/U1_MAX,  0;
        -1/U1_MAX, 0;
        0,      1/U2_MIN;
        0,      -1/U2_MIN] %should be 2*m * m


%%%%%%%% Model of the Quanser System

%define model parameters of the quancer system Identified from experimental
%data and manual of the system
par = [-1.0000   -2.4000   -0.0943    0.1200    0.1200   -2.5000   -0.0200    0.200    2.1000   10.0000];

p1 = par(1);
p2 = par(2);
p3 = par(3);
p4 = par(4);
p5 = par(5);
p6 = par(6);
p7 = par(7);
p8 = par(8);
p9 = par(9);
p10 = par(10);

A_c = [0  0   1 0 ;
       0  0   0 1 ; 
       0  0   0 0 ; 
       0  0   0 0 ; ];
 
B_c = [0    0; 
       0    0; 
       p8   p8; 
       p9   -p9; ];


sys1=c2d(ss(A_c,B_c,eye(n),zeros(n,m)),t1);

A_t1=sys1.a;
B_t1=sys1.b;

syst2=c2d(ss(A_c,B_c,eye(n),zeros(n,m)),t2);

A_t2=syst2.a;
B_t2=syst2.b;
%% this section is for paper
format short
%state variables
n_com = 6;
%control inputs
m_com = 2;  

K_f = 0.1188;
L_a = 0.660;
L_w = 0.470;
L_h = 0.178;
M_f = 0.713;
M_w = 1.87;
g_param = 9.81;

pA62 = - (L_w*M_w - 2*L_a*M_f)*g_param / (M_w * L_w^2 + 2* M_f * L_h^2 + 2*M_f*L_a^2)
pB41 = L_a * K_f / (2*M_f * L_a^2 + M_w * L_w^2)
pB42 = L_a * K_f / (2*M_f * L_a^2 + M_w * L_w^2)
pB51 = 1/2 * K_f / (M_f *L_h)
pB52 = -1/2 * K_f / (M_f *L_h)

A_com=[0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 pA62 0 0 0 0;]
 
 
B_com=[0 0;
     0 0;
     0 0;
     pB41 pB42;
     pB51 pB52;
     0 0;]

sys1_com=c2d(ss(A_com,B_com,eye(n_com),zeros(n_com,m_com)),delta_t);

A_d_com=sys1_com.a
B_d_com=sys1_com.b

%% Find The Readjusted Safety Region
g_temp = g

for kk = 1:N_STEPS_FOR_COMPUTING_READJUSTED_SAFETY_REGION
    kk
    g_temp = findCompressedInvRegion(F, g_temp, A_c, B_c, Hu, delta_t/N_STEPS_FOR_COMPUTING_READJUSTED_SAFETY_REGION);
end

gCompressed = g_temp


%% Find the Invariant Region

% C = findMaxInvPolyMPT2(n, m, A, B, F, gCompressed, 10/8*Hu, maxD, minD, MAX_ITERATIONS_INVARIANT_REGION)
% C = findMaxInvPolyMPT2(n, m, A, B, F, g, 10/8*Hu, maxD, minD, MAX_ITERATIONS_INVARIANT_REGION)
% C = findMaxInvPolyMPT2(n, m, A, B, F, g, Hu, maxD, minD, MAX_ITERATIONS_INVARIANT_REGION)

hx = gCompressed
Hx = F

  opt.abs_tol=1e-10;
  opt.rel_tol=1e-10;
  opt.verbose=2;
  opt.facecolor=[0.2 0.4 0.2];

for i = 1:10
    
    C = findMaxInvPolyMPT2(n, m, A_t1, B_t1,  Hx, hx, Hu, maxD, minD, MAX_ITERATIONS_INVARIANT_REGION)
    
    [Hx hx]=double(C)
    P1 = polytope(Hx, hx)
    
    C = findMaxInvPolyMPT2(n, m, A_t2, B_t2,  Hx, hx, Hu, maxD, minD, MAX_ITERATIONS_INVARIANT_REGION)
    [Hx hx]=double(C)
    P2 = polytope(Hx, hx)
    
       if le(P2,P1,opt)
           hello = 12
           break
       end
end

[Hx hx]=double(C)

%% test the controller with the computed invariant region
% Uncomment the following for testing the Invariant region

% 
% DISC_SIM_MODE = 0;
% 
% CPolyhedron = Polyhedron(Hx, hx) %computes extreme points of a polytobe
% V = CPolyhedron.V;
% 
% Chi = polytope(Hx, hx);
% % V = lcon2vert(Hx,hx);
% 
% dima = 1;
% dimb = 2;
% 
% size_v = size(V,1)
% 
% for i = 1:size_v
% 
%     disp(['************************************************ vertice number: ', num2str(i)])
%     x_next = 0.99999*V(i,:)';
% 
%     figure 
%     projChi = Chi.projection([dima, dimb]);
%     projChi.plot();
% 
%     hu=ones(2*m,1)
%     H1=[Hx*B; Hu];
% 
%     hold on
%     u= [];
%     
%     for j=0:100
%         
%       plot(x_next(dima),x_next(dimb),'*b')
%       x_curr = x_next;
% 
%       if DISC_SIM_MODE == 1
%           
%           h=[hx-Hx*A*x_curr; hu];
%           U = lcon2vert(H1,h)
%           size_u = size(U);
%           if (size_u(1) > 0)
% 
%             if (j < STEPS_FOR_CONTROLLED_OUTCOME)
%                 j;
%                 x_curr;
%                 u = [];
% 
%                 u_mid = (max(U) + min(U))/2;
%                 x_temp = A*x_curr + B*u_mid';
%                 h_temp=[hx-Hx*A*x_temp; hu];
%                 U_temp = lcon2vert(H1,h_temp);
%                 
%                 if size(U_temp,1) > 0
%                     display('mid')
%                    u = u_mid';
%                 end   
%                        
%                 u_max = max(U);
%                 x_temp = A*x_curr + B*u_max';
%                 h_temp=[hx-Hx*A*x_temp; hu];
%                 U_temp = lcon2vert(H1,h_temp);
%                 if size(U_temp,1) > 0
%                     display('max')
%                    u = u_max';
%                 end
% 
%                 u_min = min(U);
%                 x_temp = A*x_curr + B*u_min';
%                 h_temp=[hx-Hx*A*x_temp; hu];
%                 U_temp = lcon2vert(H1,h_temp);
%                 
%                 if size(U_temp,1) > 0
%                    display('min')
%                    u = u_min';
%                 end
% 
%             else
%                 u = min(U)/2 + max(U)/2;
%             end
%           else
%             disp('There is no contorl input!')
%             j
%             i
%             u
%             x_curr
%             return
%           end          
%       elseif DISC_SIM_MODE == 0
%           
%         if mod(j,DISCRETE_STEP_SIZE) == 0
%           h=[hx-Hx*A*x_curr; hu];
%           U = lcon2vert(H1,h);
%           size_u = size(U);
%           if (size_u(1) > 0)
% 
%             if (j < STEPS_FOR_CONTROLLED_OUTCOME*DISCRETE_STEP_SIZE)
%                 j;
%                 u;
%                 x_curr;
%                 u = [];
% 
% 
%                 u_mid = (max(U) + min(U))/2;
%                 x_temp = A*x_curr + B*u_mid';
%                 h_temp=[hx-Hx*A*x_temp; hu];
%                 U_temp = lcon2vert(H1,h_temp);
%                 if size(U_temp,1) > 0
% %                     display('mid')
%                    u = u_mid';
%                 end   
%                        
%                 u_max = max(U);
%                 x_temp = A*x_curr + B*u_max';
%                 h_temp=[hx-Hx*A*x_temp; hu];
%                 U_temp = lcon2vert(H1,h_temp);
%                 if size(U_temp,1) > 0
% %                     display('max')
%                    u = u_max';
%                 end
% 
%                 u_min = min(U);
%                 x_temp = A*x_curr + B*u_min';
%                 h_temp=[hx-Hx*A*x_temp; hu];
%                 U_temp = lcon2vert(H1,h_temp);
%                 
%                 if size(U_temp,1) > 0
% %                    display('min')
%                    u = u_min';
%                 end
% 
%             else
%                 u = (min(U)/2 + max(U)/2)';
%             end
%           else
%             disp('There is no contorl input!')
%             j
%             i
%             u
%             x_curr
%             return
%           end  
%           
%         end
% 
%       end 
%       
%       
%      if DISC_SIM_MODE == 1
%       %discrete simulation of the path
%        x_next = A*x_curr + B*u;
%      else
%       %continuous simulation of the path
%       x_next = x_curr + (A_c*x_curr + B_c*u)*DT_CONTINUOUS;
%      end 
%      
%      if (F*x_next < g ~= ones(size(F,1),1))
%          display('System went outside of the boundaries')
%             j
%             i
%             u
%             x_next
%             return
%          break
%      end
%      
%     end  
%     
% end
