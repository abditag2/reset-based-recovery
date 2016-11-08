
% In order to run this 

% 1. Run ContQuanser with desired parameters first
% 2. Alternativelyt, in order to generate the figures of the paper, load 
% the matrices stored in the "matrices" folder to the workspace. 


SIM_STEP = 0.001
sys_sim=c2d(ss(A_c,B_c,eye(n),zeros(n,m)),SIM_STEP);

A_sim=sys_sim.a;
B_sim=sys_sim.b;


dim_names = {'$\epsilon$ (rad)'; '$\rho$ (rad)'; '$\dot{\epsilon}$ (rad/s)'; '$\dot{\rho}$ (rad/s)'};




CPolyhedron = Polyhedron(Hx, hx) %computes extreme points of a polytobe
% V = CPolyhedron.V;

% init_point = [-0.141034871577966,-0.00,-0.0281111111111133,0.0513333333333605]

init_point = [-0.141034871577966,0.00256666666666817,-0.0281111111111131,-0.0513333333333617]

coeff = 0.999

Chi = polytope(Hx, hx);
% V = CPolyhedron.V;
%%
dim_y = 1;
dim_x = 2;
hold on

subplot(4,4,14)

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')
%%
dim_y = 1;
dim_x = 3;
hold on
subplot(4,4,15)

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')
%%
dim_y = 1;
dim_x = 4;
hold on
subplot(4,4,16)

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')
%%
dim_y = 2;
dim_x = 3;
hold on
subplot(4,4,11)

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')
%%
dim_y = 2;
dim_x = 4;
hold on
subplot(4,4,12)

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')
%%
dim_y = 3;
dim_x = 4;
hold on

subplot(4,4,8)

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')
