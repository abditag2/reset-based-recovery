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

coeff = 0.99

Chi = polytope(Hx, hx);
% V = CPolyhedron.V;
%%

dim_y = 1;
dim_x = 2;

disp(['************************************************ vertice number: ', num2str(i)])

x_next = coeff*init_point';
hold on


subplot(2,2,1)
% h = figure('position',[100 100 200 200])

% h = figure('position',[300 100 200 200])
% h = figure('position',[500 100 200 200])
% h = figure('position',[700 100 200 200])

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')



hold on
u= [];

Xs = []

for j=0:300
    
    Xs = [Xs; [x_next(dim_x) x_next(dim_y)]]
    if (j == 0 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',10)
    elseif ( j < 50 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',1)
    elseif ( j == 50)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',10)
    elseif ( j < 300)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',1)
    else
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',10)
    end
    x_curr = x_next;
    
    if (j == 0)
        hu=ones(2*m,1)
        
        H1=[Hx*B_t1; Hx*B_t2; Hu];
        h=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
        
        U = lcon2vert(H1,h)
        size_u = size(U);
        if (size_u(1) > 0)
            
            u_mean= mean(U);
            u = u_mean';
            %             if (j < STEPS_FOR_CONTROLLED_OUTCOME)
            %                 j;
            %                 x_curr;
            %                 u = [];
            %                 u_mid = (max(U) + min(U))/2;
            %                 x_temp = A_t1*x_curr + B_t1*u_mid';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('mid')
            %                     u = u_mid';
            %                 end
            %
            %                 u_max = max(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %                 if size(U_temp,1) > 0
            %                     display('max')
            %                     u = u_max';
            %                 end
            %
            %                 u_min = min(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('min')
            %                     u = u_min';
            %                 end
            %
            %             else
            %                 u = min(U)/2 + max(U)/2;
            %             end
            
        else
            disp('There is no contorl input!')
            j
            i
            u
            x_curr
            return
        end
    end
    
    %discrete simulation of the path
    %     x_next = x_curr + (A_c*x_curr + B_c*u)*DT_CONTINUOUS;
    x_next = A_sim*x_curr + B_sim*u;
    
    
    if (F*x_next < g ~= ones(size(F,1),1))
        display('System went outside of the boundaries')
        j
        i
        u
        x_next
        return
        break
    end
    
end



plot(Xs(:,1), Xs(:,2), '-b')

axis([-1e-3 8e-3 -0.1415 -0.134])

% %%
% dim_y = 1;
% dim_x = 3;
% 
% disp(['************************************************ vertice number: ', num2str(i)])
% 
% x_next = coeff*init_point';
% hold on
% 
% % h = figure('position',[100 100 200 200])
% h = figure('position',[300 100 200 200])
% % h = figure('position',[500 100 200 200])
% % h = figure('position',[700 100 200 200])
% 
% projChi = Chi.projection([dim_x, dim_y]);
% plot(projChi);
% xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
% ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')
% 
% 
% 
% hold on
% u= [];
% 
% Xs = []
% 
% for j=0:300
%     
%     Xs = [Xs; [x_next(dim_x) x_next(dim_y)]]
%     if (j == 0 )
%         plot(x_next(dim_x),x_next(dim_y),'-mo',...
%             'LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor',[1 1 1],...
%             'MarkerSize',10)
%     elseif ( j < 50 )
%         plot(x_next(dim_x),x_next(dim_y),'-mo',...
%             'LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor',[0 0 1],...
%             'MarkerSize',1)
%     elseif ( j == 50)
%         plot(x_next(dim_x),x_next(dim_y),'-mo',...
%             'LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor',[0 0 1],...
%             'MarkerSize',10)
%     elseif ( j < 300)
%         plot(x_next(dim_x),x_next(dim_y),'-mo',...
%             'LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor',[.49 1 .63],...
%             'MarkerSize',1)
%     else
%         plot(x_next(dim_x),x_next(dim_y),'-mo',...
%             'LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor',[.49 1 .63],...
%             'MarkerSize',10)
%     end
%     x_curr = x_next;
%     
%     if (j == 0)
%         hu=ones(2*m,1)
%         
%         H1=[Hx*B_t1; Hx*B_t2; Hu];
%         h=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
%         
%         U = lcon2vert(H1,h)
%         size_u = size(U);
%         if (size_u(1) > 0)
%             
%             u_mean= mean(U);
%             u = u_mean';
%             %             if (j < STEPS_FOR_CONTROLLED_OUTCOME)
%             %                 j;
%             %                 x_curr;
%             %                 u = [];
%             %                 u_mid = (max(U) + min(U))/2;
%             %                 x_temp = A_t1*x_curr + B_t1*u_mid';
%             %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
%             %                 U_temp = lcon2vert(H1,h_temp);
%             %
%             %                 if size(U_temp,1) > 0
%             %                     display('mid')
%             %                     u = u_mid';
%             %                 end
%             %
%             %                 u_max = max(U);
%             %                 x_temp = A_t1*x_curr + B_t1*u_max';
%             %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
%             %                 U_temp = lcon2vert(H1,h_temp);
%             %                 if size(U_temp,1) > 0
%             %                     display('max')
%             %                     u = u_max';
%             %                 end
%             %
%             %                 u_min = min(U);
%             %                 x_temp = A_t1*x_curr + B_t1*u_max';
%             %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
%             %                 U_temp = lcon2vert(H1,h_temp);
%             %
%             %                 if size(U_temp,1) > 0
%             %                     display('min')
%             %                     u = u_min';
%             %                 end
%             %
%             %             else
%             %                 u = min(U)/2 + max(U)/2;
%             %             end
%             
%         else
%             disp('There is no contorl input!')
%             j
%             i
%             u
%             x_curr
%             return
%         end
%     end
%     
%     %discrete simulation of the path
%     %     x_next = x_curr + (A_c*x_curr + B_c*u)*DT_CONTINUOUS;
%     x_next = A_sim*x_curr + B_sim*u;
%     
%     
%     if (F*x_next < g ~= ones(size(F,1),1))
%         display('System went outside of the boundaries')
%         j
%         i
%         u
%         x_next
%         return
%         break
%     end
%     
% end
% 
% 
% 
% plot(Xs(:,1), Xs(:,2), '-b')
% 


% axis([-0.012 0.056 -0.1415 -0.134])
%%
dim_y = 1;
dim_x = 4;

disp(['************************************************ vertice number: ', num2str(i)])

x_next = coeff*init_point';
hold on

subplot(2,2,2)

% h = figure('position',[100 100 200 200])
% h = figure('position',[300 100 200 200])
% h = figure('position',[500 100 200 200])
% h = figure('position',[700 100 200 200])

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')



hold on
u= [];

Xs = []

for j=0:300
    
    Xs = [Xs; [x_next(dim_x) x_next(dim_y)]]
    if (j == 0 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',10)
    elseif ( j < 50 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',1)
    elseif ( j == 50)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',10)
    elseif ( j < 300)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',1)
    else
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',10)
    end
    x_curr = x_next;
    
    if (j == 0)
        hu=ones(2*m,1)
        
        H1=[Hx*B_t1; Hx*B_t2; Hu];
        h=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
        
        U = lcon2vert(H1,h)
        size_u = size(U);
        if (size_u(1) > 0)
            
            u_mean= mean(U);
            u = u_mean';
            %             if (j < STEPS_FOR_CONTROLLED_OUTCOME)
            %                 j;
            %                 x_curr;
            %                 u = [];
            %                 u_mid = (max(U) + min(U))/2;
            %                 x_temp = A_t1*x_curr + B_t1*u_mid';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('mid')
            %                     u = u_mid';
            %                 end
            %
            %                 u_max = max(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %                 if size(U_temp,1) > 0
            %                     display('max')
            %                     u = u_max';
            %                 end
            %
            %                 u_min = min(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('min')
            %                     u = u_min';
            %                 end
            %
            %             else
            %                 u = min(U)/2 + max(U)/2;
            %             end
            
        else
            disp('There is no contorl input!')
            j
            i
            u
            x_curr
            return
        end
    end
    
    %discrete simulation of the path
    %     x_next = x_curr + (A_c*x_curr + B_c*u)*DT_CONTINUOUS;
    x_next = A_sim*x_curr + B_sim*u;
    
    
    if (F*x_next < g ~= ones(size(F,1),1))
        display('System went outside of the boundaries')
        j
        i
        u
        x_next
        return
        break
    end
    
end



plot(Xs(:,1), Xs(:,2), '-b')



axis([-0.012 0.056 -0.1415 -0.134])
%%
dim_y = 2;
dim_x = 3;


disp(['************************************************ vertice number: ', num2str(i)])

x_next = coeff*init_point';
hold on

% h = figure('position',[100 100 200 200])
% h = figure('position',[300 100 200 200])
% h = figure('position',[500 100 200 200])
% h = figure('position',[700 100 200 200])


subplot(2,2,3)

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')



hold on
u= [];

Xs = []

for j=0:300
    
    Xs = [Xs; [x_next(dim_x) x_next(dim_y)]]
    
    if (j == 0 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',10)
    elseif ( j < 50 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',1)
    elseif ( j == 50)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',10)
    elseif ( j < 300)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',1)
    else
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',10)
    end
    x_curr = x_next;
    
    if (j == 0)
        hu=ones(2*m,1)
        
        H1=[Hx*B_t1; Hx*B_t2; Hu];
        h=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
        
        U = lcon2vert(H1,h)
        size_u = size(U);
        if (size_u(1) > 0)
            
            u_mean= mean(U);
            u = u_mean';
            %             if (j < STEPS_FOR_CONTROLLED_OUTCOME)
            %                 j;
            %                 x_curr;
            %                 u = [];
            %                 u_mid = (max(U) + min(U))/2;
            %                 x_temp = A_t1*x_curr + B_t1*u_mid';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('mid')
            %                     u = u_mid';
            %                 end
            %
            %                 u_max = max(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %                 if size(U_temp,1) > 0
            %                     display('max')
            %                     u = u_max';
            %                 end
            %
            %                 u_min = min(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('min')
            %                     u = u_min';
            %                 end
            %
            %             else
            %                 u = min(U)/2 + max(U)/2;
            %             end
            
        else
            disp('There is no contorl input!')
            j
            i
            u
            x_curr
            return
        end
    end
    
    %discrete simulation of the path
    %     x_next = x_curr + (A_c*x_curr + B_c*u)*DT_CONTINUOUS;
    x_next = A_sim*x_curr + B_sim*u;
    
    
    if (F*x_next < g ~= ones(size(F,1),1))
        display('System went outside of the boundaries')
        j
        i
        u
        x_next
        return
        break
    end
    
end



plot(Xs(:,1), Xs(:,2), '-b')



axis([-0.05 0.07 -0.001 0.01])

%%
dim_y = 2;
dim_x = 4;


disp(['************************************************ vertice number: ', num2str(i)])

x_next = coeff*init_point';
hold on


subplot(2,2,4)

% h = figure('position',[100 100 200 200])
% h = figure('position',[300 100 200 200])
% h = figure('position',[500 100 200 200])
% h = figure('position',[700 100 200 200])

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')



hold on
u= [];

Xs = []

for j=0:300
    
    Xs = [Xs; [x_next(dim_x) x_next(dim_y)]]
    
    if (j == 0 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',10)
    elseif ( j < 50 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',1)
    elseif ( j == 50)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',10)
    elseif ( j < 300)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',1)
    else
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',10)
    end
    x_curr = x_next;
    
    if (j == 0)
        hu=ones(2*m,1)
        
        H1=[Hx*B_t1; Hx*B_t2; Hu];
        h=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
        
        U = lcon2vert(H1,h)
        size_u = size(U);
        if (size_u(1) > 0)
            
            u_mean= mean(U);
            u = u_mean';
            %             if (j < STEPS_FOR_CONTROLLED_OUTCOME)
            %                 j;
            %                 x_curr;
            %                 u = [];
            %                 u_mid = (max(U) + min(U))/2;
            %                 x_temp = A_t1*x_curr + B_t1*u_mid';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('mid')
            %                     u = u_mid';
            %                 end
            %
            %                 u_max = max(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %                 if size(U_temp,1) > 0
            %                     display('max')
            %                     u = u_max';
            %                 end
            %
            %                 u_min = min(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('min')
            %                     u = u_min';
            %                 end
            %
            %             else
            %                 u = min(U)/2 + max(U)/2;
            %             end
            
        else
            disp('There is no contorl input!')
            j
            i
            u
            x_curr
            return
        end
    end
    
    %discrete simulation of the path
    %     x_next = x_curr + (A_c*x_curr + B_c*u)*DT_CONTINUOUS;
    x_next = A_sim*x_curr + B_sim*u;
    
    
    if (F*x_next < g ~= ones(size(F,1),1))
        display('System went outside of the boundaries')
        j
        i
        u
        x_next
        return
        break
    end
    
end



plot(Xs(:,1), Xs(:,2), '-b')




axis([-0.012 0.056 -0.001 0.01])
%%
dim_y = 3;
dim_x = 4;


disp(['************************************************ vertice number: ', num2str(i)])

x_next = coeff*init_point';
hold on

% h = figure('position',[100 100 200 200])
% h = figure('position',[300 100 200 200])
% h = figure('position',[500 100 200 200])
h = figure('position',[700 100 200 200])

projChi = Chi.projection([dim_x, dim_y]);
plot(projChi);
xlabel(dim_names(dim_x),'FontSize',13,'FontWeight','bold','interpreter','latex')
ylabel(dim_names(dim_y),'FontSize',13,'FontWeight','bold','interpreter','latex')



hold on
u= [];

Xs = []

for j=0:300
    
    Xs = [Xs; [x_next(dim_x) x_next(dim_y)]]
    
    if (j == 0 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',10)
    elseif ( j < 50 )
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',1)
    elseif ( j == 50)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0 1],...
            'MarkerSize',10)
    elseif ( j < 300)
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',1)
    else
        plot(x_next(dim_x),x_next(dim_y),'-mo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[.49 1 .63],...
            'MarkerSize',10)
    end
    x_curr = x_next;
    
    if (j == 0)
        hu=ones(2*m,1)
        
        H1=[Hx*B_t1; Hx*B_t2; Hu];
        h=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
        
        U = lcon2vert(H1,h)
        size_u = size(U);
        if (size_u(1) > 0)
            
            u_mean= mean(U);
            u = u_mean';
            %             if (j < STEPS_FOR_CONTROLLED_OUTCOME)
            %                 j;
            %                 x_curr;
            %                 u = [];
            %                 u_mid = (max(U) + min(U))/2;
            %                 x_temp = A_t1*x_curr + B_t1*u_mid';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('mid')
            %                     u = u_mid';
            %                 end
            %
            %                 u_max = max(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %                 if size(U_temp,1) > 0
            %                     display('max')
            %                     u = u_max';
            %                 end
            %
            %                 u_min = min(U);
            %                 x_temp = A_t1*x_curr + B_t1*u_max';
            %                 h_temp=[hx-Hx*A_t1*x_curr; hx-Hx*A_t2*x_curr; hu];
            %                 U_temp = lcon2vert(H1,h_temp);
            %
            %                 if size(U_temp,1) > 0
            %                     display('min')
            %                     u = u_min';
            %                 end
            %
            %             else
            %                 u = min(U)/2 + max(U)/2;
            %             end
            
        else
            disp('There is no contorl input!')
            j
            i
            u
            x_curr
            return
        end
    end
    
    %discrete simulation of the path
    %     x_next = x_curr + (A_c*x_curr + B_c*u)*DT_CONTINUOUS;
    x_next = A_sim*x_curr + B_sim*u;
    
    
    if (F*x_next < g ~= ones(size(F,1),1))
        display('System went outside of the boundaries')
        j
        i
        u
        x_next
        return
        break
    end
    
end



plot(Xs(:,1), Xs(:,2), '-b')




axis([-0.012 0.056 -0.05 0.07])
