clear ''all'' 
close all 
clc
% 1 function to change here l45
%% figure start
figure(1) 
title('Chemin vers le minimum')
hold on
%% k : 1
starting_point = [4;-8];
image(1,1)=f_2(starting_point');
dim = length(starting_point);
alphaini = 1;
c1 =0.0001;
c2 = 0.9;
p =1/2;
alphamax=10000;
QNM = 0;
%%%% calcule du gradient %%%%%%%%%%%
x_sym = sym(zeros(1, dim));
for k=1:dim
    x_sym(k) = sym(sprintf('x%d', k));
end


F_1 = f_2(x_sym);
FX = gradient(F_1);
%%% calcule du gradient au point x_0 %%%
grad = zeros((dim),1);

for i=1:dim
    grad(i,1) = (subs(FX(i,1), x_sym, starting_point'));
end

%%% approche du hessien %%%%
B_0 = eye(length(starting_point),length(starting_point));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (1)
supx =abs(starting_point(1,1));
supy =abs(starting_point(2,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(-supx:0.5:supx,-supy:0.5:supy);
F = log(1 + 3.*(Y - (X.^3 - X)).^2 + (X - 4/3).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
surf(X,Y,F)
plot3(starting_point(1,1),starting_point(2,1),image(1,1),'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
text(starting_point(1,1),starting_point(2,1)+1,image(1,1)+1,'0','Color','red','FontSize',6)
drawnow
%%%%%%%%%%%%% FIN GRAPH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIND ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

direc = B_0*grad;
syms alpha_sym
phi_var = starting_point-alpha_sym*direc;
phi =f_2(phi_var');

alpha = linesearch1(phi,alpha_sym,alphaini,alphamax,c1,c2,QNM)

%% new point
newcord = starting_point- alpha*(B_0*grad)
t=2;
image(t,1) = f_2(newcord');
plot3(newcord(1,1),newcord(2,1),image(2,1),'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
text(newcord(1,1),newcord(2,1)+1,image(2,1),'1','Color','red','FontSize',2)
drawnow
QNM = QNM+1;
%% FIND B_1
delta= (newcord-starting_point);
grad2 = zeros((dim),1);
for i=1:dim
    grad2(i,1) = subs(FX(i,1), x_sym, newcord');
end
gamma = (grad2-grad);
 B_2= (B_0-((delta*gamma')/(gamma'*delta)))*B_0*(B_0-((gamma*delta')/(gamma'*delta)))+((delta*delta')/(gamma'*delta)); 
 B_K=B_2;

P = -(B_2*grad2);
eps = norm(grad2);
%% algorithme
while eps > 0.001
    t=t+1;
    starting_point= newcord;
    grad = zeros((dim),1);
    
    for i=1:dim
        grad(i,1) = subs(FX(i,1), x_sym, starting_point');
    end
    
    eps = norm(grad);
    if eps<0.001
        break
    end
    
    B_k=B_K;
    direc_k = B_k*grad;

    syms alpha_sym
    phi_var = starting_point-alpha_sym*direc_k;
    phi =f_2(phi_var');

    alpha = linesearch1(phi,alpha_sym,alphaini,alphamax,c1,c2,QNM)  
    
    %% new point
    newcord = starting_point- alpha*(B_k*grad)
    image(t,1) = f_2(newcord');
    plot3(newcord(1,1),newcord(2,1),image(t,1),'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
    %text(newcord(1,1),newcord(2,1),image(t,1),num2str(QNM+1),'Color','red','FontSize',10)
    drawnow
    %% FIND B_k+1
    delta=(newcord-starting_point);
    grad2 = zeros((dim),1);
    for i=1:dim
        grad2(i,1) = subs(FX(i,1), x_sym, newcord');
    end
    gamma =(grad2-grad);
    
    %B_K = B_k+(gamma*gamma'/(gamma'*delta))-((B_k*delta)*(B_k*delta)'/(delta'*B_k*delta));
    B_K= (B_0-((delta*gamma')/(gamma'*delta)))*B_k*(B_0-((gamma*delta')/(gamma'*delta)))+((delta*delta')/(gamma'*delta));
    QNM = QNM+1
end
hold off
%%% Probleme dans l update de Bk car elle ne change pas bien la direction