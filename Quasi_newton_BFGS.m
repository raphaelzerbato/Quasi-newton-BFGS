clear ''all'' 
close all 
clc
% 4 Fonctions à changer ici (l22, l56, l81, l122)+ fonction f_2 a changer
%% figure start
figure(1) 
title('Chemin vers le minimum')
hold on
%% k : 1
x = [4;-8];
image(1,1)=f_2(x');
dim = length(x);
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
    grad(i,1) = (subs(FX(i,1), x_sym, x'));
end

%%% approche du hessien %%%%
B_0 = eye(length(x),length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (1)
supx =abs(x(1,1));
supy =abs(x(2,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(-supx:0.5:supx,-supy:0.5:supy);
F = log(1 + 3.*(Y - (X.^3 - X)).^2 + (X - 4/3).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
surf(X,Y,F)
plot3(x(1,1),x(2,1),image(1,1),'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
text(x(1,1),x(2,1)+1,image(1,1)+1,'0','Color','red','FontSize',6)
drawnow
%%%%%%%%%%%%% FIN GRAPH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIND ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

direc = B_0*grad;
syms alpha_sym
phi_var = x-alpha_sym*direc;
phi =f_2(phi_var');

alpha = linesearch1(phi,alpha_sym,alphaini,alphamax,c1,c2,QNM)

%% new point
newcord = x- alpha*(B_0*grad)
t=2;
image(t,1) = f_2(newcord');
plot3(newcord(1,1),newcord(2,1),image(2,1),'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
text(newcord(1,1),newcord(2,1)+1,image(2,1),'1','Color','red','FontSize',2)
drawnow
QNM = QNM+1;
%% FIND B_1
delta= (newcord-x);
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
    x= newcord;
    grad = zeros((dim),1);
    
    for i=1:dim
        grad(i,1) = subs(FX(i,1), x_sym, x');
    end
    
    eps = norm(grad);
    if eps<0.001
        break
    end
    
    B_k=B_K;
    direc_k = B_k*grad;

    syms alpha_sym
    phi_var = x-alpha_sym*direc_k;
    phi =f_2(phi_var');

    alpha = linesearch1(phi,alpha_sym,alphaini,alphamax,c1,c2,QNM)  
    
    %% new point
    newcord = x- alpha*(B_k*grad)
    image(t,1) = f_2(newcord');
    plot3(newcord(1,1),newcord(2,1),image(t,1),'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
    %text(newcord(1,1),newcord(2,1),image(t,1),num2str(QNM+1),'Color','red','FontSize',10)
    drawnow
    %% FIND B_k+1
    delta=(newcord-x);
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