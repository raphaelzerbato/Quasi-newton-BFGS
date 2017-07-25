clear ''all'' 
close all 
clc
% 1 Function to change here l45 
%% figure start
figure(1) 
title('Chemin vers le minimum')
hold on
%% k : 1
starting_point = [15;-20];
image_x1=f_2(starting_point);
dim = length(starting_point);
alphaini = 1;
c1 =0.0001;
c2 = 0.9;
p =1/2;
alphamax=1000;
QNM = 0;
%%%% calcule du gradient %%%%%%%%%%%
x_sym = sym(zeros(1, dim));
for k=1:dim
    x_sym(k) = sym(sprintf('x%d', k));
end
F_1 = f_2(x_sym')
FX = gradient(F_1);

%%% gradiant at starting point %%%
grad = zeros((dim),1);

for i=1:dim
    grad(i,1) = (subs(FX(i,1), x_sym, starting_point'));
end

%%% approche du hessien %%%%
B_0 =eye(length(starting_point),length(starting_point));

%%%%% Graf creation. To change graf aspect go to  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% see https://fr.mathworks.com/help/matlab/ref/meshgrid.html %%%%%%%%%
figure (1)
supx =abs(starting_point(1,1));
supy =abs(starting_point(2,1));

%%%%% FUNCTION TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(-supx-5:1:supx+5,-supy-200:10:supy+200);
F = (1-X).^2+100*(Y-X.^2).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
surf(X,Y,F)
plot3(starting_point(1,1),starting_point(2,1),image_x1,'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
text(starting_point(1,1),starting_point(2,1)+1,image_x1+1,'1','Color','red','FontSize',6)
drawnow
%%%%%%%%%%%%% FIN GRAPH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIND ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

direc = B_0*grad;
newcord=starting_point-B_0*grad;
image_grad = f_2(newcord);

syms alpha_sym
phi_var = starting_point-alpha_sym*direc;
phi =f_2(phi_var);

alpha = linesearch1(phi,alpha_sym,alphaini,alphamax,c1,c2,QNM)

%% new point
newcord = starting_point- alpha*(B_0*grad)
image_point2 = f_2(newcord);
plot3(newcord(1,1),newcord(2,1),image_point2,'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
text(newcord(1,1),newcord(2,1)+1,image_point2+1,'2','Color','red','FontSize',6)
drawnow
QNM = QNM+1;
%% FIND B_1
delta= (newcord-starting_point);
grad2 = zeros((dim),1);
for i=1:dim
    grad2(i,1) = subs(FX(i,1), x_sym, newcord');
end
gamma = (grad2-grad);
B_2 = B_0+ (((delta-B_0*gamma)*(delta-B_0*gamma)')/((delta-B_0*gamma)'*gamma));
B_K=B_2;

P = -(B_2*grad2);
eps = norm(grad2);
t=2
image(t,1) = f_2(newcord)

%% algorithme
while eps>0.001
    t = 1+t;
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
    phi =f_2(phi_var);
    alpha = linesearch1(phi,alpha_sym,alphaini,alphamax,c1,c2,QNM)  
    
    %% new point
    newcord = starting_point- alpha*(B_k*grad)
    image(t,1) = f_2(newcord);
    plot3(newcord(1,1),newcord(2,1),image(t,1),'or','MarkerSize',5,'MarkerEdgeColor','b',...
    'MarkerFaceColor','w')
    %text(newcord(1,1),newcord(2,1)+1,image(t,1)+1,num2str(QNM+2),'Color','red','FontSize',6)
    drawnow
    %% FIND B_k+1
    delta=(newcord-starting_point);
    grad2 = zeros((dim),1);
    for i=1:dim
        grad2(i,1) = subs(FX(i,1), x_sym, newcord');
    end
    gamma =(grad2-grad);
    B_K = B_k+ (((delta-B_k*gamma)*(delta-B_k*gamma)')/((delta-B_k*gamma)'*gamma));
    QNM = QNM+1;
end
hold off
