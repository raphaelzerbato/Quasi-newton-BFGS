function alpha_etoile= linesearch1(varargin)
phi = varargin{1};
alpha_sym = varargin{2};
alphain = varargin{3}; 
alphamax=  varargin{4};
c1 = varargin{5};
c2 = varargin{6};
QNM = varargin{7};
bracketing1 = 0;
%% calcule du gradient de phi et passage interval [0,1] %%%
dphi = diff(phi,alpha_sym);

phi_0 =  double(subs(phi,alpha_sym,0));
dphi_0 = double(subs(dphi,alpha_sym,0));

phi_1 = double(subs(phi,alpha_sym,alphain));
dphi_1= double(subs(dphi,alpha_sym,alphain));

if phi_1 > phi_0+c1*alphain*dphi_0 || QNM>0 && phi_1>= phi_0
    bracketing1 = 1;
    alpha_etoile = Zoom(phi,alpha_sym,c1,c2,0,alphain,QNM);
end
if bracketing1 == 0
    if abs(dphi_1)<=abs(c2*dphi_0)%-c2*dphi_0 
        bracketing1 = 2;
        alpha_etoile= alphain;
    end
end
if bracketing1 == 0
    if dphi_1>= 0
        bracketing1 = 3;
        alpha_etoile= Zoom(phi,alpha_sym,c1,c2,0,alphain,QNM);
    end
end
%% line search autre boucle
if bracketing1==0
    %bracketing2 = 0
    
    alphaold = alphain;
    alphanew =  alphaold*10;
    
    while alphanew<alphamax
        %incrementation des alpha
        phi_old =  double(subs(phi,alpha_sym,alphaold));
        dphi_old = double(subs(dphi,alpha_sym,alphaold));
        
        phi_new = double(subs(phi,alpha_sym,alphanew));
        dphi_new= double(subs(dphi,alpha_sym,alphanew));
        
        if phi_new > phi_old+c1*alphanew*dphi_old||QNM>0 && phi_new>= phi_old
            bracketing2 = 1;
            alpha_etoile = Zoom(phi,alpha_sym,c1,c2,alphaold,alphanew,QNM);
            return
        end
        if abs(dphi_new)<= abs(c2*dphi_old)%-c2*dphi_lo 
            bracketing2 = 2;
            alpha_etoile= alphanew;
            return
        end
        if dphi_new >= 0
            bracketing2 = 3;
            alpha_etoile= Zoom(phi,alpha_sym,c1,c2,alphaold,alphanew,QNM);
            return
        end
        
        alphaold =alphanew;
        alphanew= alphaold*10;
%         a = phi_old;
%         b = dphi_old;
%         dalpha = alphanew-alphaold;
%         c = (phi_new - phi_old - dalpha*dphi_old)/dalpha^2;
%         alphanew = alphaold - 0.5*b/c
    end
end
end
function alpha_etoile = Zoom(varargin)
phi = varargin{1};
alpha_sym = varargin{2}; 
c1 = varargin{3};
c2 = varargin{4};
alphalo= varargin{5};
alphahi = varargin{6};
QNM= varargin{7};
dphi = diff(phi,alpha_sym);
iter_cnt = 0;
zoom = 0;

while 1
    % calcule des phi
    phi_lo =  double(subs(phi,alpha_sym,alphalo));
    dphi_lo = double(subs(dphi,alpha_sym,alphalo));
    
    phi_hi = double(subs(phi,alpha_sym,alphahi));
    dphi_hi= double(subs(dphi,alpha_sym,alphahi));
    
    %quadratic interpolation
    a = phi_lo;
    b = dphi_lo;
    dalpha = alphahi-alphalo;
    c = (phi_hi - phi_lo - dalpha*dphi_lo)/dalpha^2;
    % iter_cnt
    if ( ( c <= 0 ) || ( mod(iter_cnt,3) == 2 ) )
        % Use bisection
        alpha1 = alphalo + 0.5*dalpha
       b=1;
    else
        % Use min of quadratic
        alpha1 = alphalo - 0.5*b/c
        q = 1;
    end
    
    %calcule de phi a alpha_1
    phi_alpha1= double(subs(phi,alpha_sym,alpha1));
    dphi_alpha1 = double(subs(dphi,alpha_sym,alpha1));
    
    if phi_alpha1 > phi_lo+c1*alpha1*dphi_lo || phi_alpha1 >= phi_lo
        zoom = 1;
        alphahi = alpha1;        
    elseif abs(dphi_alpha1 ) <= abs(c2*dphi_lo)%-c2*dphi_lo
        zoom = 2;
        alpha_etoile = alpha1;
        break
    elseif dphi_alpha1*(alphahi-alphalo)>=0
        zoom= 3;
        alphahi = alphalo;
        alphalo = alpha1;
        
    end
    if zoom == 0
        alphalo = alpha1;
        zoom=4;
    end
 iter_cnt = iter_cnt+1;   
 zoom = 0;
end
end
