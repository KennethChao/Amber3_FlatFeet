function [a_opt,vhip,p_opt] = find_transition(a0,xb_plus,xb_minus,cur,Jcur)




% a_opt = zeros(1,7);
alpha = zeros(4,7);


x_plus = xb_plus([3:11 14:end]);
x_minus = xb_minus([3:11 14:end]);


y0 = cur(x_plus,xb_plus);
dy0 = Jcur(x_plus,xb_plus)*[xb_plus(12:end);zeros(11,1)];

yf = cur(x_minus,xb_minus);
dyf = Jcur(x_minus,xb_minus)*[xb_minus(12:end);zeros(11,1)];


p0 = pe_com_mat(x_plus,xb_plus);
pf = pe_com_mat(x_minus,xb_minus);
p_opt = [p0(1);pf(1)];

vhip = (p_opt(2)-p_opt(1));
dtdp = 1/vhip;
Jhip0 = Je_com_mat(x_plus,xb_plus);
Jhipf = Je_com_mat(x_minus,xb_minus);


    
    alpha(:,[2 4 6]) = ones(4,1) * a0;
    alpha(:,[1 3 5 7]) = eye(4);
    Y = [yf;dyf;y0;dy0];
    N = 4;
    yd0 = zeros(N,1);
    ydf = zeros(N,1);
    dyd0 = zeros(N,1);
    dydf = zeros(N,1);
    for i = 1:N
        yd0(i) = y_ExtCwf(0, alpha(i,:));
        ydf(i) = y_ExtCwf(1, alpha(i,:));
    
        dydt0 = dy_ExtCwf(0, alpha(i,:));
        dyd0(i) = dydt0 * (Jhip0(1,:)*dtdp) * xb_plus(12:end);
    
        dydtf = dy_ExtCwf(1, alpha(i,:));
        dydf(i) = dydtf * (Jhipf(1,:)*dtdp) * xb_minus(12:end);
    end
    h_alpha_rem = [ydf';dydf';yd0';dyd0'] \ Y;
    a_opt = [h_alpha_rem(1),a0(1),h_alpha_rem(2),...
        a0(2),h_alpha_rem(3),a0(3),h_alpha_rem(4)];
    

    yd = zeros(1,51);
    for k = 1:51
        t = (k-1)/50;
        yd(k) = y_ExtCwf(t, a_opt);
    end
    plot(yd)
end