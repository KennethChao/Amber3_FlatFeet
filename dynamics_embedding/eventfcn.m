function [value, isterminal, direction] = eventfcn(t, x, domain,ref)
direction = -1;
isterminal = 1;

switch domain.type
    case 1
        %         calc = ref.h.calc;
        %         value = calc.Fv_ns - 1e-5;
%         q_nsf = y_NonStanceFoot(x);
%         value = 0.6 + q_nsf;
%         disp(['domain 1:',num2str(value)]);
        %         q_sf = y_StanceFoot(x);
%         if q_sf < 0 % ignore early scuffing
            value = h_nsheel_sca(x);
%         else
%             value = 1;
%         end
        disp(['domain 1:',t]);
        
    case 2
%         q_sf = y_StanceFoot(x);
%         if q_sf < 0 % ignore early scuffing
            value = h_nsheel_sca(x);
%         else
%             value = 1;
%         end
% num2str(value)
        disp(['domain 2:',t]);
    case 3
        value = h_nstoe_sca(x);
        disp(['domain 3:',num2str(value)]);
    otherwise
        error('Unknown Domain Type!');
end


end