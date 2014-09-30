

load('figs/simulation_data');
nSteps = size(calcs_set,1)/2;
ndomains = 2;
nb = 2;
ndof = 7;


ll = 1;
fig_hl = figure(1000);
set(fig_hl,'Position', [20 80 800 800])
    clf;
    plot([-1.5 2.2], [0 0], 'b:');
    hold on;
    sleg = plot([0 0], 'r','LineWidth',2);
    nleg = plot([0 0], 'b','LineWidth',2);
    torso = plot([0 0], 'g','LineWidth',2);
    pZMP = plot([0 0], 'ro','LineWidth',2);    
    hold off;
    anim_axis=[-0.6 2.4 -0.1 0.5];
    axis off
    axis(anim_axis)
    axis equal
    grid
    
    tic;
    calcs_full = [calcs_set{:}];
    out = horzcat_fields_domains(calcs_full);
    px=0;
    pxtemp=0;
     for k = 1:nsteps
%          for j = 1:ndomains
            calcs = calcs_set{k};
            ts = calcs.t;
            ts = ts - ts(1);
            qs = calcs.q_anim;
            oZMP = calcs.ZMP;
             [Et,Ex] = even_sample(ts,qs,10);
             [Et,ZMP] = even_sample(ts,oZMP,10);             
             ncalcs = length(Et);
             for i = 1:ncalcs%5:428
                % position
%                  qe = calcs.qe(1:9, i);
                 qe = Ex(1:end,i);
                pos = jpos_mat(qe)+[ones(1,13);zeros(1,13)]*px;
                
                posStack(:,i)=[pos(1,:)';pos(2,:)'];
                set(sleg, 'XData', pos(1,1:8), 'YData', pos(2,1:8), ...
                    'erasemode', 'normal');
                set(nleg, 'XData', pos(1,8:end), 'YData', pos(2,8:end), ...
                    'erasemode', 'normal');
                set(torso, 'XData', pos(1,6:7), 'YData', pos(2,6:7), ...
                    'erasemode', 'normal');
                set(pZMP, 'XData', ZMP(i)+px, 'YData', 0, ...
                    'erasemode', 'normal');                
                
%                new_axis=anim_axis+[pos(1,6) pos(1,6) 0 0];
%                 axis(new_axis)
    
    
                drawnow();
                pause(0.01);
                if 1
                    % Regulate frameIndex rate, use timeFactor to slow down / speed up
                    pauseTime = Et(i) - toc();
                    if pauseTime >= 0
                        pause(pauseTime);
                        % Shouldn't be non-zero, except in maybe rare slowdowns
                    end
                    % Skip frames otherwise?
                else
                    % Mini pause for funsies
                    pause(0.00001);
                end
                M(ll)= getframe();
                frame = getframe();
                if ll == 1
                    % Setup - not sure if there's a better way to actually have
                    % it work
                    [gifImg, gifMap] = rgb2ind(frame.cdata, 255, 'nodither');
                    % Set length... somehow
                    % Just use small increment since we're dynamically sizing
                    % Gonna be inefficient, but oh wells
                    gifImg(1, 1, 1, 2) = 0;
                else
                    gifImg(:, :, 1, ll) = rgb2ind(M(ll).cdata, gifMap, 'nodither');
                end
                ll = ll+1;
                
                
             end
             if mod(k,2)==1
                pxtemp = pos(1,13)-pos(1,1);
             else
                px= pxtemp;
             end
             
          end
        % Update coordinates
%         qe = calcs.qe(:, end);
%         q = qe(nb + (1:ndof));
%         qb = pe_nsf_mat(q, qb);
%      end
    
%     movie2avi(M,'proxi_slip','fps',60,'compression','None','quality',100);
    imwrite(gifImg, gifMap, ['figs/proxi_slip.gif'], 'DelayTime', 1/60, 'LoopCount', Inf);