%% Post Process LinAnalysis


% name states
st = linavg.x_desc;
for j = 1:length(st)
    if contains(st(j),'blade 1')
        st{j} = strrep(st{j},'blade 1', 'coning mode');
    elseif contains(st(j),'blade 2')
        st{j} = strrep(st{j},'blade 2', 'cosine-cyclic mode');
    elseif contains(st(j),'blade 3')
        st{j} = strrep(st{j},'blade 3', 'sine-cyclic mode');
    end
    if contains(st(j), 'First time derivative')
        st{j} = strrep(st{j},'First time derivative', 'd/dt');
    end
    
    sta = strfind(st{j},'(');
    if contains(st(j), 'm/s')
        en = length(st{j}) - 5;
    elseif contains(st(j), 'rad/s')
        en = length(st{j}) - 7;
    elseif contains(st(j), 'rad')
        en = length(st{j}) - 5;
    elseif contains(st(j), 'm')
        en = length(st{j}) - 3;
    end
    
    
    st{j}(sta-1:en) = [];
    st{j}(1:3) = [];
    st{j} = strcat(st{j},' - ',num2str(j));
end
% for j = 1:length(linout)
    

%% some plots
st2 = flipud(st);
cats = categorical(st,st2);
ctrmag = fliplr(abs(linavg_t.SVD.U(1,:)));
ctrmag2 = fliplr(abs(linavg_b.SVD.U(1,:)));
ctrmag3 = fliplr(abs(linavg_tb.SVD.U(1,:)));
fig = figure(1)
b = barh(cats,[ctrmag; ctrmag2; ctrmag3]')
fig.Children.FontSize = 20;
leg = legend('\tau','\beta','\tau\beta')
leg.FontSize = 20;


figure(2)
myplot(1:length(linavg_tb.SVD.S),linavg_tb.SVD.S)
title('Singular Values for Collective Pitch + Torque Control')
%%
% st2 = flipud(st);
% cats = categorical(st,st2);
% mag = fliplr(abs(linavg.SVD.U(1,:)));
% fig = figure
% b = barh(cats,mag')
% % legend('$\tau$','$\beta$','$yaw$')
% fig.Children.FontSize = 20


