function patch_downstroke(h_fig, t, Euler_R_dot)

figure(h_fig);

index_downstroke=find(Euler_R_dot(1,:)>0);
index_upstroke=find(Euler_R_dot(1,:)<0);
%plot(t(index_downstroke), Euler_R(1,index_downstroke), 'b.');

index_down_end=[index_downstroke(find(diff(index_downstroke)>1)) index_downstroke(end)];
index_down_begin=[index_upstroke(find(diff(index_upstroke)>1)) index_upstroke(end)];
%plot(t(index_down_end), Euler_R(1,index_down_end), 'r.');
%plot(t(index_down_begin), Euler_R(1,index_down_begin), 'ro');

for i=1:length(index_down_begin)   
    patch_index=[index_down_begin(i), ...
        index_down_end(min(find(index_down_end>index_down_begin(i))))];
   
    if size(patch_index,2) > 1
        Ylim=ylim;
        patch(t([patch_index flip(patch_index)]), kron(ylim,[1 1]), 0.5*[1 1 1],'FaceAlpha',0.2,'LineStyle','none');
    end
end
hold on;

end

