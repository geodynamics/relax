function [sfl]=plot_faults(faults,dim,xlim,ylim);
%fid=fopen('saf-pkfd.flt','wt');
[sfl,l]=size(dim);
for i=1:sfl
    if i==1
        S1=1;
    else
        S1=S2+1;
    end
    S2=S1+dim(i)-1;
    ab=faults(S1:S2,1);
    or=faults(S1:S2,2);
    if min(ab)>xlim(1) & max(ab)<xlim(2) & min(or)>ylim(1) & max(or)<ylim(2)
        line(ab,or,'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',0.5); hold on
        %fprintf(fid,'> \n');
        %     for j=1:length(ab)
        %       fprintf(fid,'%9.4f %8.4f \n',ab(j),or(j));
        %     end
    end
end
%fclose(fid);

