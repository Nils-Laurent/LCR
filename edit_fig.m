figs_A = {'lin_CR568_SNRs', 'lin_CR1k2_SNRs', 'lin_CR3k5_SNRs', 'cosine_SNRs',...
    'lin_CR568_SNRs_cas2', 'lin_CR1k2_SNRs_cas2', 'lin_CR3k5_SNRs_cas2', 'cosine_SNRs_cas2',...
    'lin_2r_SNRs'};

for fig_name=figs_A
    fpath = "./fig/"+fig_name;
    fprintf(fpath+"\n");
    
    fig = openfig(fpath);
    f_lgd = findobj('type','legend');
    f_lgd.FontSize = 34;
    f_lgd.Location ='northwest';
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 30);
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 30);
    pbaspect([1 1 1]);
    set(gcf, 'PaperPosition', [0 0 25 25]);
    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
      lines(i).LineWidth = 2.0;
      lines(i).MarkerSize = 10.0;
      if lines(i).LineStyle == ':'
          lines(i).LineStyle = '-.';
      end
    end
    savefig(fpath)
    saveas(fig,fpath,'epsc')
    close all;
end

figs_B = {'cosine_TFR', 'lin_TFR_all', 'lin_2r_TFR'};
for fig_name=figs_B
    fpath = "./fig/"+fig_name;
    fprintf(fpath+"\n");
    
    fig = openfig(fpath);
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 30);
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 30);
    pbaspect([1 1 1]);
    set(gcf, 'PaperPosition', [0 0 25 25]);
    savefig(fpath)
    saveas(fig,fpath,'epsc')
    close all;
end

figs_C = {'lin_2r_diff_l2'};
for fig_name=figs_C
    fpath = "./fig/"+fig_name;
    fprintf(fpath+"\n");
    
    fig = openfig(fpath);
    f_lgd = findobj('type','legend');
    f_lgd.FontSize = 32;
    f_lgd.Location ='northeast';
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 30);
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 30);
    pbaspect([1 1 1]);
    set(gcf, 'PaperPosition', [0 0 25 25]);
    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
      lines(i).LineWidth = 2.0;
      lines(i).MarkerSize = 10.0;
      if lines(i).LineStyle == ':'
          lines(i).LineStyle = '-.';
      end
    end
    savefig(fpath)
    saveas(fig,fpath,'epsc')
    close all;
end

figs_C = {'estimators'};
for fig_name=figs_C
    fpath = "./fig/"+fig_name;
    fprintf(fpath+"\n");
    
    fig = openfig(fpath);
    f_lgd = findobj('type','legend');
    f_lgd.FontSize = 28;
    f_lgd.Location ='northwest';
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 30);
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 30);
    %pbaspect([1 1 1]);
    %set(gcf, 'PaperPosition', [0 0 25 25]);
    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
      lines(i).LineWidth = 2.0;
      lines(i).MarkerSize = 10.0;
    end
    savefig(fpath)
    saveas(fig,fpath,'epsc')
    close all;
end
