function [] = disp_2D(xf_2, Yt, xp, model, mean_t, n, Cycle_count, Acycles)
    figure(1);
    subplot(4,2,[1 3]);
    plot(xf_2(:,1),xf_2(:,2),'og','LineWidth',3,'MarkerSize',6);
    hold on;
    plot(Yt(1),Yt(2),'^b','LineWidth',3,'MarkerSize',10);
    hold off;
    title('prior');
    
    subplot(4,2,[5 7]);
    plot(xp(:,1),xp(:,2),'og','LineWidth',3,'MarkerSize',6);
    hold on;
    plot(Yt(1),Yt(2),'^b','LineWidth',3,'MarkerSize',10);
    plot(model.x_real(1),model.x_real(2),'pr','LineWidth',3,'MarkerSize',10);    
    plot(mean_t(1,n, Cycle_count),mean_t(2,n, Cycle_count),'Or','LineWidth',3,'MarkerSize',10);
    hold off;
    title('posterior');
    
    subplot(4,2,2);
    plot(model.x_real(1)*ones(1,Acycles(n)-Acycles(1)+1),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(1,1:n, Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 1');
    
    subplot(4,2,4);
    plot(model.x_real(2)*ones(1,Acycles(n)-Acycles(1)+1),'-*','LineWidth',2,'MarkerSize',6);
    hold on;
    plot(mean_t(2,1:n, Cycle_count),'-or','LineWidth',2,'MarkerSize',6);
    hold off;
    title('variable 2');  
end