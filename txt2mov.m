clear all;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DG1D Movie
% load ('data.txt');
% N = 20;
% time = 0.1;
% Nt = 1000;
% 
% count=1;
% dt = time/Nt;
% x = linspace(-1,1,N+1);
% figure
% 
% for i = 1:size(data,1)
%     plot(x,data(i,:),'k')
%     axis([-1,1,-1,1])
%     set(gca,'fontname','Times New Roman','fontsize',10,'fontweight','bold')
%     S = sprintf('Time=%.3f units',dt*i);
%     text(0.6,0.9,S,'fontname','Times New Roman','fontweight','bold')
%     M(count)=getframe;
%     count=count+1;
% end
% 
% % movie(M,1,1000);
% movie2avi(M, 'wave.avi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % DG2D Movie
% load ('data.txt');
% 
% N = 10;
% M = 5;
% kmax = size(data,1)/(N+1);
% 
% for k = 1:kmax
%     for i = 1:N+1
%         for j = 1:M+1
%             U(i,j,k) = data(i+(k-1)*(N+1),j);
%         end
%     end
% end
% 
% figure
% 
% x = 1:1:N+1;
% y = 1:1:M+1;
% [X,Y] = meshgrid(y,x);
% 
% count = 1;
% for k = 1:50
%     surf(U(:,:,k)','edgecolor','none')
%     % view(45, 30);
%     % contourf(Y,X,U(:,:,k),10,'LineStyle','none')
%     xlim([1 N])
%     ylim([1 M])
%     zlim([-0.2 0.2])
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     % colormap(bone)
%     set(gca,'fontname','Times New Roman','fontsize',10,'fontweight','bold')
%     S = sprintf('Iterations=%d',k);
%     text(0.6,1.9,S,'fontname','Times New Roman','fontweight','bold')
%     colorbar
%     caxis([-0.2 0.2])
%     W(count)=getframe;
%     count = count + 1;
% end
%  movie(W,1,kmax);
%  movie2avi(W, 'wave.avi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DGSEM1D Movie
load ('data.txt');

K = 5;
N = 14;
Npts = K*N;
Nt = 30+1;
x = reshape(data(:,1),[Npts,Nt]);
y = reshape(data(:,2),[Npts,Nt]);

time = 3.0;
dt = time/(Nt-1);
% x = linspace(0,K-1/N,Npts);
figure
count = 1;

for i = 1:size(y,2)
    plot(x(:,i),y(:,i),'k')
    axis([0,K,-1,2])
    %axis([0,K,-1,4])
    grid on;
    set(gca,'fontname','Times New Roman','fontsize',10,'fontweight','bold')
    S = sprintf('Time=%.1f units',(i-1)*dt);
    text(3.5,1.7,S,'fontname','Times New Roman','fontweight','bold')
    M(count)=getframe(gcf);
    count=count+1;
end

movie2avi(M, 'wave_DGSEM1D.avi');
