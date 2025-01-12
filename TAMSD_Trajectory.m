function [ ListeTAMSD ] = TAMSD_Trajectory( Trajectory,lagtime_max )
%TAMSD pour une trajectoire
NombrePoint=size(Trajectory,1);
M=size(Trajectory,2);
Resultat=zeros(1,2);
ListeTAMSD=zeros(lagtime_max,M);

for TailleFenetre=1:lagtime_max
    TailleListeMoyenne=NombrePoint-TailleFenetre;
%     Somme=0;
%     for i=1:TailleListeMoyenne
%         Somme=Somme+(Trajectory(i+TailleFenetre)-Trajectory(i))^2;
%     end;
Somme=sum((Trajectory((1:TailleListeMoyenne)+TailleFenetre,:)-Trajectory((1:TailleListeMoyenne),:)).^2,1,'omitnan');
    R=Somme./(TailleListeMoyenne);
    ListeTAMSD(TailleFenetre,:)=R;
end
%loglog(ListeTAMSD)
% xlabel('\Delta')
% ylabel('\delta^2(\Delta)')
% fit=polyfit(log(1:lagtime_max),log(ListeTAMSD(:,1)'),1);
% Resultat(1,1)=fit(1,1); %exposant alpha
% Resultat(1,2)=exp(fit(1,2));% D


