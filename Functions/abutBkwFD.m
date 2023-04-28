function [bilFDunit,refFDunit,kcm,press]=abutBkwFD(H,type,plotter)
%calculate force-displacement relationship of abutment-backwall systems
%according to shamsabadi 2007

%THE OUTPUT LAWS ARE RELATED TO A unit-width (1 METER-portion) OF WALL

if strcmp(type,'coh')
    %Cohesive
    maxdef=0.1;ymax=H*maxdef*100; %[cm]
    temp=0:0.01:ymax;
    Fvect=temp./(.00191+.00204.*temp); %[kN per meter of wall]
    
elseif strcmp(type,'gra')
    %Granular
    maxdef=0.05;ymax=H*maxdef*100; %[cm]
    temp=0:0.01:ymax;
    Fvect=temp./(.00383+.00204.*temp); %[kN per meter of wall]
end
    

refFDunit=[temp'/100,Fvect']%[kN,m] per meter of the wall

yave=interp1(Fvect,temp,max(Fvect)/2);
kcm=max(Fvect)/2/yave; %equivalent stiffness of the wall per meter (kN/cm/m)
press=max(Fvect)/H %pressure in kPa

km=kcm*100; %[k per meter of wall]
Fmax=max(Fvect); %kN
deltay=Fmax/km;
bilFDunit=[0 0; deltay Fmax; ymax/100 Fmax]

%output
if strcmp(plotter,'plot')
    figure
    hold on
    plot(refFDunit(:,1),refFDunit(:,2),'-k')
    plot(bilFDunit(:,1),bilFDunit(:,2),'--k')    
end

end

