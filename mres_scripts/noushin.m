clear;

 

path=input('what is the path of the folder= ','s');

title_figure=input('what is the structure of the device= ','s');

files=dir(strcat(path,'/','*.csv'));

cell_files = cell(size(files));

 

for index=1:length(files)

    

data=readmatrix(strcat(files(index).folder,'/',files(index).name));

    cell_files{index}.struct=title_figure;

    cell_files{index}.folder=files(index).folder;

    cell_files{index}.name=files(index).name;

    cell_files{index}.iter=index;

    cell_files{index}.r=data(:,1);

    cell_files{index}.v=data(:,3);

    cell_files{index}.I=data(:,4);

    cell_files{index}.t=data(:,5);

    cell_files{index}.D=GetElectrodeDiameter(strcat(files(index).folder,'/',files(index).name));

    plot_IV(cell_files{index})

    %plot_IV(cell_files{index}.v,cell_files{index}.I,cell_files{index}.t,cell_files{index}.r,cell_files{index}.D)

 

end

 %%

function c = GetElectrodeDiameter(namefile)

a=extractBefore(namefile,'mm');

b=a((length(a)-3):end);

b(2)='.';

c=str2double(b);

end

 

function plot_IV(file)

 

 

counter=0;

for i=1:length(file.r)

    if ( file.r(i)==1)

        counter=i;

    else

        break;

    end

end

 

 

v=file.v(1:counter);

I=file.I(1:counter);

A= 10^-2* pi*(file.D/2)^2;

 

%it finds the index at which the voltage is maximum in the first cycle.

for i=1:length(v)

    if ( v(i)>v(i+1)&& v(i+1)>v(i+2)&& v(i+2)>v(i+3) )

        index_max=i ;

        break;   

    end

    if ((i+3)>=length(v))

        file.name

    end

        

end

 

%it finds the index at which the voltage passes zero in the first cycleis.

for i=index_max : length(v)

    if ( v(i)<0 && v(i+1)<0 && v(i+2)<0)

        index_zero=i;

        break

    end

end

 

file.name

 

%it finds the index at which the voltage is minimum in the first cycle.

for i=index_zero : length(v)

    if ( v(i)<v(i+1)&& v(i+1)<v(i+2)&& v(i+2)<v(i+3))

        index_min=i;

        break

    end

end

 

%it finds the index at which the voltage passes zero for second time in the first cycleis.

% for i=index_min : length(v)

%     if ( v(i)>0 && v(i+1)>0 && v(i+2)>0)

%         index_zero2=i;

%         break

%     end

% end

 

 

% catagorize data with the aid of the indexes

v1=v(1:index_max);

 

v2=v(index_max:index_zero);

 

v3=v(index_zero:index_min);

 

v4=v(index_min:length(v));

 

I1=I(1:index_max) ;

I2= I(index_max:index_zero) ;

I3=I(index_zero:index_min) ;

I4=I(index_min:length(v)) ;

 

 

 

 

% I1= 10^3* I(1:index_max)/A ;

% I2= 10^3*I(index_max:index_zero)/A ;

% I3=10^3* I(index_zero:index_min)/A ;

% I4=10^3* I(index_min:length(v))/A ;

 

area=string(round(A*10^2,3));

Area=append('elecrtode area = ',area,'mm^2');

 

scan_rate=string(round((v(1)-v(5))/(file.t(1)-file.t(5)),2));

scan=append('scan rate = ',scan_rate,' V/s');

 

%create plot of the current density versus voltage of different sequences

figure(1)

plot(v1,I1,'-o','LineWidth',2.1,'MarkerSize',0.3)

hold on

plot(v2,I2,'--o','LineWidth',2.1,'MarkerSize',0.3)

hold on

plot(v3,I3,'--o','LineWidth',2.1,'MarkerSize',0.3)

hold on

plot(v4,I4,'-o','LineWidth',2.1,'MarkerSize',0.3)

hold off

 legend({'1','2','3','4'}, 'Location','southeast','FontSize',15)

set(legend,'Box','off')

set(gca,'FontSize',25.5);

 

xlim([v(index_min)*1.1, v(index_max)*1.1])

% ylim([-5500,5500])

title(file.struct)

xlabel('Voltage (V)','FontSize',21)

ylabel('Current (A)','FontSize',21)

 text(0.08,0.88,Area,'Units','normalized','Color','black','FontSize',16)

 text(0.08,0.94,scan,'Units','normalized','Color','black','FontSize',16)

set(gca,'FontSize',27.5);

saveas(gcf, strcat(file.folder,'\JV\',extractBefore(file.name,".csv"),".png"))

 

%create plot of the logarithm of the current density versus voltage of different sequences

 

figure(2)

semilogy(v1,abs(I1),'-o','LineWidth',2.1,'MarkerSize',0.3)

hold on

semilogy(v2,abs(I2),'--o','LineWidth',2.1,'MarkerSize',0.3)

hold on

semilogy(v3,abs(I3),'--o','LineWidth',2.1,'MarkerSize',0.3)

hold on

semilogy(v4,abs(I4),'-o','LineWidth',2.1,'MarkerSize',0.3)

hold off

     legend({'1','2','3','4'}, 'Location','southeast')

     set(legend,'Box','off')

% xlim([v(index_min)*1.1, v(index_max)*1.1])

 

title(file.struct)

xlabel('Voltage (V)')

ylabel('Current (A)')

set(gca,'FontSize',27);

 

% text(0.08,0.84,scan,'Units','normalized','Color','black','FontSize',16)

% text(0.08,0.92,Area,'Units','normalized','Color','black','FontSize',16)

saveas(gcf, strcat(file.folder,'\JV\',extractBefore(file.name,".csv"),'semilogy',".png"))

end