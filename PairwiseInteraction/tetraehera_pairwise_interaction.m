function tetraehera_pairwise_interaction
%written by Jiahui Li 
clf;
%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%%%%%% 
ionic_strength=[20 22 25 30 50 60 70];%mM
%%%%%%%%%%%%constants%%%%%%%%%%%%%%%%%%%%%%%%%%
H=1.2*10^-19; %harmaker constant 10^-19J
kb=1.38*10^-23;
kbT=kb*298;
elementary_charge=1.6*10^-19;%coulombs
Na = 6.02*10^23;%Avogadro constant
e0=8.85*10^-12;%vaccum permittivity F/m
e_water=e0*78.4;%water permittivity F/m
viscosity = 0.8921*10^(-3);  %Pa s=Ns/m2 298K
zeta_potential = 51.5; %mv
lb=0.7; %nm
ligand_density = 1; %number/nm^2
zb = 0.0615;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
diameter = 0.667;
tetrasize = 1; %1: small 2: mid 3: large
switch tetrasize 
    case 3
        L = 97.2;
        LT = 14.3;
        LE = L-2*LT;
        tetradist = 14.7;
    case 2
        L = 64.51;
        LT = 10.68;
        LE = L-2*LT;
        tetradist = 9;
    case 1
        L = 45.74;
        LT = 6.82;
        LE = L-2*LT;
        tetradist = 8.1;
    otherwise
        disp('Not included')
end

filenameatom = ['obj_to_coordinate/perpendicular_' num2str(L) '_' num2str(LT) '_atom_coordinate.txt'];
filenameligand = ['obj_to_coordinate/perpendicular_' num2str(L) '_' num2str(LT) '_0.6tr_ligand_coordinate.txt'];
atoms = importdata(filenameatom);
ligands = importdata(filenameligand);

units1 = gpuArray(atoms);            %get the postions and the diameter of each divided spheres
ligands1 = gpuArray(ligands);
center1 = [0 0 0];
uc = 0;

disp('Tetrahedra loaded!');

xangle = 180;
yangle = 60;
zangle = 30;
zangle2 = 60;
zangle3 = 0; %rotation

Rx = [1 0 0
    0 cosd(xangle) -sind(xangle);
    0 sind(xangle) cosd(xangle)];

Ry = [cosd(yangle) 0 sind(yangle);
    0 1 0;
    -sind(yangle) 0 cosd(yangle)];

Rz = [cosd(zangle) -sind(zangle) 0;
     sind(zangle) cosd(zangle) 0;
     0 0 1];
 
Rz2 = [cosd(zangle2) -sind(zangle2) 0;
     sind(zangle2) cosd(zangle2) 0;
     0 0 1];
 
Rz3 = [cosd(zangle3) -sind(zangle3) 0;
     sind(zangle3) cosd(zangle3) 0;
     0 0 1];
%not moving tetra
units1 = units1*Rz;
ligands1 = ligands1*Rz;
%moving tetra in upper layer
center2_0 = center1;
units2_0 = units1*Rx;
ligands2_0 = ligands1*Rx;
units2_0 = units2_0*Rz2;
ligands2_0 = ligands2_0*Rz2;

units2_0(:,1) = units2_0(:,1)+(LE+LT)/(2*sqrt(3))+tetradist/3;
units2_0(:,2) = units2_0(:,2)+(LE+LT)/2+tetradist/sqrt(3);
units2_0(:,3) = units2_0(:,3)+L*sqrt(6)/6 -LT*sqrt(6)/3;
ligands2_0(:,1) = ligands2_0(:,1)+(LE+LT)/(2*sqrt(3))+tetradist/3;
ligands2_0(:,2) = ligands2_0(:,2)+(LE+LT)/2+tetradist/sqrt(3);
ligands2_0(:,3) = ligands2_0(:,3)+L*sqrt(6)/6 -LT*sqrt(6)/3;
center2_0(1) = center2_0(1)+(LE+LT)/(2*sqrt(3))+tetradist/3;
center2_0(2) = center2_0(2)+(LE+LT)/2+tetradist/sqrt(3);
center2_0(3) = center2_0(3)+L*sqrt(6)/6 -LT*sqrt(6)/3;

units1 = units1*Rz3;
ligands1 = ligands1*Rz3;

%moving tetra in lower layer
center3_0 = center1;
units3_0 = units1;
ligands3_0 = ligands1;
units3_0(:,1) = units3_0(:,1)+(LE+LT)*sqrt(3)/2+tetradist;
units3_0(:,2) = units3_0(:,2)+(LE+LT)/2+tetradist/sqrt(3);
ligands3_0(:,1) = ligands3_0(:,1)+(LE+LT)*sqrt(3)/2+tetradist;
ligands3_0(:,2) = ligands3_0(:,2)+(LE+LT)/2+tetradist/sqrt(3);
center3_0(:,1) = center3_0(:,1)+(LE+LT)*sqrt(3)/2+tetradist;
center3_0(:,2) = center3_0(:,2)+(LE+LT)/2+tetradist/sqrt(3);

vector = [-(LE+LT)/sqrt(3)-tetradist*2/3,0];
vector2 = [-(LE+LT)*sqrt(3)/2-tetradist,(LE+LT)/2+tetradist/sqrt(3)];

index = [0:0.01:0.6];


figure(1)
clf;
hold on
axis equal
xlim([-60 120])
ylim([-60 90])
zlim([-20 75])

scatter3(units1(:,1),units1(:,2),units1(:,3),1,'r')
scatter3(units2_0(:,1),units2_0(:,2),units2_0(:,3),1,'b')
scatter3(units3_0(:,1),units3_0(:,2),units3_0(:,3),1,'g')

scatter3(ligands1(:,1),ligands1(:,2),ligands1(:,3),1,[0.9 0.1 0.1])
scatter3(ligands2_0(:,1),ligands2_0(:,2),ligands2_0(:,3),1,[0.1 0.1 0.9])
scatter3(ligands3_0(:,1),ligands3_0(:,2),ligands3_0(:,3),1,[0.1 0.9 0.1])
%% 

E_sumss1=zeros(length(ionic_strength),length(index));
E_vdwss1=zeros(size(index));
E_elss1=zeros(length(ionic_strength),length(index));
E_sumss2=zeros(length(ionic_strength),length(index));
E_vdwss2=zeros(size(index));
E_elss2=zeros(length(ionic_strength),length(index));
E_totalss=zeros(length(ionic_strength),length(index));
ctcd=zeros([1 length(index)]);

overlapp = 0;
zzi = 1;
for ind = index
    cla;
    tic
    E_sum1= 0;
    E_sum2= 0;
    E_total= 0;

    disp(ind);
    units2=units2_0;
    ligands2 = ligands2_0;
    center2 = center2_0;
    units2(:,1) = units2(:,1)+vector(1)*ind;
    units2(:,2) = units2(:,2)+vector(2)*ind;
    ligands2(:,1) = ligands2(:,1)+vector(1)*ind;
    ligands2(:,2) = ligands2(:,2)+vector(2)*ind;
    center2(1) = center2(1)+vector(1)*ind;
    center2(2) = center2(2)+vector(2)*ind;

    units3=units3_0;
    ligands3 = ligands3_0;
    units3(:,1) = units3(:,1)+vector2(1)*ind;
    units3(:,2) = units3(:,2)+vector2(2)*ind;
    ligands3(:,1) = ligands3(:,1)+vector2(1)*ind;
    ligands3(:,2) = ligands3(:,2)+vector2(2)*ind;
    clf;
    hold on
    scatter3(units1(:,1),units1(:,2),units1(:,3),1,'r')
    scatter3(units2(:,1),units2(:,2),units2(:,3),1,'b')
    scatter3(units3(:,1),units3(:,2),units3(:,3),1,'g')
    scatter3(ligands1(:,1),ligands1(:,2),ligands1(:,3),1,'r','filled')
    scatter3(ligands2(:,1),ligands2(:,2),ligands2(:,3),1,'b','filled')
    scatter3(ligands3(:,1),ligands3(:,2),ligands3(:,3),1,'g','filled')
    xlabel("x")
    ylabel("y")
    zlabel("z")
    view(0,90)
    axis equal
    drawnow

    ctcd(1,zzi) = sqrt(sum((center2-center1).^2));


    %vdw
    gpu = gpuDevice;
    max_szie = gpu.AvailableMemory;

    max_number_double  =floor((max_szie/8)^0.5);

    max_chunk_size = floor(max_number_double*0.4);
    unit_chunk_number = ceil(length(units1)/max_chunk_size);
    disp(['Units dividing points into ',num2str(unit_chunk_number),' chunks...'])

    bins = 1 : max_chunk_size : length(units1);
    bins(end+1) = length(units1);
    ua_cell = {};
    ub1_cell = {};
    ub2_cell = {};

    for i = 1 : unit_chunk_number
        if i == unit_chunk_number
            ua_cell{i} = units1(bins(i):bins(i+1),:);
            ub1_cell{i} = units2(bins(i):bins(i+1),:);
            ub2_cell{i} = units3(bins(i):bins(i+1),:);
        else
            ua_cell{i} = units1(bins(i):bins(i+1)-1,:);
            ub1_cell{i} = units2(bins(i):bins(i+1)-1,:);
            ub2_cell{i} = units3(bins(i):bins(i+1)-1,:);
        end
    end

    uc1 = 0;
    uc2 = 0;

    for i = 1 :length(ua_cell)
        for j = 1: length(ub1_cell)
            uc1 = uc1 + sum(1./(pdist2(ua_cell{i},ub1_cell{j}).^6),'all');
            uc2 = uc2 + sum(1./(pdist2(ua_cell{i},ub2_cell{j}).^6),'all');
        end
    end
    E_vdw1 = -H*(diameter^6)*uc1/((pi^2) );
    E_sum1 = E_sum1 + E_vdw1;
    E_vdwss1(1,zzi) =E_vdw1;        
    E_vdw2 = -H*(diameter^6)*uc2/((pi^2) );
    E_sum2 = E_sum2 + E_vdw2;
    E_vdwss2(1,zzi) =E_vdw2;
    disp(['vdw',num2str(E_vdw1/kbT),num2str(E_vdw2/kbT)]);


%el
    gpu = gpuDevice;
    max_szie = gpu.AvailableMemory;

    max_number_double  =floor((max_szie/8)^0.5);

    max_chunk_size = floor(max_number_double*0.4);
    ligand_chunk_number = ceil(length(ligands1)/max_chunk_size);
    disp(['Ligands dividing points into ',num2str(ligand_chunk_number),' chunks...'])

    bins = 1 : max_chunk_size : length(ligands1);
    bins(end+1) = length(ligands1);
    la_cell = {};
    lb1_cell = {};
    lb2_cell = {};

    for i = 1 : ligand_chunk_number
        if i == ligand_chunk_number
            la_cell{i} = ligands1(bins(i):bins(i+1),:);
            lb1_cell{i} = ligands2(bins(i):bins(i+1),:);
            lb2_cell{i} = ligands3(bins(i):bins(i+1),:);
        else
            la_cell{i} = ligands1(bins(i):bins(i+1)-1,:);
            lb1_cell{i} = ligands2(bins(i):bins(i+1)-1,:);
            lb2_cell{i} = ligands3(bins(i):bins(i+1)-1,:);
        end
    end

    E_el1=zeros([1,length(kwater)]);
    E_el2=zeros([1,length(kwater)]);

    for i = 1 :length(la_cell)
        for j = 1: length(lb1_cell)
            Dl=pdist2(lb1_cell{j},la_cell{i});
            D2=pdist2(lb2_cell{j},la_cell{i});               
            examine_D1 = Dl<0.7;
            examine_D2 = D2<0.7;
            sum_examine_D1 = sum(examine_D1, 'all');
            sum_examine_D2 = sum(examine_D2, 'all');         
            if sum_examine_D1 == 0
            else
                if overlapp == 0
                    overlapp = ind;
                    disp(['Overlap:', num2str(overlapp)]);
                else
                    disp(['Overlap:', num2str(overlapp)]);
                end
            end

            for kw = 1: length(kwater)
                E_el1(kw) = E_el1(kw) + zb^2*lb*kbT*sum((1./Dl).*(exp(-kwater(kw)*Dl)),'all');
                E_el2(kw) = E_el2(kw) + zb^2*lb*kbT*sum((1./D2).*(exp(-kwater(kw)*D2)),'all');
            end
        end
    end

    for kw = 1: length(kwater)
        E_elss1(kw,zzi) =E_el1(kw);
        E_elss2(kw,zzi) =E_el2(kw);
        disp(['ionic',num2str(ionic_strength(kw)),'E_el ',num2str(E_el1(kw)/kbT),' ',num2str(E_el2(kw)/kbT)]);
        E_sum1_i=E_sum1+E_el1(kw);
        E_sum2_i=E_sum2+E_el2(kw);
        E_sumss1(kw,zzi)=E_sum1_i;
        E_sumss2(kw,zzi)=E_sum2_i;
        disp(['E_sum',num2str(E_sum1_i/kbT),' ',num2str(E_sum2_i/kbT)]);
        E_total = E_sum1_i+E_sum2_i*2;
        E_totalss(kw,zzi)=E_total;
    end

    zzi = zzi+1;
    toc
end
disp('Finished!');
disp([num2str(overlapp)]);

E_sumsskbt1=zeros([length(ionic_strength),length(index)]);
E_elsskbt1=zeros([length(ionic_strength),length(index)]);
E_sumsskbt2=zeros([length(ionic_strength),length(index)]);
E_elsskbt2=zeros([length(ionic_strength),length(index)]);
E_totalsskbt=zeros([length(ionic_strength),length(index)]);

%plot by ionic strength vdw
figure(5)
clf;
hold on
E_vdwsskbt1=double(gather(E_vdwss1/kbT));
E_vdwsskbt2=double(gather(E_vdwss2/kbT));

plot(index,E_vdwsskbt1,'r');
plot(index,E_vdwsskbt2,'b');
xlabel('Offset')
ylabel('$E_{vdW} (k_{B}T)$','interpreter','latex')
title('vdW')
legend('Upper', 'Lower')
ylim([-200 0])

drawnow

summin = zeros([1 length(ionic_strength)]);
for i = 1: length(ionic_strength)

    E_sumsskbt1(i,:)=double(gather(E_sumss1(i,:)/kbT));
    E_elsskbt1(i,:)=double(gather(E_elss1(i,:)/kbT));
    E_sumsskbt2(i,:)=double(gather(E_sumss2(i,:)/kbT));
    E_elsskbt2(i,:)=double(gather(E_elss2(i,:)/kbT));
    E_totalsskbt(i,:)=double(gather(E_totalss(i,:)/kbT));
    if i == 1
        figure(6)
        clf;
        hold on
        plot(index,E_vdwsskbt1,'--','Color','r');
        plot(index,E_elsskbt1(i,:),'--','Color','b');
        plot(index,E_sumsskbt1(i,:),'--','Color','g');
        plot(index,E_vdwsskbt2,':','Color','r');
        plot(index,E_elsskbt2(i,:),':','Color','b');
        plot(index,E_sumsskbt2(i,:),':','Color','g');
        plot(index,E_totalsskbt(i,:),'Color','k');
        xlabel('Offset')
        ylabel('$E_{sum} (k_{B}T)$','interpreter','latex')
        yline(0,'k')
        title('Sum')
        
        drawnow
    end
end
