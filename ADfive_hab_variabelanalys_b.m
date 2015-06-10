function ADsim(var1,var2)


% Use a matrix zz to represent lineages within the population
%A NN matrix represents the abunances for lineages (rows) in the habitats
%(collons)
%The KK matrix represents the K0 value of a lineage (traitvalue) in the
%different habitats

%This script loops over a vector with elements representing a parameter in
%the model, all other parameters are constant

%Set seeds for randfunctions
rand('seed',fix(sum(1000*clock)+var2))
randn('seed',fix(sum(2000*clock)+var2))



sigA_vector=[1.25 1 0.75 0.66 0.5 0.25 0.17 0.125 0.1]; % sigmaA values to loop over
Ropt_matrix=[-0.4 -0.2 0 0.2 0.4;
             -1.2 -0.6 0 0.6 1.2;
             -2 -1 0 1 2 ;
             -2.8 -1.4 0 1.4 2.8;
             -3.6 -1.8 0 1.8 3.6;
             -4.4 -2.2 0 2.2 4.4;
             -5.2 -2.6 0 2.6 5.2;
             -6 -3 0 3 6 ;
             -6.8 -3.4 0 3.4 6.8];

         

Ropt_values=Ropt_matrix(:,4);
%save parameter vector 
save sigmaA.mat sigA_vector 
save Ropt.mat Ropt_values

for iter=1:49
run=str2num([num2str(var1) num2str(var2) num2str(iter)]); %the name of this run

r = 1; % intrinsic growth rate
K0 = 1000; % Maximum K-value (carrying capacity)
nHab = 5; % Number of habitats

%Different resource dist. in differnt habitat
Ropt=Ropt_matrix(var1,:);

%Parameters of the model (comment when looped over vector above)
sigma_K = 1; % width of the resource base
sigma_a = sigA_vector(var2); % competition niche breadth
P_mut = 1e-3; % Mutation probability
sigma_mut = 0.02; % standard deviation of mutations
P_disp=1e-3; %Dispersal probability


%Probability for dispersal to one or the other habitat
% P_disp2 = [ 0 0.5 0;
%             1 0   1;
%             0 0.5 0];
P_disp2=[0 0.5 0 0 0; 
         1 0 0.5 0 0; 
         0 0.5 0 0.5 0;
         0 0 0.5 0 1;
         0 0 0 0.5 0];
P_disp_cum = cumsum(P_disp2);

%Seed the system
zz = 0; %The trait values in vector zz
NN = [0 0 10 0 0]; %Abundance in habitats in matrix NN
KK = K0*exp(-(zz-Ropt).^2/2/sigma_K^2); %K-values matching zz 

figure(run), clf

U=[]; %matrix for node information 

%For clustering function 
clus=zeros(1,length(zz))';
clusid=1;
lo=[0 length(zz)];

t=0; %counter for time
asym=[]; %for the asymptot evaluation

Uevo=[]; %matrix for saving of coll:runid,specid,time,meanZ,habitat/abund



while 1  %run the model until (x-timesteps) untill popsize reach asymptot
   
    % Do some plotting and evaluate popsize every 100th time step:
    if rem(t,100)==0
        doplot(t,zz,NN,KK,r,sigma_a,run);
        % Stop run if asymptote in popsize is reached
           asym(end+1)=sum(NN(:));
           if length(asym)>100
               %if t>40000
               if abs(mean(asym(end-50:end))-mean(asym(end-99:end-49)))<100
                   break
               end
               %end
           end
    end

    %New vectors ans matrices for the next generation
    newNN = zeros(size(NN));
    newzz = zz;
    newKK = KK;
    newclusind=clus;
    
    for i=1:length(zz)
        % For each lineage:
        %   Calculate fitness in each habitat (where it exists)
        fitness = zeros(1,nHab);
        for h = find(NN(i,:)>0)
            fitness(h) = threehabfitfunc(zz(i), KK(i,h), zz, NN(:,h), r, sigma_K, sigma_a); %into (zp, Kp, zz, N, r, sigma_K, sigma_a)
        end

        %   Calculate the number of offspring
        number_offspring = poissrnd(fitness.*NN(i,:));
        % Mutations:
        number_mutants = min(number_offspring,poissrnd(number_offspring*P_mut));
        for h = find(number_mutants>0)
            for m = 1:number_mutants(h)
                newzz(end+1,1) = zz(i) + sigma_mut*randn;
                newclusind(end+1,1)=clus(i);
                newNN(end+1,h) = 1;
                newKK(end+1,:) = K0*exp(-(newzz(end)-Ropt).^2/2/sigma_K^2);
            end
        end
        newNN(i,:) = number_offspring - number_mutants; 
    end
    NN = newNN;
    KK = newKK;
    zz = newzz;
    clus=newclusind;
    % remove extinct lineages:
    alive = sum(NN,2)>0;
    zz = zz(alive);
    NN = NN(alive,:);
    KK = KK(alive,:);
    clus=clus(alive);
    
    %   Dispersal
    for i=1:length(zz)
        for h = find(NN(i,:)>0)
            number_dispersers = binornd(NN(i,h), P_disp);
            for d = 1:number_dispersers
                new_habitat = find(rand < P_disp_cum(:,h), 1, 'first');
                NN(i,h) = NN(i,h)-1;
                NN(i,new_habitat) = NN(i,new_habitat)+1;
            end
        end
    end
    
    %Cluster phenotypic values into clusters(species)
       
    if length(zz)>1; %No use to start untill some variation exicist 
            new_u=[]; %species matrix for next generation 
        
        u=[[1:length(zz)]' zz clus]; %matrix with index, zz-values and cluster identity
        
        %Save a summary of u and NN every 100 generations
        if rem(t,100)==0
            Uevo_tmp=[];
            Uevo_tmp=[u NN]; %Matrix with coll:index, z-value, specname, abundance in hab 1-5
            %loop over unique speciesnames
            uniqspec=unique(Uevo_tmp(:,3));
            for i=uniqspec'
                specind=find(Uevo_tmp(:,3)==i);
                Uevo(end+1,1)=run;
                Uevo(end,2)=t;
                Uevo(end,3)=i;
                Uevo(end,4)=mean(Uevo_tmp(specind,2));
                
                if length(specind)==1
                    Uevo(end,5:9)=Uevo_tmp(specind,4:end);
                else
                    Uevo(end,5:9)=sum(Uevo_tmp(specind,4:end));
                end
                
                
            end
                %Uevo(end+1,:)=NaN;
        end
        for cluster= unique(u(:,3))'   %loop over the species to see if they has branched, one at a time 
            clear delta_u
            
            u_sort=u(find(u(:,3)==cluster),:); %pick the cluster to investigate
            u_sort = sortrows(u_sort,-2); %sort u_sort based on zz-values
            delta_u(length(zz)-1) = 0; 

            for i=1:length(u_sort(:,1))-1;
                delta_u(i)=abs(u_sort(i,2)-u_sort(i+1,2)); %differnce between zz from high to low
            end
            
            newlo=[];
            newlo= find(delta_u >0.1); %breakes in zz-values            

            %Warn the user if there are more then one break in a cluster
            if length(newlo)>1; 
                disp('WARNING more than one brake in cluster');
                disp(newlo);
                disp(u_sort);
            end
                       if length(newlo)>0

                        if newlo(1)>=2 & newlo(end)<=length(u_sort(:,1))-2 %only break when the new cluster is dimorf

                            %Matrix for information of branshing
                            for j=1:2
                                U(end+1,1)=u_sort(newlo(1),3);
                                U(end,2)=t;
                                U(end,3)=mean(u_sort(:,2));
                            end
                            %New cluster id after branching
                            u_sort(1:newlo(1),3)=clusid;
                            u_sort(newlo(1)+1:end,3)=clusid+1;
                            clusid=clusid+2; 

                        end
                       end


                new_u=[new_u; u_sort]; %New specis matrix with new clusternames
        end %end of cluster loop       
       
       %Find clusternames to put into next generation 
       u=sortrows(new_u,1);
       clus=(u(:,3));
    
    end
        t=t+1; 
       drawnow % This is a command which makes this program behave better.    
end


%find end community from U-matrix and save in Uend-matrix
[m n]=size(U);
Uend=[];
for i = m:-1:1
    indend=find(u(:,3)==i);
    if indend>0
        Uend(end+1,1)=i;
        Uend(end,2)=mean(u(indend,2));
    end
end

[m n]=size(Uend);
if m>1 %if only one species in community no use in doing tree and sample
    Uend=sortrows(Uend,-2);

    %Pair vise distance between species in end community
    [m n]=size(Uend);
    dist=zeros(m,m);

    for i=1:m
        node1=Uend(i,1);
        for j=1:m
            node2=Uend(j,1);
            if i~=j
                %line up node history for species
                while node1(end)>0
                    node1(end+1)=U(node1(end),1);
                end
                while node2(end)>0
                    node2(end+2)=U(node2(end),1);
                end
                %find first common node
                cnode=intersect(node1,node2);            
                cnodeind=find(U(:,1)==cnode(end));
                dist(i,j)= (t-U(cnodeind(1),2));
            end
        end
    end

      %U
      %Uend

    %Create phylogenetic tree and find the newick string

    %Extract the pairvise distances from the dist matrix
    dist2=dist; dist2(:,1)=[]; dist2(end,:)=[];
    PVdist=[];
    m=length(dist2);

    for i=1:m
        add=dist2(i,i:end);
        PVdist=[PVdist add];
    end

    %create the tree also ad names for leafs in the tree
    name={};
    namebase1 = 'species';
    for i=1:length(Uend(:,1))
        namebase2 = num2str(Uend(i,1));
        namebase=[namebase1 namebase2];  
        name(i)={num2str(namebase)}; 
    end
    
%    tree = seqlinkage(PVdist,'average',name);   %Use for large analysis (no saving)
    treename=['treedata' num2str(run) '.mat'];
    save(treename,'PVdist','name');
    %tree(var,iter) = seqlinkage(PVdist,'average',name); %use when you want to save tree in tree.mat
    %view(tree(var,iter))

%     newick=getnewickstr(tree); %input for phylocom (saved below)
%     % save treefile (newickformat)
%     fil=fopen(['newick' num2str(run)],'w');
%     fprintf (fil,newick);  
%     fclose(fil);


    %Create output file for input to phylocom
    
    phylocom_sample={};
    
    for i = 1:nHab
        nameHab=['habitat' num2str(i) num2str(run)];

        all_spec= unique(u(find(NN(:,i)>0),3));
        for ii=all_spec'
            sp_ind=find(u(:,3)==ii);
            abun= sum(NN(sp_ind,i));

            nameSpec=['species' num2str(ii)];
            
            phylocom_sample(end+1,1)={nameHab};
            phylocom_sample(end,2)={abun};
            phylocom_sample(end,3)={nameSpec};
            
        end
    end
else
    phylocom_sample(end+1,1)={'only_one'};
    phylocom_sample(end,2)={'species_in'};
    phylocom_sample(end,3)={'end_com'};
end
    
% 
%     %samplefile; three columns with habitat,abundance,species
%     %Extract info from model and write to file
%     fil=fopen(['phylosample' num2str(run)],'w');
%     for i = 1:nHab
%         nameHab=['habitat' num2str(i) num2str(run)];
% 
%         all_spec= unique(u(find(NN(:,i)>0),3));
%         for ii=all_spec'
%             sp_ind=find(u(:,3)==ii);
%             abun= sum(NN(sp_ind,i));
% 
%             nameSpec=['species' num2str(ii)];
%             fprintf(fil,'%s\t%i\t%s\t\n',nameHab,abun,nameSpec);
%         end
%     end
%     fclose(fil);
% else
%     disp('WARNING')
%     disp(['only one species in end community if run' num2str(run)])    
% end



%Save figure for this run
saveas(run,['figure' num2str(run)]);
%Save nodeinfo (U) and end community info (Uend) as MATLAB variables
filename = ['clusterinfo' num2str(run)];
save(filename, 'U', 'Uend','Uevo')
samplename=['phylocom_sample' num2str(run)]
save(samplename,'phylocom_sample')
close

end %end of iterationloop

%Save all trees from this session in one variable matrix 
%(row=same parameters, coll.=differnt parameters(the loop above))
        %save tree.mat tree         %Do not use for large analyses
 

   

%The fitness function 
function f = threehabfitfunc(zp, Kp, zz, N, r, sigma_K, sigma_a)

%Beräkna Neff, dvs det effektiva antalet konkurrenter
alphaij = exp(-(zp-zz).^2/2/sigma_a^2);
Neff = sum(alphaij.*N);

%Beräkna fitness
f = 1 + r*(1-Neff/Kp);
f = max(f,0);

%The Plotting function 
function doplot(t,zz,NN,KK,r,sigma_a,run)
set(0,'currentfigure',run) % Instead of "figure(1)"
subplot(221)
% Plot the position of each lineage along a timeline
sel = find(NN(:,1)>5);
plot(t+10*ones(length(sel),1),zz(sel),'.','markersize',10,'col','r')
hold on
sel = find(NN(:,2)>5);
plot(t+20*ones(length(sel),1),zz(sel),'.','markersize',10,'col','k')
sel = find(NN(:,3)>5);
plot(t+30*ones(length(sel),1),zz(sel),'.','markersize',10,'col','b')
sel = find(NN(:,4)>5);
plot(t+40*ones(length(sel),1),zz(sel),'.','markersize',10,'col','m')
sel = find(NN(:,5)>5);
plot(t+50*ones(length(sel),1),zz(sel),'.','markersize',10,'col','g')
% % Add a red dot for population mean:
% plot(t, mean(zz(2,:)), 'r.')
% ylabel('z')
% xlabel('time')

subplot(223)
% Plot total population size over time
plot(t,sum(NN(:)),'.'), hold on
ylabel('population size')
xlabel('time')

% subplot(222)
% % Plot the current fitness landscape:
% zplot = -15:0.1:15;
% fplot = zeros(size(zplot));
% for i=1:length(fplot)
%     fplot(i) = ADsimfitnessfunc(zplot(i),zz,K0,r,sigma_K,sigma_a);
% end
% plot(zplot,fplot);
% title('Fitness landscape')
% xlabel('z')
% ylabel('fitness function')

% Also plot the current distribution of z-values:
subplot(224)
zplot = -25:0.1:25;
bar(zz,sum(NN,2));
title('Genetic distribution')
xlabel('z')

%Plot the spatial distribution of individuals
subplot(222)
bar(1:size(NN,2),sum(NN,1))
title('Spatial distribution')
xlabel('habitat')

drawnow


