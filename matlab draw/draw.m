function draw

     data=importdata('E:\ADI_FDTD\ADI_FDTD_steam\result\platform_matel.txt');
     
     Z=data(:,1)+1;
  
     X=data(:,2)+1;
     C=data(:,3);
     Cmn=0.8;
     maxz=max(Z);
     a=0.02286;
     z=0.1;
     b=0.01016;
     Y=0;  dy=0.02;
        dz=z/maxz;
       
     maxx=max(X);
        dx=a/maxx; 
        
        
        m=1;n=0;
     img=zeros(maxz,maxx);
     
     for i=1:length(Z)
       f=[Z(i),X(i)]
       img(Z(i),X(i))=  C(i);
     end
     %imshow(img')
     subplot(1,2,1)
     surfc(img) 
     subplot(1,2,2)
     contourf(img)
     hold on

     
    
     
    
end


