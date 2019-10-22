function In=normarea2(I,A)
%%%% In:  I = Binary Image, A = area of normalization (px)     
%%%% Out:  In = Normalized Image 
Area=regionprops(I,'Area');
Area=Area.Area;
k=sqrt(A/Area);
In=imresize(I,k);


    
    

