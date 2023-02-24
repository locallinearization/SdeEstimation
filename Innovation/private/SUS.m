function[Index]=SUS(SelCant,Sumvector) 

aux=rand(1,1);

toadd=1/SelCant;

cant=1;

while cant<=SelCant
 count=1;
  while aux>Sumvector(count) 
   count=count+1;
  end
  Index(cant)=count;
  cant=cant+1;
  aux=aux+toadd;
  if aux>1
   aux=aux-1;
  end
end

return

