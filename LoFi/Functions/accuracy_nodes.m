function [Nodes,Elements] = accuracy_nodes(Accuracy,Nodes,Elements,Element_number)
%  [Nodes,Elements] = accuracy_nodes(accuracy,Nodes,Elements,Element_number)
%   accuracy_nodes creates extra nodes and elements in order to increase
%   the model accuracy
%
%   Nodes : [ Nodes ID , x , y , z] 
%   Elements : [ElID , Node 1 , Node 2 , D1 , D2 , t1 , t2 , E , rho , nu , Ca , Cd , marine growth , Type of element (1 = leg, 2 = brace , 3 = tower, 4 = pile sleeve), Material ID] 
%   Accuracy : The maximum lenght for each elements
%   Element_number : The element that is going to be divided


Element_line = find(Elements(:,1)==Element_number);
Node_line_1 = find(Nodes(:,1)==Elements(Element_line,2));
Node_line_2 = find(Nodes(:,1)==Elements(Element_line,3));
ctp_elem = size(Elements,1)+1;
ctp_node=size(Nodes,1)+1;
Node_accu(1,:) = Nodes(Node_line_1,:);
Node_accu(1,1) = 1;
Node_accu(2,:) = Nodes(Node_line_2,:);
Node_accu(2,1) = 2;
L_tot= norm(Node_accu(1,2:4)-Node_accu(2,2:4)); 
variation_diameter = (Elements(Element_line,5) - Elements(Element_line,4))/L_tot;
origin_diameter = Elements(Element_line,4);
variation_thickness = (Elements(Element_line,7) - Elements(Element_line,6))/L_tot;
origin_thickness = Elements(Element_line,6);
     j=1;
      while j<size(Node_accu,1)
        while norm(Node_accu(j,2:4)-Node_accu(j+1,2:4))>Accuracy 
            for k=j+1:size(Node_accu,1)
              Node_inter(k+1,:)= [ k+1, Node_accu(k,2), Node_accu(k,3), Node_accu(k,4)];
          end
              Node_inter(j+1,:) = [j+1, (Node_accu(j,2)+Node_accu(j+1,2))/2 , (Node_accu(j,3)+Node_accu(j+1,3))/2 , (Node_accu(j,4)+Node_accu(j+1,4))/2];
          cpt = j+1;
          for k=j+1:size(Node_inter,1)
               Node_accu(cpt,:) = Node_inter(k,:);
               cpt=cpt+1;
           end
        end
        j=j+1;
      end
    
for i=1:size(Node_accu,1)-1  
    L_elem(i) = norm(Node_accu(i,2:4)-Node_accu(i+1,2:4));
end
    %The new elements and nodes are created 
x=0;
if size(Node_accu,1)>2
    ctp_accu = 1;
    Elements_accu(ctp_accu,1:10) = [ ctp_elem , Elements(Element_line,2) , ctp_node , origin_diameter , variation_diameter*(x+L_elem(1)) + origin_diameter , origin_thickness , variation_thickness*(x+L_elem(1)) + origin_thickness , Elements(Element_line,8) , Elements(Element_line,9) , Elements(Element_line,10)];
    Elements_accu(ctp_accu,15) = Elements(Element_line,15);
    Elements_accu(ctp_accu,14) = Elements(Element_line,14);
    ctp_elem = ctp_elem + 1;
    ctp_accu = ctp_accu + 1;
    x=x+L_elem(1);
    for i =2:size(Node_accu,1)-2
        Elements_accu(ctp_accu,1:10) = [ ctp_elem , ctp_node , ctp_node+1 , variation_diameter*x + origin_diameter , variation_diameter*(x+L_elem(i)) + origin_diameter , variation_thickness*x + origin_thickness , variation_thickness*(x+L_elem(i)) + origin_thickness , Elements(Element_line,8) , Elements(Element_line,9) , Elements(Element_line,10)];
        Elements_accu(ctp_accu,15) = Elements(Element_line,15);
        Elements_accu(ctp_accu,14) = Elements(Element_line,14);
        Nodes(ctp_node,:) = Node_accu(i,:);
        Nodes(ctp_node,1) = ctp_node;
        ctp_node = ctp_node + 1; 
        ctp_elem = ctp_elem + 1;
        ctp_accu = ctp_accu + 1;
        x=x+L_elem(i);
    end
        Nodes(ctp_node,:) = Node_accu(size(Node_accu,1)-1,:);
        Nodes(ctp_node,1) = ctp_node;

        Elements_accu(ctp_accu,1:10) = [ ctp_elem , ctp_node , Elements(Element_line,3) , variation_diameter*x + origin_diameter , variation_diameter*L_tot + origin_diameter , variation_thickness*x + origin_thickness , variation_thickness*L_tot + origin_thickness , Elements(Element_line,8) , Elements(Element_line,9) , Elements(Element_line,10)];
        Elements_accu(ctp_accu,15) = Elements(Element_line,15);
        Elements_accu(ctp_accu,14) = Elements(Element_line,14);
        
    if Element_line==1
        Elements_inter (1:size(Elements_accu,1),:) = Elements_accu(:,:); 
        Elements_inter (size(Elements_accu,1)+1:size(Elements_accu,1)+size(Elements,1)-1,:) = Elements(2:end,:); 
    elseif  Element_line==size(Elements,1)    
        Elements_inter (1:size(Elements,1)-1,:) = Elements (1:size(Elements,1)-1,:); 
        Elements_inter (Element_line:size(Elements_accu,1) + size(Elements,1)-1,:) = Elements_accu; 
    else  
        Elements_inter (1:Element_line-1,:) = Elements (1:Element_line-1,:); 
        Elements_inter (Element_line:Element_line+size(Elements_accu,1)-1,:) = Elements_accu(:,:); 
        Elements_inter (Element_line+size(Elements_accu,1):size(Elements_accu,1) + size(Elements,1)-1,:) = Elements (Element_line+1:end,:); 
    end
    Elements = Elements_inter;
else
    
    Nodes = Nodes;
    Elements = Elements;
end
end

