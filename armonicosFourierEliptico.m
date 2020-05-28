
% % Create hexagon polyshape
% pgon = nsidedpoly(6,'Center',[100,100],'Radius',50);
% % Generate 100x100 binary image with hexagon
% [xGrid,yGrid] = meshgrid(1:200,1:200);
% BW = isinterior(pgon,xGrid(:),yGrid(:));
% BW = reshape(BW,size(xGrid));
function cons = armonicosFourierEliptico(armon,img_route)

    N = armon; %Numero de armonicos
    im=imread(img_route);

    %im = imread('im_pruebas/a1.jpg');
    %im = imrotate(im,0,'bicubic');
    %imshow(im);
    BW = imbinarize(im(:,:,1));
    [B,L] = bwboundaries(BW);
    %figure,imshow(label2rgb(L, @jet, [.5 .5 .5]))
    %imshow(a_bw);
    for k = 1:size(B,1)
       boundary = B{k};
       %figure,plot(boundary(:,2),boundary(:,1))
    end

    initial_point = boundary(1,:);
    shape = size(BW);
    actual_point = initial_point;
    length = size(boundary,1);
    chain_code = "";
    chain_code_array = zeros(size(boundary,1),1);
    count = 1;

    while count <= size(boundary,1)
        ac_bin_i = and(boundary(:,1)==actual_point(1,1),boundary(:,2)==actual_point(1,2));
        actual_indexes = find(ac_bin_i==1);
        rel_indexes = find_relational_indexes(actual_indexes,length);
        rel_points = boundary(rel_indexes,:);
        rel_codes = obtain_codes(rel_points,actual_point);
        [rel_codes, pos] = unique(rel_codes);
        rel_points = rel_points(pos,:);

        if(actual_point == initial_point)
            %if(all(ismember(rel_codes,[0,1,2,6,7]))) %Cuadrante derecho
            [code,pos]=max(rel_codes);
            actual_point = rel_points(pos,:);
            chain_code = strcat(chain_code,int2str(code));
            chain_code_array(count) = code;
            %end
        else
%             if(count == 117)
%                chain_code 
%             end
            next_code = obtain_codes_no_proc(rel_codes,chain_code_array(count-1));
            actual_point = rel_points(rel_codes == next_code(1),:);
            chain_code = strcat(chain_code,int2str(next_code(1)));
            chain_code_array(count) = next_code(1);
        end

        count = count + 1;
    end



    T = get_tp(chain_code_array,size(chain_code_array,1));
    t = 1:round(T);

    % t = zeros(size(chain_code_array,1),1);
    % for p = 1:size(chain_code_array,1)
    %     t(p) = get_tp(chain_code_array,p);
    % end

    [yN,xN,excen,cons] = get_both_xN_yN(chain_code_array,t,N);
    %xN = get_function_xN(chain_code_array,t,N);
    %yN = get_function_yN(chain_code_array,t,N);

    % tp = get_tp(chain_code_array,1);
    % get_xp(chain_code_array,272);
    % get_yp(chain_code_array,3)

   % rec_im = reconstructe_im(shape,chain_code_array,initial_point);
    %figure,imshow(rec_im);
    %figure,plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
    %hold on
    %figure,plot(xN,yN);
    %cons = cons';
    cons = cons(:);
end

%%%%%%Functions
function rel_indexes = find_relational_indexes(actual_indexes,length)
    if(actual_indexes == 1 | actual_indexes == length)
        rel_indexes = actual_indexes;
        rel_indexes(actual_indexes == 1) = 2;
        rel_indexes(actual_indexes == length) = length - 1;
    else
        rel_indexes = [actual_indexes-1;actual_indexes+1];
    end
end

function codes = obtain_codes(rel_points,actual_point)
    codes = zeros(size(rel_points,1),1);
    codes_cases = [0,1;1,1;1,0;1,-1;0,-1;-1,-1;-1,0;-1,1];
    for i = 1:size(rel_points,1)
        temp_exp = rel_points(i,:)-actual_point;
        if temp_exp == codes_cases(1,:)
            codes(i) = 0;
        elseif temp_exp == codes_cases(2,:)
            codes(i) = 1;
        elseif temp_exp == codes_cases(3,:)
            codes(i) = 2;
        elseif temp_exp == codes_cases(4,:)
            codes(i) = 3;
        elseif temp_exp == codes_cases(5,:)
            codes(i) = 4;
        elseif temp_exp == codes_cases(6,:)
            codes(i) = 5;
        elseif temp_exp == codes_cases(7,:)
            codes(i) = 6;
        elseif temp_exp == codes_cases(8,:)
            codes(i) = 7;
        end
    end
end

function codes_wop = obtain_codes_no_proc (rel_codes , past_code)
    codes_procedence = [4,5,6,7,0,1,2,3];
    codes_wop = rel_codes;
    if(size(codes_wop,1)>1)
        for i = 1:size(rel_codes,1)
            if(codes_procedence(past_code + 1)==rel_codes(i))
                codes_wop(i) = [];
            end
        end
        if(codes_procedence(past_code + 1) == 4)
            codes_wop = codes_wop(size(codes_wop,1));
        end
    end
end

function rec_img = reconstructe_im (shape,chain_code_array,initial_point)
    rec_img = zeros(shape);
    actual_point = initial_point;
    rec_img(actual_point(1),actual_point(2)) = 1;
    for j = 1:size(chain_code_array,1)
       actual_point(1) = actual_point(1) + get_delta_y(chain_code_array(j));
       actual_point(2) = actual_point(2) + get_delta_x(chain_code_array(j));
       rec_img(actual_point(1),actual_point(2)) = 1;
    end
end

function delta_ti = get_delta_ti(chain_code_array) 
    delta_ti = 1+((sqrt(2)-1)/2)*(1-(-1).^chain_code_array);    
end

function tp = get_tp(chain_code_array,p)
    if(p~=0)
        tp = sum(get_delta_ti(chain_code_array(1:p)));
    else
        tp = 0;
    end
end

function delta_x = get_delta_x(chain_code) 
    delta_x = sgn_Z(6 - chain_code) * sgn_Z(2 - chain_code);
end

function delta_y = get_delta_y(chain_code) 
    delta_y = sgn_Z(4 - chain_code) * sgn_Z(chain_code);
end

function sgn = sgn_Z(n)
    if(n>0)
        sgn = 1;
    elseif(n<0)
        sgn = -1;
    else
        sgn = 0;
    end
end

function xp = get_xp(chain_code_array,p) 
    xp = 0;
    for i = 1:p
        xp = xp + get_delta_x(chain_code_array(i));
    end
end

function yp = get_yp(chain_code_array,p) 
    yp = 0;
    for i = 1:p
        yp = yp + get_delta_y(chain_code_array(i));
    end
end

function a_n = get_a_n (T,chain_code_array,n)
    a_n = 0;
    
    for p = 1:size(chain_code_array,1)
        a_n = a_n + get_delta_x(chain_code_array(p))/get_delta_ti(chain_code_array(p)) *(cos(2*pi*n*get_tp(chain_code_array,p)/T)-cos(2*pi*n*get_tp(chain_code_array,p-1)/T));
    end
    
    a_n = a_n * (T/(2*n^2*pi^2));
end

function b_n = get_b_n (T,chain_code_array,n)
    b_n = 0;
    
    for p = 1:size(chain_code_array,1)
        b_n = b_n + get_delta_x(chain_code_array(p))/get_delta_ti(chain_code_array(p)) *(sin(2*pi*n*get_tp(chain_code_array,p)/T)-sin(2*pi*n*get_tp(chain_code_array,p-1)/T));
    end
    
    b_n = b_n * (T/(2*n^2*pi^2));
end

function c_n = get_c_n (T,chain_code_array,n)
    c_n = 0;
    
    for p = 1:size(chain_code_array,1)
        c_n = c_n + get_delta_y(chain_code_array(p))/get_delta_ti(chain_code_array(p)) *(cos(2*pi*n*get_tp(chain_code_array,p)/T)-cos(2*pi*n*get_tp(chain_code_array,p-1)/T));
    end
    
    c_n = c_n * (T/(2*n^2*pi^2));
end

function d_n = get_d_n (T,chain_code_array,n)
    d_n = 0;
    
    for p = 1:size(chain_code_array,1)
        d_n = d_n + get_delta_y(chain_code_array(p))/get_delta_ti(chain_code_array(p)) *(sin(2*pi*n*get_tp(chain_code_array,p)/T)-sin(2*pi*n*get_tp(chain_code_array,p-1)/T));
    end
    
    d_n = d_n * (T/(2*n^2*pi^2));
end

function A0 = get_A0(T,chain_code_array)
    A0 = 0;
    for p = 1:size(chain_code_array,1)
       A0 = A0 + ((get_delta_x(chain_code_array(p))/(2*get_delta_ti(chain_code_array(p)))) ...
        * (get_tp(chain_code_array,p)^2-get_tp(chain_code_array,p-1)^2)) + get_epsilon_p(p,chain_code_array)*...
         (get_tp(chain_code_array,p)-get_tp(chain_code_array,p-1));
    end
    
    A0 = 1/T * A0;
end

function C0 = get_C0(T,chain_code_array)
    C0 = 0;
    for p = 1:size(chain_code_array,1)
       C0 = C0 + ((get_delta_y(chain_code_array(p))/(2*get_delta_ti(chain_code_array(p)))) ...
        * (get_tp(chain_code_array,p)^2-get_tp(chain_code_array,p-1)^2)) + get_Delta_p(p,chain_code_array)*...
         (get_tp(chain_code_array,p)-get_tp(chain_code_array,p-1));
    end
    
    C0 = 1/T * C0;
end

function ep_p = get_epsilon_p (p,chain_code_array)
    if(p>1)
        sum_x = get_xp(chain_code_array,p-1);
        ep_p = sum_x - get_tp(chain_code_array,p-1)*(get_delta_x(chain_code_array(p))/get_delta_ti(chain_code_array(p))); 
    else
        ep_p = 0;
    end
end

function DEL_p = get_Delta_p (p,chain_code_array)
    if(p>1)
        sum_y = get_yp(chain_code_array,p-1);
        DEL_p = sum_y - (get_delta_y(chain_code_array(p))/get_delta_ti(chain_code_array(p)))*get_tp(chain_code_array,p-1); 
    else
        DEL_p = 0;
    end
end

function xN = get_function_xN(chain_code_array,t_array,N)
    T = get_tp(chain_code_array,size(chain_code_array,1)); 
    A0 = get_A0(T,chain_code_array);
    xN = zeros(size(t_array));
    for n = 1:N
       a_n = get_a_n(T,chain_code_array,n);
       b_n = get_b_n(T,chain_code_array,n);
       xN = xN + a_n * (cos((2*n*pi).*t_array/T)) + b_n * (sin((2*n*pi).*t_array/T));
    end
    xN = A0 + xN;
end

function yN = get_function_yN(chain_code_array,t_array,N)
    T = get_tp(chain_code_array,size(chain_code_array,1)); 
    C0 = get_C0(T,chain_code_array);
    yN = zeros(size(t_array));
    
    for n = 1:N
       c_n = get_c_n(T,chain_code_array,n);
       d_n = get_d_n(T,chain_code_array,n);
       yN = yN + c_n * (cos(2*n*pi.*t_array/T)) + d_n * (sin(2*n*pi.*t_array/T));
    end
    
    yN = C0 + yN;
end

function [yN,xN,excen,cons] = get_both_xN_yN(chain_code_array,t_array,N)
    T = get_tp(chain_code_array,size(chain_code_array,1)); 
    C0 = get_C0(T,chain_code_array);
    A0 = get_A0(T,chain_code_array);
    excen = zeros(N,1);
    xN = zeros(size(t_array));
    yN = zeros(size(t_array));
    cons = zeros(N,4);
    
    for n = 1:N
       a = get_a_n(T,chain_code_array,n);
       b = get_b_n(T,chain_code_array,n);

       c = get_c_n(T,chain_code_array,n);
       d = get_d_n(T,chain_code_array,n);
       
        if( n == 1)
           theta1 = 1/2 * atan2((2*(a*b+c*d)),(a^2+c^2-b^2-d^2));
           
           if(theta1<0)
              theta1 = theta1 + pi;
              
           end
           
           ast_val = [cos(theta1),sin(theta1);-sin(theta1),cos(theta1)]*[a,c;b,d];
           phi1 = atan2(ast_val(1,2),ast_val(1,1));
          % if(phi1<0)
           %   phi1 = phi1 + 2*pi; 
           %end
           E = sqrt(ast_val(1,2)^2+ast_val(1,1)^2); %Magnitud del semieje mayor de la primera elipse
        end
       theta2 = 1/2 * atan2((2*(a*b+c*d)),(a^2+c^2-b^2-d^2));
           
       if(theta2<0)
          theta2 = theta2 + pi;

       end
           
       ast_val = [cos(theta2),sin(theta2);-sin(theta2),cos(theta2)]*[a,c;b,d];
       major = sqrt(ast_val(1,2)^2+ast_val(1,1)^2); %Magnitud del semieje mayor de la elipse actual
       minor = sqrt(ast_val(2,2)^2+ast_val(2,1)^2); %Magnitud del semieje menor de la elipse actual
       excen(n) = sqrt(major^2-minor^2)/major; %Excentricidad de la elipse actual
       
       res_1=[cos(phi1),sin(phi1);-sin(phi1),cos(phi1)]*[a,b;c,d]*[cos(n*theta1),-sin(n*theta1);sin(n*theta1),cos(n*theta1)];
       %res_2 = res_1;
       res_2 = res_1/E;
       cons(n,:) = res_2(:);
       xN = xN + res_2(1,1) * (cos((2*n*pi).*t_array/T)) + res_2(1,2) * (sin((2*n*pi).*t_array/T));
       yN = yN + res_2(2,1) * (cos(2*n*pi.*t_array/T)) + res_2(2,2) * (sin(2*n*pi.*t_array/T));
       %figure,plot(res_2(1,1) * (cos((2*n*pi).*t_array/T)) + res_2(1,2) * (sin((2*n*pi).*t_array/T)),res_2(2,1) * (cos(2*n*pi.*t_array/T)) + res_2(2,2) * (sin(2*n*pi.*t_array/T)))
    end
end

