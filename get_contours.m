function [cont,n_cont] = get_contours(theta,phi,mars_topo,elevation)
% author: Marc Hesse
% date: Dec 8, 2021
% description: Extracts all contour for a given elevation and sorts them by
% length
    C = contourc(theta,phi,mars_topo,elevation*[1 1]); %hold on
%     close
    [row,col] = size(C);
    k=1; i = 1;

    while k <= col
        len = C(2,k);
        x = C(1,k+1:k+len);
        y = C(2,k+1:k+len);
        cont(i).x = x';
        cont(i).y = y';
        cont(i).N = len;
        Nvec(i) = len;
        k = k+1+len;
        i = i+1;
    end
    % Sort contours by length
    [Nsort,i_sort] = sort(Nvec,'descend');
    cont = cont(i_sort);
    n_cont = length(cont);
end