function [matches, P_init] = posMatching7(tar_feat,ref_feat,Indexs, Scores, match_id,C,th1)

Xall = tar_feat(1:2,:)';
Yall = ref_feat(1:2,:)';
[Xall, Yall] = normr(Xall, Yall);
Xall = [[1:size(Xall,1)]',Xall];
Yall = [[1:size(Yall,1)]',Yall];

Xi = Xall(match_id(1,:),:);

Xo = Xall; Xo(match_id(1,:),:) = [];

Kn1 = 10;
kdtreeXi = vl_kdtreebuild(Xi(:,2:3)');
[~, distanceX] = vl_kdtreequery(kdtreeXi, Xi(:,2:3)', Xo(:,2:3)','NumNeighbors', Kn1);
K_dist1 = distanceX(Kn1,:);
epsilon = 0.005;
X_index = find(K_dist1<epsilon);

Xoi = Xo(X_index,:);
M = 15;
[K, ~] = FourierKernel(Xoi(:,2:3), M);
Xta = K*C;

match_pos = []; k = size(Scores,2);
for i = 1:length(X_index)
    in_xoi = Xoi(i,1);
    rescore = Scores(in_xoi,:)./[Scores(in_xoi,1),Scores(in_xoi,1:k-1)];
    ii = find(rescore>th1)-1;
    if ~isempty(ii)
        radius = (rescore(ii(1)+1)-1);
        index_i = Indexs(in_xoi,1:ii(1));
        desNN = setdiff(index_i, match_id(2,:));
        if ~isempty(desNN)
            distdiff = sqrt(sum((repmat(Xta(i,:),length(desNN),1)- Yall(desNN,2:3)).^2,2))';
            matchpos_y = desNN(distdiff<radius);
            matchp = [in_xoi*ones(1,length(matchpos_y)); matchpos_y];
            match_pos = [match_pos,matchp];
        end
    end
end

matches = [match_id, match_pos];
P_init = [ ones(size(match_id,2),1);0.00001*ones(size(match_pos,2),1)];
end