function matches = PSC(tar_feat, ref_feat, tar_desc, ref_desc)

% Fourier-base Progressive Matching - version.2
% by Yifan Xia, 15.4.2022
k =10;  th0 =1.5; th1 = 1.1;
[matches, Indexs, Scores] = desdist_match(tar_desc, ref_desc, k, th0);

X = tar_feat(1:2,matches(1,:))';
Y = ref_feat(1:2,matches(2,:))';

[Xn, Yn] = normr(X, Y);
conf = CRC_LLT_init([]);

Xm = Xn; Ym = Yn;
iters = 0; n_matid = 0;
P_init=[];C =[];

while true
    [idx, ~, C] = CRC_LLT(Xm, Ym, conf, P_init,C);
    %     [idx, ~, C] = CRC(Xm, Ym, conf);
    match_id = matches(:,idx); 
    if   ~(size(match_id,2)>size(n_matid,2)) % compact representation is not increasing
        matches = match_id; % fprintf('stop iterations %.0f, no extension.\n',iters)
        break;     
    end

    n_matid = match_id;
    [matches, P_init] = posMatching7(tar_feat,ref_feat,Indexs, Scores,match_id,C,th1);
    
    if  ~(size(matches,2)>size(match_id,2)) % posMatching is unusable
        % fprintf('stop iterations %.0f, no pos-matching.\n',iters+1)
        break; 
    end
    
    iters = iters+1;
    Xm = tar_feat(1:2,matches(1,:))';Ym = ref_feat(1:2,matches(2,:))';
    [Xm, Ym] = normr(Xm, Ym);
end 

