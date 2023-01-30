function [M0, Indexs, Scores] = desdist_match(tar_desc, ref_desc, k, ratio0)

D_o = vl_alldist(single(tar_desc),single(ref_desc));

[Scores, Indexs] = mink(D_o, k, 2);

indd = (Scores(:,2)./Scores(:,1))>ratio0;
xind = find(indd);
yind = Indexs(xind,1);
M0 = [xind'; yind'];  

end