function [tar_feat,tar_desc, ref_feat, ref_desc] = sift_process(img1,img2)

if size(img1,3)>1
    img1 = rgb2gray(img1);
end
if size(img2,3)>1
    img2 = rgb2gray(img2);
end

IMG1 = im2single(img1); IMG2 = im2single(img2);
[tar_feat, tar_desc] = vl_sift(IMG1);
[ref_feat, ref_desc] = vl_sift(IMG2);

end
