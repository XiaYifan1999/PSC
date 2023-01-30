function [ inliers ] = ground_truth_verification( feat1, feat2, matches, ground_truth, threshold )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ground_truth = single(ground_truth);
db_point = feat1(1:2,:);
db_point(3,:) = 1;

warped_points = ground_truth * db_point;
warped_points(1,:) = warped_points(1,:)./ warped_points(3,:);
warped_points(2,:) = warped_points(2,:)./ warped_points(3,:);
warped_points(3,:) = warped_points(3,:)./ warped_points(3,:);

query_point = feat2(1:2,:);
query_point(3,:) = 1;


X = warped_points(:,matches(1,:));
Y = query_point(:,matches(2,:));
dist = X - Y;
dist = dist(1,:) .* dist(1,:) + dist(2,:) .* dist(2,:);
dist = sqrt(dist);

inliers = dist <= threshold;

end

