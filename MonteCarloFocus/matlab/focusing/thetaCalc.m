function [theta] = thetaCalc(y,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


theta=atan(y/(2*f-sqrt(f^2-y^2)));
end

