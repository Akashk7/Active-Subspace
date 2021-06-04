function [R2,RMSE]=metrics(y_actual,y_predicted)
y_mean = mean(y_actual);
n=size(y_actual,1);
% R^2 prediction
% sum of square due to error
  SSE=(sum((y_actual-y_predicted).^2));
% sum of square about mean
  SST=sum((y_actual-y_mean).^2);
  R2=1-(SSE/SST);
   
%  RMSE
RMSE = ((sum((y_actual-y_predicted).^2))/n).^0.5;