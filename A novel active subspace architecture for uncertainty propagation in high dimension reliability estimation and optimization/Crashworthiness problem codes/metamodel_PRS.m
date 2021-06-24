function [srgt_PRS, PRESSRMS_PRS, eXV_PRS, srgtOPT_PRS, Y_hat_PRS, predvar_PRS, R2_pred] = metamodel_PRS(x,y)

X_PRS=x;
Y_PRS= y;
[srgt_PRS, PRESSRMS_PRS, eXV_PRS, srgtOPT_PRS] = build_PRS_SRGT(X_PRS,Y_PRS); 
Y_hat_PRS = srgtsPRSEvaluate(X_PRS, srgt_PRS);
predvar_PRS = srgtsPRSPredictionVariance(X_PRS, X_PRS, srgt_PRS);

%  R^2 predicted
% sum of square about mean
   SST=sum((Y_hat_PRS-mean(Y_hat_PRS)).^2);
%  sum of square due to error
   SSE_pred=((sum(eXV_PRS.^2)));
   R2_pred=1-(SSE_pred/SST);
   