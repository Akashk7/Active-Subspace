function [srgt_KRG, PRESSRMS_KRG, eXV_KRG, srgtOPT_KRG, Y_hat_KRG, predvar_KRG, R2_pred] = metamodel_KRG(x,y)

X_KRG=x;
Y_KRG= y;
[srgt_KRG, PRESSRMS_KRG, eXV_KRG, srgtOPT_KRG] = build_KRG_SRGT(X_KRG,Y_KRG); 
Y_hat_KRG = srgtsKRGEvaluate(X_KRG, srgt_KRG);
predvar_KRG = srgtsKRGPredictionVariance(X_KRG, srgt_KRG);

%  R^2 predicted
% sum of square about mean
   SST=sum((Y_hat_KRG-mean(Y_hat_KRG)).^2);
%  sum of square due to error
   SSE_pred=((sum(eXV_KRG.^2)));
   R2_pred=1-(SSE_pred/SST);
   