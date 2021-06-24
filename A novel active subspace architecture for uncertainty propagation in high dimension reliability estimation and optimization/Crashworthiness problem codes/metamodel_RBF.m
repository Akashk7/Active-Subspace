function [srgt_RBF, PRESSRMS_RBF, eXV_RBF, srgtOPT_RBF, Y_hat_RBF, R2_pred] = metamodel_RBF(x,y)

X_RBF=x;
Y_RBF= y;
[srgt_RBF, PRESSRMS_RBF, eXV_RBF, srgtOPT_RBF] = build_RBF_SRGT(X_RBF,Y_RBF); 
Y_hat_RBF = srgtsRBFEvaluate(X_RBF, srgt_RBF);
%predvar_RBF = srgtsRBFPredictionVariance(X_RBF, srgt_RBF);

%  R^2 predicted
% sum of square about mean
   SST=sum((Y_hat_RBF-mean(Y_hat_RBF)).^2);
%  sum of square due to error
   SSE_pred=((sum(eXV_RBF.^2)));
   R2_pred=1-(SSE_pred/SST);
   